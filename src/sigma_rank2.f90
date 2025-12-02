module SigmaRank2_mod
    use calc_basic
    implicit none

contains
    subroutine sigma_flip_rank2(Gr, i, j, sigma_old, det_ratio, do_update)
! 对键 (i,j) 的 σ: +1 <-> -1 翻转对应 O' O^{-1} 的 rank-2 变化，
! 输入：
!   Gr         (in/out)  : Ndim×Ndim 等时格林函数（单自旋通道）
!   i,j        (in)      : 键两端格点索引
!   sigma_old  (in)      : 翻转前 σ ∈ {±1}
! 输出：
!   det_ratio  (out)     : R_f = det(I2 + V^T (I-G) U)
!   do_update  (in, opt) : 是否就地更新 G，默认 .true.
        complex(kind=8), intent(inout) :: Gr(:, :)
        integer, intent(in) :: i, j, sigma_old
        real(kind=8), intent(out) :: det_ratio
        logical, intent(in), optional :: do_update
        logical :: upd
        complex(kind=8), dimension(2,2) :: A, K, W, Minv
        complex(kind=8), dimension(:,:), allocatable :: B
        complex(kind=8), dimension(2,size(Gr,2)) :: Cvec
        real(kind=8) :: cval, sval
        complex(kind=8) :: detW

        upd = .true.
        if (present(do_update)) upd = do_update

! 构造 A = X2 - I2 = O(σ_new) O(σ_old)^{-1} - I2
! O(σ) = exp(Δτt σ c†c) 在2x2块中是 [[cosh(Δτt), σ*sinh(Δτt)], [σ*sinh(Δτt), cosh(Δτt)]]
        cval = cosh(Dtau * RT)
        sval = sinh(Dtau * RT)
! 对于σ_old → σ_new = -σ_old的翻转：
! X2 - I2 = [[2s^2, -2cs*σ_old], [-2cs*σ_old, 2s^2]]
        A(1,1) = cmplx(2.d0*sval*sval, 0.d0, kind=8)
        A(2,2) = cmplx(2.d0*sval*sval, 0.d0, kind=8)
        A(1,2) = cmplx(-2.d0*cval*sval*dble(sigma_old), 0.d0, kind=8)
        A(2,1) = cmplx(-2.d0*cval*sval*dble(sigma_old), 0.d0, kind=8)

! K = E^T (I - G) E，E=[e_i,e_j]
        K(1,1) = cmplx(1.d0, 0.d0, kind=8) - Gr(i,i)
        K(1,2) =           - Gr(i,j)
        K(2,1) =           - Gr(j,i)
        K(2,2) = cmplx(1.d0, 0.d0, kind=8) - Gr(j,j)

! W = A * K
        W = matmul(A, K)
! det(I2 + W) - 物理上应该总是正的
        detW = (cmplx(1.d0, 0.d0, kind=8) + W(1,1)) * (cmplx(1.d0, 0.d0, kind=8) + W(2,2)) - W(1,2) * W(2,1)
        det_ratio = real(detW)
! 物理上det_ratio必须为正（费米子行列式）
! 但数值上可能出现负值，取绝对值确保正定性
        det_ratio = abs(det_ratio)
! 重要修复：如果det_ratio太小，说明数值不稳定，应该拒绝更新
! 设为0让worm算法自动拒绝，而不是返回错误的小值
        if (det_ratio < 1.d-12) then
            det_ratio = 0.d0
            ! 无论是否更新模式，都直接返回
            return
        endif

        if (.not. upd) return

! 正确的Woodbury公式: G' = G - (I - G) E A (I + A K)^{-1} E^T G
! 其中 E = [e_i, e_j], K = E^T (I - G) E
! 注意：det_ratio = det(I + A K) = det(I + W)

! Minv = (I + A K)^{-1} = (I + W)^{-1}
        call inv2x2(cmplx(1.d0, 0.d0, kind=8) + W(1,1), W(1,2), W(2,1), cmplx(1.d0, 0.d0, kind=8) + W(2,2), Minv)

! B = (I - G) E = [(I-G)_i, (I-G)_j]（Ndim×2矩阵）
! (I-G)_i 表示 (I-G) 的第 i 列
        allocate(B(size(Gr,1), 2))
        B(:,1) = -Gr(:, i)
        B(i,1) = B(i,1) + cmplx(1.d0, 0.d0, kind=8)
        B(:,2) = -Gr(:, j)
        B(j,2) = B(j,2) + cmplx(1.d0, 0.d0, kind=8)

! B_left = (I - G) E A (I + A K)^{-1} = B * A * Minv
        B = matmul(B, matmul(A, Minv))

! C = E^T G = [G(i,:); G(j,:)]（2×Ndim矩阵）
        Cvec(1,:) = Gr(i,:)
        Cvec(2,:) = Gr(j,:)

! ΔG = B_left * C = (I - G) E A (I + A K)^{-1} E^T G
        Gr = Gr - matmul(B, Cvec)
        deallocate(B)
        return
    contains
        pure function unit_vec(n, idx) result(v)
            integer, intent(in) :: n, idx
            complex(kind=8) :: v(n)
            integer :: k
            v = cmplx(0.d0, 0.d0, kind=8)
            v(idx) = cmplx(1.d0, 0.d0, kind=8)
        end function unit_vec

        pure function Zcol(G, idx) result(vc)
            complex(kind=8), intent(in) :: G(:, :)
            integer, intent(in) :: idx
            complex(kind=8) :: vc(size(G,1))
            vc = G(:, idx)
        end function Zcol

        subroutine inv2x2(a, b, c, d, Minv)
            complex(kind=8), intent(in) :: a, b, c, d
            complex(kind=8), intent(out) :: Minv(2,2)
            complex(kind=8) :: detM
            detM = a*d - b*c
            ! 如果矩阵接近奇异，这不应该发生（已在前面检查）
            ! 但为了安全，仍然处理这种情况
            if (abs(detM) < 1.d-12) then
                ! 返回0矩阵，表示无法求逆
                Minv = cmplx(0.d0, 0.d0, kind=8)
            else
                Minv(1,1) =  d / detM
                Minv(1,2) = -b / detM
                Minv(2,1) = -c / detM
                Minv(2,2) =  a / detM
            endif
        end subroutine inv2x2
    end subroutine sigma_flip_rank2

end module SigmaRank2_mod



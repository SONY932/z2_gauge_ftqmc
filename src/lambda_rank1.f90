module LambdaRank1_mod
    implicit none
    
! 数值稳定性阈值：当 |r1| 小于此值时拒绝更新
    real(kind=8), parameter :: LAMBDA_RANK1_THRESHOLD = 1.d-12

contains
    subroutine lambda_flip_rank1(Gr, i, det_ratio, do_update)
! 末片 λ_i 翻转的 rank-1 快更新
!
! 理论背景：
!   当末片 λ_i 翻转时，格林函数需要进行 rank-1 更新。
!   设 r1 = 2*G_ii - 1，则：
!     - 单通道 det_ratio = |r1|^2
!     - G' = G + (2/r1) * G(:,i) * (e_i - G(i,:))
!
! 数值稳定性：
!   当 G_ii ≈ 0.5 时，r1 ≈ 0，更新公式不稳定。
!   此时 det_ratio ≈ 0，应该拒绝更新（概率为0）。
!
        complex(kind=8), intent(inout) :: Gr(:, :)
        integer, intent(in) :: i
        real(kind=8), intent(out) :: det_ratio
        logical, intent(in), optional :: do_update
        logical :: upd
        complex(kind=8) :: r1, factor
        complex(kind=8), dimension(size(Gr,1)) :: Gei
        complex(kind=8), dimension(size(Gr,2)) :: delta_row
        integer :: n, k
        real(kind=8) :: abs_r1

        upd = .true.
        if (present(do_update)) upd = do_update

! 计算 r1 = 2*G_ii - 1
        r1 = dcmplx(2.d0, 0.d0) * Gr(i, i) - dcmplx(1.d0, 0.d0)
        abs_r1 = abs(r1)
        
! det_ratio = |r1|^2（单通道）
        det_ratio = abs_r1 * abs_r1
        
! 数值稳定性检查：如果 r1 太小，拒绝更新
        if (abs_r1 < LAMBDA_RANK1_THRESHOLD) then
            det_ratio = 0.d0
            return
        endif
        
        if (.not. upd) return

! 执行 rank-1 更新：G' = G + factor * G(:,i) * (e_i - G(i,:))^T
! 其中 factor = 2/r1
        factor = dcmplx(2.d0, 0.d0) / r1
        
! 提取第 i 列
        Gei = Gr(:, i)
        
! 构造 (e_i - G(i,:))，即第 i 行取负号，第 i 个元素变成 1 - G_ii
        delta_row = -Gr(i, :)
        delta_row(i) = dcmplx(1.d0, 0.d0) - Gr(i, i)
        
! 执行 rank-1 更新：G' = G + factor * outer(Gei, delta_row)
! 为了数值稳定性，逐元素计算而不是使用 matmul
        do n = 1, size(Gr, 2)
            do k = 1, size(Gr, 1)
                Gr(k, n) = Gr(k, n) + factor * Gei(k) * delta_row(n)
            enddo
        enddo
        
        return
    end subroutine lambda_flip_rank1

end module LambdaRank1_mod



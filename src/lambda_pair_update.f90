module LambdaPair_mod
    use calc_basic
    use GaugeFields_mod
    use GaugeAction_mod
    use LambdaRank1_mod
    use MyLattice
    implicit none

contains
    logical function check_nan_in_row(Mat, row) result(has_nan)
! 检查矩阵的第row行和列是否包含NaN
        complex(kind=8), intent(in) :: Mat(:, :)
        integer, intent(in) :: row
        integer :: k
        real(kind=8) :: val
        has_nan = .false.
        do k = 1, size(Mat, 2)
            val = real(Mat(row, k)) + aimag(Mat(row, k))
            if (val /= val) then  ! NaN check
                has_nan = .true.
                return
            endif
            val = real(Mat(k, row)) + aimag(Mat(k, row))
            if (val /= val) then
                has_nan = .true.
                return
            endif
        enddo
        return
    end function check_nan_in_row

    subroutine lambda_pair_flip_v2(GrU, GrD, Gauge, i, j, Latt, iseed, acc)
! 【修复版】λ成对翻转，正确处理投影的撤回
! 关键点：撤回投影必须使用旧的λ值，而不是新的λ值
! 
! 设 G_proj = P_old G P_old（投影后的Green函数，输入）
! 如果翻转被接受：
!   1. 先用旧λ撤回投影：G = P_old G_proj P_old
!   2. 更新λ
!   3. 再用新λ施加投影：G_proj_new = P_new G P_new
! 这等价于：G_proj_new = P_new P_old G_proj P_old P_new = D G_proj D
! 其中 D = diag(1,...,-1,...,-1,...,1)，在i和j位置为-1
!
! 所以正确的更新是：对G_proj的第i和第j行列取负（除了对角元）
        complex(kind=8), intent(inout) :: GrU(:, :), GrD(:, :)
        class(GaugeConf), intent(inout) :: Gauge
        integer, intent(in) :: i, j
        class(SquareLattice), intent(in) :: Latt
        integer, intent(inout) :: iseed
        logical, intent(out) :: acc
        real(kind=8), external :: ranf
        real(kind=8) :: dS_Gauss, prob, thr
        real(kind=8) :: gam
        integer :: s_i_0, s_i_M, s_j_0, s_j_M
        integer :: k

        acc = .false.
        if (i == j) return
        
! Gauss边界项：计算翻转λ_i和λ_j的作用量变化
        gam = temporal_gamma()
        s_i_0 = star_product(Gauge, Latt, i, 1)
        s_i_M = star_product(Gauge, Latt, i, Ltrot)
        s_j_0 = star_product(Gauge, Latt, j, 1)
        s_j_M = star_product(Gauge, Latt, j, Ltrot)
        dS_Gauss = -2.d0 * gam * dble( Gauge%lambda(i) * s_i_0 * s_i_M + &
                                        Gauge%lambda(j) * s_j_0 * s_j_M )
        
! 数值保护：防止exp溢出
! dS_Gauss > 700 时 exp(-dS_Gauss) ≈ 0，直接拒绝
! dS_Gauss < -700 时 exp(-dS_Gauss) 溢出，直接接受
        if (dS_Gauss > 700.d0) then
            acc = .false.
            return
        endif
        if (dS_Gauss < -700.d0) then
            prob = 1.d0  ! 直接接受
        else
            prob = exp(-dS_Gauss)
        endif
        thr = ranf(iseed)
        
        if (min(1.d0, prob) > thr) then
! 接受：更新Green函数和λ
! 正确的更新是：G_proj_new = D G_proj_old D
! 其中 D(k,k) = -1 if k=i or k=j, else 1
! 这意味着：对第i和j行列取负（除了(i,i)和(j,j)以及(i,j)和(j,i)）

! 检查Green函数是否有效（防止传播NaN）
            if (check_nan_in_row(GrU, i) .or. check_nan_in_row(GrU, j) .or. &
                check_nan_in_row(GrD, i) .or. check_nan_in_row(GrD, j)) then
                acc = .false.
                return
            endif
            
! 更新 GrU
            do k = 1, Ndim
                if (k /= i .and. k /= j) then
                    GrU(i, k) = -GrU(i, k)
                    GrU(k, i) = -GrU(k, i)
                    GrU(j, k) = -GrU(j, k)
                    GrU(k, j) = -GrU(k, j)
                endif
            enddo
! (i,j) 和 (j,i) 位置：(-1)*(-1) = 1，不变
! (i,i) 和 (j,j) 位置：(-1)*(-1) = 1，不变
            
! 更新 GrD（同样的操作）
            do k = 1, Ndim
                if (k /= i .and. k /= j) then
                    GrD(i, k) = -GrD(i, k)
                    GrD(k, i) = -GrD(k, i)
                    GrD(j, k) = -GrD(j, k)
                    GrD(k, j) = -GrD(k, j)
                endif
            enddo
            
! 翻转λ
            Gauge%lambda(i) = -Gauge%lambda(i)
            Gauge%lambda(j) = -Gauge%lambda(j)
            acc = .true.
        endif
        return
    end subroutine lambda_pair_flip_v2

    subroutine lambda_pair_flip(GrU, GrD, Gauge, i, j, Latt, iseed, acc)
! 随机成对翻转 λ_i 与 λ_j，保持 ∏λ=Q
! 接受率包含：费米子因子 Rf + Gauss边界项 exp(-ΔS_Gauss)
        complex(kind=8), intent(inout) :: GrU(:, :), GrD(:, :)
        class(GaugeConf), intent(inout) :: Gauge
        integer, intent(in) :: i, j
        class(SquareLattice), intent(in) :: Latt
        integer, intent(inout) :: iseed
        logical, intent(out) :: acc
        real(kind=8), external :: ranf
        real(kind=8) :: rU_i, rD_i, rU_j, rD_j
        real(kind=8) :: Rf, dS_Gauss, prob, thr
        real(kind=8) :: gam
        integer :: s_i_0, s_i_M, s_j_0, s_j_M
        complex(kind=8), allocatable :: GUtmp(:, :), GDtmp(:, :)

        acc = .false.
! 若 i==j 则跳过
        if (i == j) return
        
! Gauss边界项
        gam = temporal_gamma()
        s_i_0 = star_product(Gauge, Latt, i, 1)
        s_i_M = star_product(Gauge, Latt, i, Ltrot)
        s_j_0 = star_product(Gauge, Latt, j, 1)
        s_j_M = star_product(Gauge, Latt, j, Ltrot)
        dS_Gauss = -2.d0 * gam * dble( Gauge%lambda(i) * s_i_0 * s_i_M + &
                                        Gauge%lambda(j) * s_j_0 * s_j_M )
        
! 计算两点 rank-1 的乘积（先不更新，确保回滚容易）
        allocate(GUtmp(Ndim, Ndim)); allocate(GDtmp(Ndim, Ndim))
        GUtmp = GrU; GDtmp = GrD
        call lambda_flip_rank1(GUtmp, i, rU_i, .false.)
        call lambda_flip_rank1(GDtmp, i, rD_i, .false.)
        call lambda_flip_rank1(GUtmp, j, rU_j, .false.)
        call lambda_flip_rank1(GDtmp, j, rD_j, .false.)
        rU_i = max(rU_i, 0.d0)
        rD_i = max(rD_i, 0.d0)
        rU_j = max(rU_j, 0.d0)
        rD_j = max(rD_j, 0.d0)
        Rf = max(rU_i, 1.d-300) * max(rD_i, 1.d-300) * max(rU_j, 1.d-300) * max(rD_j, 1.d-300)
        
! 总接受概率：p = min{1, exp(-ΔS_Gauss) * Rf}
        prob = exp(-dS_Gauss) * Rf
        thr = ranf(iseed)
        if (min(1.d0, prob) > thr) then
! 接受：就地更新 G 并翻转 λ
            call lambda_flip_rank1(GrU, i, rU_i, .true.)
            call lambda_flip_rank1(GrD, i, rD_i, .true.)
            call lambda_flip_rank1(GrU, j, rU_j, .true.)
            call lambda_flip_rank1(GrD, j, rD_j, .true.)
            Gauge%lambda(i) = - Gauge%lambda(i)
            Gauge%lambda(j) = - Gauge%lambda(j)
            acc = .true.
        endif
        deallocate(GUtmp); deallocate(GDtmp)
        return
    end subroutine lambda_pair_flip

end module LambdaPair_mod



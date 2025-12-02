module LambdaPair_mod
    use calc_basic
    use GaugeFields_mod
    use GaugeAction_mod
    use LambdaRank1_mod
    use MyLattice
    implicit none

contains
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



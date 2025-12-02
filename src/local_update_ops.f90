module LocalUpdateOps_mod
    use calc_basic
    use ProcessMatrix
    use MyLattice
    use BondList_mod
    use GaugeFields_mod
    use LocalSigma_mod
    use OpMu_mod
    implicit none
contains

    subroutine prop_pre_step(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt)
        call opMu_mmult_R(PropU%UUR, -1)
        call opMu_mmult_R(PropD%UUR, -1)
    end subroutine prop_pre_step

    subroutine left_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
! 左扫时的预处理（空操作，因为现在用完整层传播）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        ! 什么都不做，等待 post_mu 做完整传播
    end subroutine left_step_prefix

    subroutine left_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 左扫σ更新后（空操作）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        ! 什么都不做
    end subroutine left_step_after_sigma_pre_mu

    subroutine left_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 累积 UUL: UUL <- UUL * B(nt)（右乘正向）
! 这样 UUL(nt) = B(Ltrot) * ... * B(nt+1)，与 stab_green 一致
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 累积 UUL: UUL <- UUL * B = UUL * B_gauge * e^{-μ}
        call apply_trotter_layer_R(PropU%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call opMu_mmult_R(PropU%UUL, -1)
        call opMu_mmult_R(PropD%UUL, -1)
    end subroutine left_step_after_sigma_post_mu

    subroutine right_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫预处理（空操作）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
    end subroutine right_step_prefix

    subroutine right_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫σ更新后（空操作）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
    end subroutine right_step_after_sigma_pre_mu

    subroutine right_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫（空操作，传播在 propagate_right_step 中完成）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
    end subroutine right_step_after_sigma_post_mu

    subroutine right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
! 累积 UUR: UUR <- UUR * B
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
    end subroutine right_step_finalize

end module LocalUpdateOps_mod

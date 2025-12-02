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
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_L(PropU%Gr, Latt, Bonds, Gauge, nt)
        call apply_half_trotter_L(PropD%Gr, Latt, Bonds, Gauge, nt)
    end subroutine left_step_prefix

    subroutine left_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_L(PropU%Gr, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_L(PropD%Gr, Latt, Bonds, Gauge, nt, reverse=.true.)
    end subroutine left_step_after_sigma_pre_mu

    subroutine left_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call opMu_mmult_L(PropU%UUL, +1)
        call opMu_mmult_L(PropD%UUL, +1)
        call apply_half_trotter_L(PropU%UUL, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_L(PropD%UUL, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_L(PropU%UUL, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_L(PropD%UUL, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
    end subroutine left_step_after_sigma_post_mu

    subroutine right_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
    end subroutine right_step_prefix

    subroutine right_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
    end subroutine right_step_after_sigma_pre_mu

    subroutine right_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_L(PropU%Gr, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_L(PropD%Gr, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_L(PropU%Gr, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_L(PropD%Gr, Latt, Bonds, Gauge, nt, reverse=.true.)
    end subroutine right_step_after_sigma_post_mu

    subroutine right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_R(PropU%UUR, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_R(PropD%UUR, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_R(PropU%UUR, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_R(PropD%UUR, Latt, Bonds, Gauge, nt, reverse=.true.)
    end subroutine right_step_finalize

end module LocalUpdateOps_mod

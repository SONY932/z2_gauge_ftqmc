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
! 预热阶段累积 UUR: UUR <- B(nt) * UUR
! B(nt) = exp(-μ) * B_gauge
! 使用左乘累积：UUR <- exp(-μ) * (B_gauge * UUR)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! UUR <- B_gauge * UUR (左乘)
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
! UUR <- exp(-μ) * UUR (左乘)
        call opMu_mmult_R(PropU%UUR, -1)
        call opMu_mmult_R(PropD%UUR, -1)
    end subroutine prop_pre_step

! ===== 左扫步骤 =====
    subroutine left_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        ! 空操作
    end subroutine left_step_prefix

    subroutine left_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        ! 空操作
    end subroutine left_step_after_sigma_pre_mu

    subroutine left_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 左扫：G(nt-1) = B(nt)^{-1} * G(nt) * B(nt)
! 同时累积 UUL: UUL <- UUL * B(nt)
! B(nt) = exp(-μ) * B_gauge
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt

! G 传播：G <- B(nt)^{-1} * G * B(nt)
! B^{-1} = B_gauge^{-1} * exp(μ)
! 步骤 1: G <- B_gauge^{-1} * G (左乘逆)
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤 2: G <- exp(μ) * G (左乘)
        call opMu_mmult_R(PropU%Gr, +1)
        call opMu_mmult_R(PropD%Gr, +1)
! 步骤 3: G <- G * B_gauge (右乘)
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤 4: G <- G * exp(-μ) (右乘)
        call opMu_mmult_L(PropU%Gr, -1)
        call opMu_mmult_L(PropD%Gr, -1)

! UUL 累积：UUL <- UUL * B(nt) = UUL * exp(-μ) * B_gauge
! 使用右乘累积
! 步骤 1: UUL <- UUL * exp(-μ) (右乘)
        call opMu_mmult_L(PropU%UUL, -1)
        call opMu_mmult_L(PropD%UUL, -1)
! 步骤 2: UUL <- UUL * B_gauge (右乘)
        call apply_trotter_layer_L(PropU%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_L(PropD%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
    end subroutine left_step_after_sigma_post_mu

! ===== 右扫步骤 =====
    subroutine right_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        ! 空操作
    end subroutine right_step_prefix

    subroutine right_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        ! 空操作
    end subroutine right_step_after_sigma_pre_mu

    subroutine right_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫：G(nt) = B(nt)^{-1} * G(nt-1) * B(nt)
! B(nt) = exp(-μ) * B_gauge
! B^{-1} = B_gauge^{-1} * exp(μ)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt

! G 传播：G <- B(nt)^{-1} * G * B(nt)
! 步骤 1: G <- B_gauge^{-1} * G (左乘逆)
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤 2: G <- exp(μ) * G (左乘)
        call opMu_mmult_R(PropU%Gr, +1)
        call opMu_mmult_R(PropD%Gr, +1)
! 步骤 3: G <- G * B_gauge (右乘)
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤 4: G <- G * exp(-μ) (右乘)
        call opMu_mmult_L(PropU%Gr, -1)
        call opMu_mmult_L(PropD%Gr, -1)
    end subroutine right_step_after_sigma_post_mu

    subroutine right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
! UUR 累积：UUR <- B(nt) * UUR
! 使用左乘累积
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 步骤 1: UUR <- B_gauge * UUR (左乘)
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤 2: UUR <- exp(-μ) * UUR (左乘)
        call opMu_mmult_R(PropU%UUR, -1)
        call opMu_mmult_R(PropD%UUR, -1)
    end subroutine right_step_finalize

end module LocalUpdateOps_mod

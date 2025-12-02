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
! 预热累积：与原始代码保持一致
! 使用 apply_trotter_layer_R（实际是右乘）和 opMu_mmult_L（实际是右乘）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt)
        call opMu_mmult_L(PropU%UUR, -1)
        call opMu_mmult_L(PropD%UUR, -1)
    end subroutine prop_pre_step

    subroutine left_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，不再需要
    end subroutine left_step_prefix

    subroutine left_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，不再需要
    end subroutine left_step_after_sigma_pre_mu

    subroutine left_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 左扫传播：G(nt-1) = B(nt)^{-1} * G(nt) * B(nt)
! 其中 B = B_gauge * exp(-μ)，B^{-1} = exp(μ) * B_gauge^{-1}
! 完整公式：G' = exp(μ) * B_gauge^{-1} * G * B_gauge * exp(-μ)
! 同时累积 UUL: UUL <- UUL * B(nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 传播 G：按从右到左的顺序应用
! 步骤1：G <- G * B_gauge（右乘）
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤2：G <- G * exp(-μ)（右乘）
        call opMu_mmult_L(PropU%Gr, -1)
        call opMu_mmult_L(PropD%Gr, -1)
! 步骤3：G <- B_gauge^{-1} * G（左乘）
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤4：G <- exp(μ) * G（左乘）
        call opMu_mmult_R(PropU%Gr, +1)
        call opMu_mmult_R(PropD%Gr, +1)
! 累积 UUL: UUL <- UUL * B(nt)
! 步骤1：UUL <- UUL * B_gauge（右乘）
        call apply_trotter_layer_R(PropU%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤2：UUL <- UUL * exp(-μ)（右乘）
        call opMu_mmult_L(PropU%UUL, -1)
        call opMu_mmult_L(PropD%UUL, -1)
    end subroutine left_step_after_sigma_post_mu

    subroutine right_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，不再需要
    end subroutine right_step_prefix

    subroutine right_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，不再需要
    end subroutine right_step_after_sigma_pre_mu

    subroutine right_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫传播：G(nt+1) = B(nt) * G(nt) * B(nt)^{-1}
! 其中 B = B_gauge * exp(-μ)，B^{-1} = exp(μ) * B_gauge^{-1}
! 完整公式：G' = B_gauge * exp(-μ) * G * exp(μ) * B_gauge^{-1}
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 传播 G：按从右到左的顺序应用
! 步骤1：G <- G * B_gauge^{-1}（右乘）
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤2：G <- G * exp(μ)（右乘）
        call opMu_mmult_L(PropU%Gr, +1)
        call opMu_mmult_L(PropD%Gr, +1)
! 步骤3：G <- exp(-μ) * G（左乘）
        call opMu_mmult_R(PropU%Gr, -1)
        call opMu_mmult_R(PropD%Gr, -1)
! 步骤4：G <- B_gauge * G（左乘）
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
    end subroutine right_step_after_sigma_post_mu

    subroutine right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫累积：与原始代码保持一致
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call opMu_mmult_L(PropU%UUR, -1)
        call opMu_mmult_L(PropD%UUR, -1)
    end subroutine right_step_finalize

end module LocalUpdateOps_mod

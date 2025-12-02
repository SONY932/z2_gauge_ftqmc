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
! 预热累积 UUR: UUR <- B(nt) * UUR（左乘）
! 注意：由于名字混乱，apply_trotter_layer_R 实际执行左乘！
! 从 nt=1 到 Ltrot 累积后得到 UUR = B(Ltrot) * ... * B(1)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! B = B_gauge * exp(-μ)
! 左乘 B：UUR <- B * UUR = B_gauge * exp(-μ) * UUR
! 由于 exp(-μ) 是标量，可以按任意顺序应用
! apply_trotter_layer_R 执行左乘 B_gauge
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt)
! opMu_mmult_R 执行左乘 exp(-μ)
        call opMu_mmult_R(PropU%UUR, -1)
        call opMu_mmult_R(PropD%UUR, -1)
    end subroutine prop_pre_step

    subroutine left_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，暂不使用
    end subroutine left_step_prefix

    subroutine left_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，暂不使用
    end subroutine left_step_after_sigma_pre_mu

    subroutine left_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 左扫传播（关键修正！）：G(nt-1) = B(nt) * G(nt) * B(nt)^{-1}
! 其中 B = B_gauge * exp(-μ)
!
! 数学推导：
! G(τ) = (1 + B_R(τ) * B_L(τ))^{-1}
! B_R(τ-1) = B(τ)^{-1} * B_R(τ)
! B_L(τ-1) = B_L(τ) * B(τ)
! => G(τ-1) = (1 + B(τ)^{-1} * B_R * B_L * B(τ))^{-1}
!           = B(τ) * (1 + B_R * B_L)^{-1} * B(τ)^{-1}
!           = B(τ) * G(τ) * B(τ)^{-1}
!
! 步骤：
! 1. G <- B * G = B_gauge * exp(-μ) * G   （左乘 B）
! 2. G <- G * B^{-1} = G * exp(μ) * B_gauge^{-1}   （右乘 B^{-1}）
!
! 累积 UUL: UUL <- B(nt) * UUL （左乘）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt

! 传播 G：
! 步骤1：G <- B_gauge * G（左乘）
! 注意：apply_trotter_layer_R 实际执行左乘！
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤2：G <- exp(-μ) * G（左乘）
        call opMu_mmult_R(PropU%Gr, -1)
        call opMu_mmult_R(PropD%Gr, -1)
! 步骤3：G <- G * B_gauge^{-1}（右乘）
! 注意：apply_trotter_layer_L 实际执行右乘！
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤4：G <- G * exp(μ)（右乘）
        call opMu_mmult_L(PropU%Gr, +1)
        call opMu_mmult_L(PropD%Gr, +1)

! 累积 UUL: UUL <- B(nt) * UUL（左乘）
! 注意：虽然理论上 B_L(nt) = B(Ltrot) * ... * B(nt+1) 需要右乘累积，
! 但 stab_UL 内部处理了这种差异（通过 UUL^H）。
! 保持与代码其他部分的一致性，使用左乘累积。
! 使用 apply_trotter_layer_R（实际执行左乘）
! 步骤1：UUL <- B_gauge * UUL（左乘）
        call apply_trotter_layer_R(PropU%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUL, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤2：UUL <- exp(-μ) * UUL（左乘）
        call opMu_mmult_R(PropU%UUL, -1)
        call opMu_mmult_R(PropD%UUL, -1)
    end subroutine left_step_after_sigma_post_mu

    subroutine right_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，暂不使用
    end subroutine right_step_prefix

    subroutine right_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 为 σ 更新做准备，暂不使用
    end subroutine right_step_after_sigma_pre_mu

    subroutine right_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫传播（关键修正！）：G(nt) = B(nt)^{-1} * G(nt-1) * B(nt)
! 其中 B = B_gauge * exp(-μ)
!
! 数学推导：
! G(τ) = (1 + B_R(τ) * B_L(τ))^{-1}
! B_R(τ+1) = B(τ+1) * B_R(τ)
! B_L(τ+1) = B_L(τ) * B(τ+1)^{-1}
! => G(τ+1) = (1 + B(τ+1) * B_R * B_L * B(τ+1)^{-1})^{-1}
!           = B(τ+1)^{-1} * (1 + B_R * B_L)^{-1} * B(τ+1)
!           = B(τ+1)^{-1} * G(τ) * B(τ+1)
!
! 步骤：
! 1. G <- B^{-1} * G = exp(μ) * B_gauge^{-1} * G   （左乘 B^{-1}）
! 2. G <- G * B = G * B_gauge * exp(-μ)   （右乘 B）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt

! 传播 G：
! 步骤1：G <- B_gauge^{-1} * G（左乘）
! 注意：apply_trotter_layer_R 实际执行左乘！
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤2：G <- exp(μ) * G（左乘）
        call opMu_mmult_R(PropU%Gr, +1)
        call opMu_mmult_R(PropD%Gr, +1)
! 步骤3：G <- G * B_gauge（右乘）
! 注意：apply_trotter_layer_L 实际执行右乘！
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤4：G <- G * exp(-μ)（右乘）
        call opMu_mmult_L(PropU%Gr, -1)
        call opMu_mmult_L(PropD%Gr, -1)
    end subroutine right_step_after_sigma_post_mu

    subroutine right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫累积 UUR: UUR <- B(nt) * UUR（左乘）
! 与预热相同的逻辑
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call opMu_mmult_R(PropU%UUR, -1)
        call opMu_mmult_R(PropD%UUR, -1)
    end subroutine right_step_finalize

end module LocalUpdateOps_mod

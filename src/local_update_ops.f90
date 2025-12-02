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
! UUR 累积：UUR <- B(nt) * UUR = exp(-μ) * B_gauge * UUR
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

! ============ 左扫传播函数（简化版，不支持 σ 更新）============
! 左扫公式：G(τ-1) = B(τ)^{-1} * G(τ) * B(τ)
! 其中 B = exp(-μ) * B_gauge，B^{-1} = B_gauge^{-1} * exp(μ)
! UUL 累积：UUL <- UUL * B = UUL * exp(-μ) * B_gauge

    subroutine left_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
! 占位，保持接口兼容
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 不做任何操作，所有传播逻辑在 left_step_after_sigma_post_mu 中
    end subroutine left_step_prefix

    subroutine left_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 占位，保持接口兼容
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 不做任何操作
    end subroutine left_step_after_sigma_pre_mu

    subroutine left_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 【关键修复】使用完整 Trotter 层，确保与 UUR 累积一致
! G 传播：G(τ-1) = B_gauge^{-1} * G(τ) * B_gauge
! 1. G <- B_gauge^{-1} * G（左乘逆）
! 2. G <- G * B_gauge（右乘）
! UUL 累积：UUL <- UUL * B_gauge（右乘，μ 在 local_sweep 中处理）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt

! UUL 累积（右乘 B_gauge）
        call apply_trotter_layer_L(PropU%UUL, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_L(PropD%UUL, Latt, Bonds, Gauge, nt)

! G 传播：G <- B_gauge^{-1} * G * B_gauge
! 步骤1：G <- B_gauge^{-1} * G（左乘逆）
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤2：G <- G * B_gauge（右乘）
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt)
    end subroutine left_step_after_sigma_post_mu

! ============ 右扫传播函数（简化版，不支持 σ 更新）============
! 右扫公式：G(τ) = B(τ) * G(τ-1) * B(τ)^{-1}
! UUR 累积：UUR <- B(τ) * UUR

    subroutine right_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
! 占位，保持接口兼容
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 不做任何操作
    end subroutine right_step_prefix

    subroutine right_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 占位，保持接口兼容
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! 不做任何操作
    end subroutine right_step_after_sigma_pre_mu

    subroutine right_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 【关键修复】使用完整 Trotter 层
! G 传播：G(τ) = B_gauge * G(τ-1) * B_gauge^{-1}
! 1. G <- B_gauge * G（左乘）
! 2. G <- G * B_gauge^{-1}（右乘逆）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt

! G 传播：G <- B_gauge * G * B_gauge^{-1}
! 步骤1：G <- B_gauge * G（左乘）
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt)
! 步骤2：G <- G * B_gauge^{-1}（右乘逆）
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
    end subroutine right_step_after_sigma_post_mu

    subroutine right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
! UUR 累积：UUR <- B_gauge * UUR（左乘，μ 在 local_sweep 中处理）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt)
    end subroutine right_step_finalize

end module LocalUpdateOps_mod

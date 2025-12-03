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
! 预传播：UUR <- B * UUR，其中 B = exp(+μ) * B_gauge
! 【重要】μ 符号为 +1，与 CodeXun 的 propT_pre 一致
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_trotter_layer_R(PropU%UUR, Latt, Bonds, Gauge, nt)
        call apply_trotter_layer_R(PropD%UUR, Latt, Bonds, Gauge, nt)
        call opMu_mmult_R(PropU%UUR, +1)
        call opMu_mmult_R(PropD%UUR, +1)
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
! 左扫：G <- B^{-1} * G * B，UUL <- UUL * B
! B = half_forward * half_reverse（与右乘一致）
! B^{-1} = half_reverse^{-1} * half_forward^{-1}
! 对于左乘（先调用的在左边）：先 false 再 true -> half_reverse^{-1} * half_forward^{-1} * G = B^{-1} * G
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
! UUL 累积
        call opMu_mmult_L(PropU%UUL, +1)
        call opMu_mmult_L(PropD%UUL, +1)
        call apply_half_trotter_L(PropU%UUL, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_L(PropD%UUL, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_L(PropU%UUL, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_L(PropD%UUL, Latt, Bonds, Gauge, nt, reverse=.true.)
! G 逆操作：G <- B^{-1} * G
! 【关键】顺序：先 false 再 true -> half_reverse^{-1} * half_forward^{-1} * G = B^{-1} * G
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
    end subroutine left_step_after_sigma_post_mu

    subroutine left_step_after_sigma_post_mu_G_only(PropU, PropD, Latt, Bonds, Gauge, nt)
! 只处理 G：G <- B^{-1} * G
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
    end subroutine left_step_after_sigma_post_mu_G_only

    subroutine left_step_after_sigma_post_mu_UUL_only(PropU, PropD, Latt, Bonds, Gauge, nt)
! 只处理 UUL：UUL <- UUL * B（使用更新后的σ）
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
    end subroutine left_step_after_sigma_post_mu_UUL_only

    subroutine right_step_prefix(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫：G <- B * G * B^{-1}，其中 B = half_forward * half_reverse（与左扫一致）
! 步骤 1：G <- half_reverse * G（左乘正向，nflag=+1）
! 注意：apply_half_trotter_R 是左乘，所以先调用的在左边
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1, reverse=.true.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1, reverse=.true.)
    end subroutine right_step_prefix

    subroutine right_step_after_sigma_pre_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫：G <- B * G * B^{-1}
! 步骤 2：G <- half_forward * half_reverse * G = B * G（左乘正向）
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1, reverse=.false.)
        call apply_half_trotter_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1, reverse=.false.)
    end subroutine right_step_after_sigma_pre_mu

    subroutine right_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫：G <- B * G * B^{-1}
! 步骤 3：G <- G * B^{-1} = G * half_reverse^{-1} * half_forward^{-1}（右乘逆）
! B^{-1} = (half_forward * half_reverse)^{-1} = half_reverse^{-1} * half_forward^{-1}
! 对于右乘，先调用的在右边，所以先 true 再 false
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
        call apply_half_trotter_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.true.)
        call apply_half_trotter_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
        call apply_half_trotter_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1, reverse=.false.)
    end subroutine right_step_after_sigma_post_mu

    subroutine right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
! 右扫 UUR 累积：UUR <- B * UUR（左乘）
! B = half_forward * half_reverse（与 G 传播一致）
! 先调用的在左边：先 true 再 false -> half_forward * half_reverse * UUR
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        call apply_half_trotter_R(PropU%UUR, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_R(PropD%UUR, Latt, Bonds, Gauge, nt, reverse=.true.)
        call apply_half_trotter_R(PropU%UUR, Latt, Bonds, Gauge, nt, reverse=.false.)
        call apply_half_trotter_R(PropD%UUR, Latt, Bonds, Gauge, nt, reverse=.false.)
    end subroutine right_step_finalize

end module LocalUpdateOps_mod

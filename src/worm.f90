module Worm_mod
    use calc_basic
    use MyLattice
    use GaugeFields_mod
    use BondList_mod
    use GaugeAction_mod
    use SigmaRank2_mod
    use OpMu_mod
    use MultiplyChecker_mod
    implicit none

    type, public :: WormStats
        real(kind=8) :: accepted
        real(kind=8) :: proposed
    contains
        procedure :: reset => WormStats_reset
        procedure :: ratio => WormStats_ratio
    end type WormStats

contains
    subroutine WormStats_reset(this)
        class(WormStats), intent(inout) :: this
        this%accepted = 0.d0
        this%proposed = 0.d0
        return
    end subroutine WormStats_reset

    subroutine WormStats_ratio(this)
        class(WormStats), intent(inout) :: this
        if (this%proposed > 0.d0) this%accepted = this%accepted / this%proposed
        return
    end subroutine WormStats_ratio

    subroutine worm_time_string(GrU, GrD, Gauge, Latt, Bonds, nt0, bdir, iseed, acc)
! 翻转随机选择的一条键 bdir ∈ {'x','y'} 在随机区间 [τ1,τ2] 的全部 σ
        complex(kind=8), intent(inout) :: GrU(:, :), GrD(:, :)
        class(GaugeConf), intent(inout) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        integer, intent(in) :: nt0
        character(len=*), intent(in) :: bdir
        integer, intent(inout) :: iseed
        logical, intent(out) :: acc
        real(kind=8), external :: ranf
        integer :: tau1, tau2, t, tt, step
        integer :: bidx, i_src, j_dst
        integer :: x, y, xp1, xm1, yp1, ym1, ntp1, ntm1
        real(kind=8) :: dS_tot, dS_spa, dS_time, dS_bdry, log_Rf_tot
        real(kind=8) :: r_up, r_dn, thr, val, logA, gam
        integer :: boundary_sites(2)
        include 'mpif.h'
        complex(kind=8), allocatable :: GUtmp(:, :), GDtmp(:, :)

        acc = .false.
! 从 wrap 边界 nt0 出发：区间限制在该段内，避免跨段重建
        tau1 = nt0 + (nranf(iseed, max(1, Nwrap)) - 1)
        if (tau1 > nt0 + Nwrap - 1) tau1 = nt0 + Nwrap - 1
        if (tau1 > Ltrot) tau1 = Ltrot
! 以小跨度优先：长度在 [1, min(worm_max_span, 段剩余长度)] 内均匀采样
        step = min(worm_max_span, nt0 + Nwrap - 1 - tau1)
        if (step < 1) step = 1
        tau2 = min(nt0 + Nwrap - 1, tau1 + (nranf(iseed, step) - 1))
        if (tau2 > Ltrot) tau2 = Ltrot
        
! 特殊处理：如果nt0是最后一个Wrap段(>=Ltrot)，避免在边界处更新
        if (nt0 >= Ltrot) then
            acc = .false.
            return
        endif

! 随机挑一条对应方向的键
        if (bdir == 'x') then
            bidx = nranf(iseed, size(Bonds%ex_src))
            i_src = Bonds%ex_src(bidx)
            j_dst = Bonds%ex_dst(bidx)
        elseif (bdir == 'y') then
            bidx = nranf(iseed, size(Bonds%ey_src))
            i_src = Bonds%ey_src(bidx)
            j_dst = Bonds%ey_dst(bidx)
        else
            write(6,*) 'worm_time_string: illegal bdir=', bdir; return
        endif
        if (bdir == 'x') then
            boundary_sites = (/ i_src, Latt%L_bonds(i_src, 1) /)
        else
            boundary_sites = (/ i_src, Latt%L_bonds(i_src, 2) /)
        endif

! 临时 Green 以顺次低秩更新计算 det-ratio；拒绝时丢弃
        allocate(GUtmp(Ndim, Ndim)); allocate(GDtmp(Ndim, Ndim))
        GUtmp = GrU; GDtmp = GrD
        dS_tot = 0.d0; dS_spa = 0.d0; dS_time = 0.d0; dS_bdry = 0.d0
        log_Rf_tot = 0.d0
        gam = temporal_gamma()

        x = Latt%n_list(i_src, 1); y = Latt%n_list(i_src, 2)
        xp1 = npbc(x+1, Nlx); xm1 = npbc(x-1, Nlx)
        yp1 = npbc(y+1, Nly); ym1 = npbc(y-1, Nly)
        ntm1 = npbc(tau1-1, Ltrot); ntp1 = npbc(tau2+1, Ltrot)

! 计算空间作用量变化 ΔS = S_new - S_old
! 当一条边翻转时，包含它的每个 plaquette 都会改变符号
! 对每个时间片，有两个 plaquette 受影响
        do t = tau1, tau2
            if (bdir == 'x') then
                ! 翻转 σx(x,y,t) 影响两个 plaquette：以 (x,y) 和 (x,y-1) 为左下角
                ! ΔS = -J * (p_new - p_old) = -J * (-p_old - p_old) = 2J * p_old
                dS_spa = dS_spa + 2.d0 * J * dble( Gauge%sigma_x(Latt%inv_n_list(x,   y),   t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(xp1, y),   t) * &
                                                 Gauge%sigma_x(Latt%inv_n_list(x,   yp1), t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(x,   y),   t) )
                dS_spa = dS_spa + 2.d0 * J * dble( Gauge%sigma_x(Latt%inv_n_list(x,   ym1), t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(xp1, ym1), t) * &
                                                 Gauge%sigma_x(Latt%inv_n_list(x,   y),   t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(x,   ym1), t) )
            else
                ! 翻转 σy(x,y,t) 影响两个 plaquette：以 (x,y) 和 (x-1,y) 为左下角
                dS_spa = dS_spa + 2.d0 * J * dble( Gauge%sigma_x(Latt%inv_n_list(x,   y),   t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(x,   y),   t) * &
                                                 Gauge%sigma_x(Latt%inv_n_list(x,   yp1), t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(xp1, y),   t) )
                dS_spa = dS_spa + 2.d0 * J * dble( Gauge%sigma_x(Latt%inv_n_list(xm1, y),   t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(xm1, y),   t) * &
                                                 Gauge%sigma_x(Latt%inv_n_list(xm1, yp1), t) * &
                                                 Gauge%sigma_y(Latt%inv_n_list(x,   y),   t) )
            endif
        enddo

        if (tau1 <= 1 .and. 1 <= tau2) then
            dS_bdry = dS_bdry + delta_S_lambda_boundary(Gauge, Latt, boundary_sites)
        endif
        if (tau1 <= Ltrot .and. Ltrot <= tau2) then
            dS_bdry = dS_bdry + delta_S_lambda_boundary(Gauge, Latt, boundary_sites)
        endif

! 时间耦合变化：ΔS_time = -2γ * [σ(τ1)*σ(τ1-1) + σ(τ2)*σ(τ2+1)]
! 标准周期边界条件，不涉及λ（λ只影响费米子部分）
        if (bdir == 'x') then
            dS_time = dS_time - 2.d0 * gam * dble( Gauge%sigma_x(i_src, tau1) * Gauge%sigma_x(i_src, ntm1) )
            dS_time = dS_time - 2.d0 * gam * dble( Gauge%sigma_x(i_src, tau2) * Gauge%sigma_x(i_src, ntp1) )
        else
            dS_time = dS_time - 2.d0 * gam * dble( Gauge%sigma_y(i_src, tau1) * Gauge%sigma_y(i_src, ntm1) )
            dS_time = dS_time - 2.d0 * gam * dble( Gauge%sigma_y(i_src, tau2) * Gauge%sigma_y(i_src, ntp1) )
        endif

! 然后再顺次进行低秩费米子更新，计算 log_Rf_tot，并就地翻转 σ
! 先把 G 从 nt0 传播到 tau1
! 重要修复：当 nt0 = tau1 时，循环不执行，格林函数已经在正确位置
! 右扫描传播公式：G(t+1) = B(t) * G(t) * B(t)^{-1}
! 其中 B = exp(+μ) * B_gauge
        if (nt0 < tau1) then
            do t = nt0, tau1 - 1
! 步骤1：G <- B_gauge * G（左乘 B_gauge）
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
! 步骤2：G <- exp(+μ) * G
                call opMu_mmult_R(GUtmp, +1); call opMu_mmult_R(GDtmp, +1)
! 步骤3：G <- G * B_gauge^{-1}（右乘 B_gauge 的逆）
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
! 步骤4：G <- G * exp(-μ)
                call opMu_mmult_L(GUtmp, -1); call opMu_mmult_L(GDtmp, -1)
            enddo
        endif

        do t = tau1, tau2
            if (bdir == 'x') then
                ! 先测试det_ratio
                call sigma_flip_rank2(GUtmp, i_src, j_dst, Gauge%sigma_x(i_src, t), r_up, .false.)
                call sigma_flip_rank2(GDtmp, i_src, j_dst, Gauge%sigma_x(i_src, t), r_dn, .false.)
                ! 如果 det_ratio 为0或负或极小（物理错误或数值不稳定），直接拒绝更新
                ! 使用更安全的阈值：避免log(val)太小导致数值溢出
                if (r_up <= 1.d-100 .or. r_dn <= 1.d-100) then
                    acc = .false.
                    ! 回滚已经翻转的 sigma
                    do tt = tau1, t-1
                        Gauge%sigma_x(i_src, tt) = - Gauge%sigma_x(i_src, tt)
                    enddo
                    deallocate(GUtmp); deallocate(GDtmp)
                    return
                endif
                ! det_ratio可接受，执行更新
                call sigma_flip_rank2(GUtmp, i_src, j_dst, Gauge%sigma_x(i_src, t), r_up, .true.)
                call sigma_flip_rank2(GDtmp, i_src, j_dst, Gauge%sigma_x(i_src, t), r_dn, .true.)
                val = r_up * r_dn
                ! 使用更安全的下限：log(1e-100) ≈ -230，即使累积10次也不会溢出
                log_Rf_tot = log_Rf_tot + log(max(val, 1.d-100))
                Gauge%sigma_x(i_src, t) = - Gauge%sigma_x(i_src, t)
            else
                ! 先测试det_ratio
                call sigma_flip_rank2(GUtmp, i_src, j_dst, Gauge%sigma_y(i_src, t), r_up, .false.)
                call sigma_flip_rank2(GDtmp, i_src, j_dst, Gauge%sigma_y(i_src, t), r_dn, .false.)
                ! 如果 det_ratio 为0或负或极小（物理错误或数值不稳定），直接拒绝更新
                ! 使用更安全的阈值：避免log(val)太小导致数值溢出
                if (r_up <= 1.d-100 .or. r_dn <= 1.d-100) then
                    acc = .false.
                    ! 回滚已经翻转的 sigma
                    do tt = tau1, t-1
                        Gauge%sigma_y(i_src, tt) = - Gauge%sigma_y(i_src, tt)
                    enddo
                    deallocate(GUtmp); deallocate(GDtmp)
                    return
                endif
                ! det_ratio可接受，执行更新
                call sigma_flip_rank2(GUtmp, i_src, j_dst, Gauge%sigma_y(i_src, t), r_up, .true.)
                call sigma_flip_rank2(GDtmp, i_src, j_dst, Gauge%sigma_y(i_src, t), r_dn, .true.)
                val = r_up * r_dn
                ! 使用更安全的下限：log(1e-100) ≈ -230，即使累积10次也不会溢出
                log_Rf_tot = log_Rf_tot + log(max(val, 1.d-100))
                Gauge%sigma_y(i_src, t) = - Gauge%sigma_y(i_src, t)
            endif
! 若不是最后一个时片，推进到下一时片（使用更新后的 σ）
! 右扫描传播公式：G(t+1) = B(t) * G(t) * B(t)^{-1}
            if (t < tau2) then
! 步骤1：G <- B_gauge * G（左乘 B_gauge）
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', +1, 0.5d0)
                call apply_group_R(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
                call apply_group_R(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', +1, 0.5d0)
! 步骤2：G <- exp(+μ) * G
                call opMu_mmult_R(GUtmp, +1); call opMu_mmult_R(GDtmp, +1)
! 步骤3：G <- G * B_gauge^{-1}（右乘 B_gauge 的逆）
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'y', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 2, 'x', -1, 0.5d0)
                call apply_group_L(GUtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
                call apply_group_L(GDtmp, Latt, Bonds, Gauge, t, 1, 'x', -1, 0.5d0)
! 步骤4：G <- G * exp(-μ)
                call opMu_mmult_L(GUtmp, -1); call opMu_mmult_L(GDtmp, -1)
            endif
        enddo

        dS_tot = dS_spa + dS_time + dS_bdry

        thr = ranf(iseed)
! 接受判据：logA = log(Rf_tot) - dS_tot，与 log(thr) 比较，避免下溢
        logA = log_Rf_tot - dS_tot
        if (IRANK == 0) then
            open(unit=96, file='worm_debug.log', status='unknown', action='write', position='append')
            write(96,'(A,1X,I6,1X,I6,1X,I6,1X,ES13.6,1X,ES13.6,1X,ES13.6,1X,ES13.6)') &
                'time_string', nt0, tau1, tau2, log_Rf_tot, dS_spa, dS_time + dS_bdry, logA
            close(96)
        endif
        if (logA > log(max(thr, 1.d-300))) then
            acc = .true.
            GrU = GUtmp; GrD = GDtmp
        else
! 回滚 σ（再翻一次区间内的 σ 恢复初值）
            do t = tau1, tau2
                if (bdir == 'x') then
                    Gauge%sigma_x(i_src, t) = - Gauge%sigma_x(i_src, t)
                else
                    Gauge%sigma_y(i_src, t) = - Gauge%sigma_y(i_src, t)
                endif
            enddo
        endif
        deallocate(GUtmp); deallocate(GDtmp)
        return
    end subroutine worm_time_string

    subroutine worm_small_loop(GrU, GrD, Gauge, Latt, nt, x, y, iseed, acc)
! 在时片 nt 上翻转以 (x,y) 为左下角的 plaquette 四条边（最小环）
        complex(kind=8), intent(inout) :: GrU(:, :), GrD(:, :)
        class(GaugeConf), intent(inout) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: nt, x, y
        integer, intent(inout) :: iseed
        logical, intent(out) :: acc
        integer :: i_ll, i_lr, i_ul
        integer :: xp1, yp1
        integer :: t
        real(kind=8) :: dS_tot, log_Rf_tot, r_up, r_dn, thr, val, logA
        include 'mpif.h'
        real(kind=8), external :: ranf
        complex(kind=8), allocatable :: GUtmp(:, :), GDtmp(:, :)

        acc = .false.
        
! 特殊处理：避免在最后一片更新（数值不稳定）
        if (nt >= Ltrot) then
            return
        endif
        xp1 = npbc(x+1, Nlx); yp1 = npbc(y+1, Nly)
        i_ll = Latt%inv_n_list(x,   y)
        i_lr = Latt%inv_n_list(xp1, y)
        i_ul = Latt%inv_n_list(x,   yp1)

        allocate(GUtmp(Ndim, Ndim)); allocate(GDtmp(Ndim, Ndim))
        GUtmp = GrU; GDtmp = GrD
        log_Rf_tot = 0.d0

! 计算同时翻转4条边的作用量变化（使用专门的函数）
        dS_tot = delta_S_plaquette_flip(Gauge, Latt, x, y, nt)

! 四条边依次翻转：σx(x,y), σy(x+1,y), σx(x,y+1), σy(x,y)
        ! 先测试det_ratio，使用更安全的阈值
        call sigma_flip_rank2(GUtmp, i_ll, i_lr, Gauge%sigma_x(i_ll, nt), r_up, .false.)
        call sigma_flip_rank2(GDtmp, i_ll, i_lr, Gauge%sigma_x(i_ll, nt), r_dn, .false.)
        if (r_up <= 1.d-100 .or. r_dn <= 1.d-100) then
            ! 物理错误或数值不稳定，直接拒绝
            acc = .false.
            deallocate(GUtmp); deallocate(GDtmp)
            return
        endif
        ! det_ratio可接受，执行更新
        call sigma_flip_rank2(GUtmp, i_ll, i_lr, Gauge%sigma_x(i_ll, nt), r_up, .true.)
        call sigma_flip_rank2(GDtmp, i_ll, i_lr, Gauge%sigma_x(i_ll, nt), r_dn, .true.)
        val = r_up * r_dn
        log_Rf_tot = log_Rf_tot + log(max(val, 1.d-100))
        Gauge%sigma_x(i_ll, nt) = - Gauge%sigma_x(i_ll, nt)

        ! 先测试第二条边
        call sigma_flip_rank2(GUtmp, i_lr, Latt%L_bonds(i_lr, 2), Gauge%sigma_y(i_lr, nt), r_up, .false.)
        call sigma_flip_rank2(GDtmp, i_lr, Latt%L_bonds(i_lr, 2), Gauge%sigma_y(i_lr, nt), r_dn, .false.)
        if (r_up <= 1.d-100 .or. r_dn <= 1.d-100) then
            ! 回滚第一条边并拒绝
            Gauge%sigma_x(i_ll, nt) = - Gauge%sigma_x(i_ll, nt)
            acc = .false.
            deallocate(GUtmp); deallocate(GDtmp)
            return
        endif
        ! 执行更新
        call sigma_flip_rank2(GUtmp, i_lr, Latt%L_bonds(i_lr, 2), Gauge%sigma_y(i_lr, nt), r_up, .true.)
        call sigma_flip_rank2(GDtmp, i_lr, Latt%L_bonds(i_lr, 2), Gauge%sigma_y(i_lr, nt), r_dn, .true.)
        val = r_up * r_dn
        log_Rf_tot = log_Rf_tot + log(max(val, 1.d-100))
        Gauge%sigma_y(i_lr, nt) = - Gauge%sigma_y(i_lr, nt)

        ! 先测试第三条边
        call sigma_flip_rank2(GUtmp, i_ul, Latt%L_bonds(i_ul, 1), Gauge%sigma_x(i_ul, nt), r_up, .false.)
        call sigma_flip_rank2(GDtmp, i_ul, Latt%L_bonds(i_ul, 1), Gauge%sigma_x(i_ul, nt), r_dn, .false.)
        if (r_up <= 1.d-100 .or. r_dn <= 1.d-100) then
            ! 回滚前两条边并拒绝
            Gauge%sigma_y(i_lr, nt) = - Gauge%sigma_y(i_lr, nt)
            Gauge%sigma_x(i_ll, nt) = - Gauge%sigma_x(i_ll, nt)
            acc = .false.
            deallocate(GUtmp); deallocate(GDtmp)
            return
        endif
        ! 执行更新
        call sigma_flip_rank2(GUtmp, i_ul, Latt%L_bonds(i_ul, 1), Gauge%sigma_x(i_ul, nt), r_up, .true.)
        call sigma_flip_rank2(GDtmp, i_ul, Latt%L_bonds(i_ul, 1), Gauge%sigma_x(i_ul, nt), r_dn, .true.)
        val = r_up * r_dn
        log_Rf_tot = log_Rf_tot + log(max(val, 1.d-100))
        Gauge%sigma_x(i_ul, nt) = - Gauge%sigma_x(i_ul, nt)

        ! 先测试第四条边
        call sigma_flip_rank2(GUtmp, i_ll, Latt%L_bonds(i_ll, 2), Gauge%sigma_y(i_ll, nt), r_up, .false.)
        call sigma_flip_rank2(GDtmp, i_ll, Latt%L_bonds(i_ll, 2), Gauge%sigma_y(i_ll, nt), r_dn, .false.)
        if (r_up <= 1.d-100 .or. r_dn <= 1.d-100) then
            ! 回滚前三条边并拒绝
            Gauge%sigma_x(i_ul, nt) = - Gauge%sigma_x(i_ul, nt)
            Gauge%sigma_y(i_lr, nt) = - Gauge%sigma_y(i_lr, nt)
            Gauge%sigma_x(i_ll, nt) = - Gauge%sigma_x(i_ll, nt)
            acc = .false.
            deallocate(GUtmp); deallocate(GDtmp)
            return
        endif
        ! 执行更新
        call sigma_flip_rank2(GUtmp, i_ll, Latt%L_bonds(i_ll, 2), Gauge%sigma_y(i_ll, nt), r_up, .true.)
        call sigma_flip_rank2(GDtmp, i_ll, Latt%L_bonds(i_ll, 2), Gauge%sigma_y(i_ll, nt), r_dn, .true.)
        val = r_up * r_dn
        log_Rf_tot = log_Rf_tot + log(max(val, 1.d-100))
        Gauge%sigma_y(i_ll, nt) = - Gauge%sigma_y(i_ll, nt)

        thr = ranf(iseed)
        logA = log_Rf_tot - dS_tot
        ! 调试输出：检查自旋简并性（暂时关闭）
        ! if (IRANK == 0 .and. abs(r_up - r_dn) > 1.d-10) then
        !     open(unit=97, file='spin_diff.log', status='unknown', action='write', position='append')
        !     write(97,*) 'small_loop: r_up=', r_up, ' r_dn=', r_dn, ' diff=', abs(r_up - r_dn)
        !     close(97)
        ! endif
        if (IRANK == 0) then
            open(unit=96, file='worm_debug.log', status='unknown', action='write', position='append')
            write(96,'(A,1X,I6,1X,I6,1X,I6,1X,ES13.6,1X,ES13.6,1X,ES13.6)') 'small_loop', nt, x, y, log_Rf_tot, dS_tot, logA
            close(96)
        endif
        if (logA > log(max(thr, 1.d-300))) then
            acc = .true.
            GrU = GUtmp; GrD = GDtmp
        else
! 回滚四条边
            Gauge%sigma_y(i_ll, nt) = - Gauge%sigma_y(i_ll, nt)
            Gauge%sigma_x(i_ul, nt) = - Gauge%sigma_x(i_ul, nt)
            Gauge%sigma_y(i_lr, nt) = - Gauge%sigma_y(i_lr, nt)
            Gauge%sigma_x(i_ll, nt) = - Gauge%sigma_x(i_ll, nt)
        endif
        deallocate(GUtmp); deallocate(GDtmp)
        return
    end subroutine worm_small_loop

end module Worm_mod

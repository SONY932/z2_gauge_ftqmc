module LocalSweepV1_mod
    use calc_basic
    use ProcessMatrix
    use Stabilize_mod
    use MyLattice
    use BondList_mod
    use GaugeFields_mod
    use LocalSigma_mod
    use MultiplyChecker_mod
    use OpMu_mod
    use SigmaRank2_mod
    use LambdaPair_mod
    use GaugeAction_mod
    use ObserMin_mod
    use Worm_mod
    use LocalUpdateOps_mod
    implicit none
    logical, parameter :: debug_green_log_enabled = .false.
    logical, parameter :: debug_sigma_flip_enabled = .false.
    integer, save :: debug_sigma_flip_count = 0
    integer, parameter :: debug_sigma_flip_target_nt = 6
    integer, parameter :: debug_sigma_flip_limit = 128

    type, public :: LocalSweep
        type(AccCounter) :: Acc_sigma
        type(AccCounter) :: Acc_lambda
        type(AccCounter) :: Acc_worm
        type(ObserMin) :: Obs
    contains
        procedure :: init  => local_init
        procedure :: reset => local_reset
        procedure :: sweep => local_sweep
        procedure :: pre   => local_pre
    end type LocalSweep

contains
    subroutine local_init(this)
        class(LocalSweep), intent(inout) :: this
        call this%Acc_sigma%init()
        call this%Acc_lambda%init()
        call this%Acc_worm%init()
        call this%Obs%reset()
        return
    end subroutine local_init

    subroutine local_reset(this)
        class(LocalSweep), intent(inout) :: this
        call this%Acc_sigma%reset()
        call this%Acc_lambda%reset()
        call this%Acc_worm%reset()
        return
    end subroutine local_reset

    subroutine local_pre(this, PropU, PropD, WrU, WrD, Latt, Bonds, Gauge)
! 按二阶对称顺序构建 UUR 段栈并 Wrap_pre
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer :: nt
        real(kind=8), parameter :: half = 0.5d0

        call this%reset()
        WrU%ULlist = dcmplx(0.d0, 0.d0); WrU%VLlist = dcmplx(0.d0, 0.d0); WrU%DLlist = dcmplx(0.d0, 0.d0)
        WrD%ULlist = dcmplx(0.d0, 0.d0); WrD%VLlist = dcmplx(0.d0, 0.d0); WrD%DLlist = dcmplx(0.d0, 0.d0)
        call reset_debug_green()
        call reset_debug_sigma_flip()
        call reset_debug_wrap_detail()
        call reset_debug_operator()
        call Stabilize_init()
        call Wrap_pre(PropU, WrU, 0)
        call Wrap_pre(PropD, WrD, 0)
        do nt = 1, Ltrot
! UUR ← UUR * B(nt)（右乘保持与 AFM 框架一致）
            call prop_pre_step(PropU, PropD, Latt, Bonds, Gauge, nt)
            if (mod(nt, Nwrap) == 0) then
                call Wrap_pre(PropU, WrU, nt)
                call Wrap_pre(PropD, WrD, nt)
            endif
        enddo
        call Wrap_pre(PropU, WrU, Ltrot)
        call Wrap_pre(PropD, WrD, Ltrot)
        return
    end subroutine local_pre

    subroutine propagate_left_step(this, PropU, PropD, Latt, Bonds, Gauge, nt, iseed)
! 左扫完整传播：G(nt-1) = B^{-1} * G(nt) * B
! 其中 B = B_gauge * e^{-μ}，B^{-1} = e^{μ} * B_gauge^{-1}
! 所以 G' = e^{μ} * B_gauge^{-1} * G * B_gauge * e^{-μ}
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(inout) :: Gauge
        integer, intent(in) :: nt
        integer, intent(inout) :: iseed

! σ Metropolis：

! 完整传播：G' = e^{μ} * B_gauge^{-1} * G * B_gauge * e^{-μ}
! 按从右到左的顺序：
! 步骤1：G <- G * B_gauge（正向右乘）
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
! 步骤2：G <- G * e^{-μ}
        call opMu_mmult_R(PropU%Gr, -1); call opMu_mmult_R(PropD%Gr, -1)
! 步骤3：G <- B_gauge^{-1} * G（逆向左乘）
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤4：G <- e^{μ} * G
        call opMu_mmult_L(PropU%Gr, +1); call opMu_mmult_L(PropD%Gr, +1)

! 累积 UUL
        call left_step_after_sigma_post_mu(PropU, PropD, Latt, Bonds, Gauge, nt)
        return
    end subroutine propagate_left_step

    subroutine propagate_right_step(this, PropU, PropD, Latt, Bonds, Gauge, nt, iseed)
! 右扫完整传播：G(nt+1) = B * G * B^{-1}
! 其中 B = B_gauge * e^{-μ}，B^{-1} = e^{μ} * B_gauge^{-1}
! 所以 G' = B_gauge * e^{-μ} * G * e^{μ} * B_gauge^{-1}
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(inout) :: Gauge
        integer, intent(in) :: nt
        integer, intent(inout) :: iseed

! σ Metropolis：

! 完整传播：G' = B_gauge * e^{-μ} * G * e^{μ} * B_gauge^{-1}
! 按从右到左的顺序：
! 步骤1：G <- G * B_gauge^{-1}（逆向右乘）
        call apply_trotter_layer_R(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
        call apply_trotter_layer_R(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=-1)
! 步骤2：G <- G * e^{μ}
        call opMu_mmult_R(PropU%Gr, +1); call opMu_mmult_R(PropD%Gr, +1)
! 步骤3：G <- e^{-μ} * G
        call opMu_mmult_L(PropU%Gr, -1); call opMu_mmult_L(PropD%Gr, -1)
! 步骤4：G <- B_gauge * G（正向左乘）
        call apply_trotter_layer_L(PropU%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)
        call apply_trotter_layer_L(PropD%Gr, Latt, Bonds, Gauge, nt, nflag_in=+1)

! 累积 UUR: UUR <- UUR * B = UUR * B_gauge * e^{-μ}
        call right_step_finalize(PropU, PropD, Latt, Bonds, Gauge, nt)
        call opMu_mmult_R(PropU%UUR, -1); call opMu_mmult_R(PropD%UUR, -1)
        return
    end subroutine propagate_right_step

    subroutine local_sweep(this, PropU, PropD, WrU, WrD, Latt, Bonds, Gauge, iseed)
! 左扫在前、右扫在后，完全对齐 AFM-ISING CODE 的段栈读写顺序
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(inout) :: Gauge
        integer, intent(inout) :: iseed
        integer :: nt
        integer :: ng, x, y
        integer :: i_lambda, j_lambda
        logical :: acc, did_worm, did_lambda, acc_lambda
        character(len=1) :: bdir
        real(kind=8), parameter :: half = 0.5d0
        real(kind=8), external :: ranf

        did_worm = .false.
        did_lambda = .false.
        acc = .false.
        acc_lambda = .false.

        call this%reset()

! 左扫（Ltrot..1），先恢复左向缓存再传播
        do nt = Ltrot, 1, -1
            if (mod(nt, Nwrap) == 0) then
                call Wrap_L(PropU, WrU, nt)
                call Wrap_L(PropD, WrD, nt)
            endif
            call propagate_left_step(this, PropU, PropD, Latt, Bonds, Gauge, nt, iseed)
        enddo
! 左扫结束后，在nt=0处做最后一次Wrap_L
        call Wrap_L(PropU, WrU, 0)
        call Wrap_L(PropD, WrD, 0)

! 右扫（1..Ltrot），左向缓存已齐全
        call Wrap_R(PropU, WrU, 0)
        call Wrap_R(PropD, WrD, 0)
        do nt = 1, Ltrot
            call propagate_right_step(this, PropU, PropD, Latt, Bonds, Gauge, nt, iseed)
            
            if (mod(nt, Nwrap) == 0) then
                call Wrap_R(PropU, WrU, nt)
                call Wrap_R(PropD, WrD, nt)
            endif
        enddo

! 可选：worm 全局更新（时间串 / 最小环），每 sweep 做 Nglobal 次提案
        if (.false.) then  ! DEBUG: 暂时禁用worm更新来测试传播
            do ng = 1, Nglobal
! 在 Nwrap 的时片边界重建 G，然后随机选择：时间串或最小环
                nt = nranf(iseed, Ltrot / Nwrap) * Nwrap
                if (nt == 0) nt = Nwrap
                call Wrap_R(PropU, WrU, nt)
                call Wrap_R(PropD, WrD, nt)
                if (ranf(iseed) < 0.5d0) then
! 时间串：区间限制在该边界段内，worm 里按时片推进 G
                    if (ranf(iseed) < 0.5d0) then
                        bdir = 'x'
                    else
                        bdir = 'y'
                    endif
                    call worm_time_string(PropU%Gr, PropD%Gr, Gauge, Latt, Bonds, nt, bdir, iseed, acc)
                else
! 最小环：同一边界时片
                    x  = nranf(iseed, Nlx)
                    y  = nranf(iseed, Nly)
                    call worm_small_loop(PropU%Gr, PropD%Gr, Gauge, Latt, nt, x, y, iseed, acc)
                endif
                call this%Acc_worm%count(acc)
                if (acc) then
                    did_worm = .true.
                    exit
                endif
            enddo
        endif

! ====== λ 更新已禁用（仅保留σ更新用于调试） ======
! ! 施加 λ 投影后再进行 λ 成对翻与观测
!         call apply_lambda_both(PropU%Gr, Gauge)
!         call apply_lambda_both(PropD%Gr, Gauge)
! 
!         i_lambda = nranf(iseed, Lq)
!         j_lambda = nranf(iseed, Lq)
!         if (i_lambda /= j_lambda) then
!             call lambda_pair_flip(PropU%Gr, PropD%Gr, Gauge, i_lambda, j_lambda, Latt, iseed, acc_lambda)
!             call this%Acc_lambda%count(acc_lambda)
!             if (acc_lambda) did_lambda = .true.
!         endif

! 末片观测（无λ投影时直接观测）
        call this%Obs%acc_fermion(PropU%Gr, PropD%Gr, Latt)
        call this%Obs%acc_gauge(Gauge, Latt, Ltrot)
        call this%Obs%reduce_and_write()

! ! 撤回 λ 投影，恢复传播基
!         call apply_lambda_both(PropU%Gr, Gauge)
!         call apply_lambda_both(PropD%Gr, Gauge)
! 汇总并输出接受率（MPI 均值）
        call write_accept_logs(this)

! 若 worm 或 λ 接受，都必须重建段栈，确保与新的规范场一致（仿 AFM 逻辑）
        if (did_worm .or. did_lambda) then
            call this%pre(PropU, PropD, WrU, WrD, Latt, Bonds, Gauge)
        endif
        return
    end subroutine local_sweep

    subroutine sweep_sigma_dir(this, GrU, GrD, Gauge, Latt, gIdx, dir, nt, srcList, dstList, iseed, tag)
! 遍历一组键，逐键提案 rank-2 更新
        class(LocalSweep), intent(inout) :: this
        complex(kind=8), intent(inout) :: GrU(:, :), GrD(:, :)
        class(GaugeConf), intent(inout) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: gIdx(:)
        character(len=*), intent(in) :: dir
        integer, intent(in) :: nt
        integer, intent(in) :: srcList(:), dstList(:)
        integer, intent(inout) :: iseed
        character(len=*), intent(in) :: tag
        integer :: k, bidx, i, j
        integer :: sigma_val
        real(kind=8) :: dS, r_up, r_dn, Rf, thr
        real(kind=8), external :: ranf
        logical :: accepted

        do k = 1, size(gIdx)
            bidx = gIdx(k)
            i = srcList(bidx); j = dstList(bidx)
            dS = delta_S_sigma_flip(Gauge, Latt, dir, i, nt)
            if (dir == 'x') then
                sigma_val = Gauge%sigma_x(i, nt)
                call debug_log_sigma_flip('pre', tag, dir, nt, k, i, j, sigma_val, GrU, GrD)
                call sigma_flip_rank2(GrU, i, j, sigma_val, r_up, .false.)
                call sigma_flip_rank2(GrD, i, j, sigma_val, r_dn, .false.)
                r_up = max(r_up, 0.d0)
                r_dn = max(r_dn, 0.d0)
                Rf = max(r_up, 1.d-300) * max(r_dn, 1.d-300)
                thr = ranf(iseed)
                accepted = (min(1.d0, Rf * exp(-dS)) > thr)
                if (accepted) then
                    call this%Acc_sigma%count(.true.)
                    call sigma_flip_rank2(GrU, i, j, sigma_val, r_up, .true.)
                    call sigma_flip_rank2(GrD, i, j, sigma_val, r_dn, .true.)
                    Gauge%sigma_x(i, nt) = - Gauge%sigma_x(i, nt)
                    call debug_log_sigma_flip('post_accept', tag, dir, nt, k, i, j, Gauge%sigma_x(i, nt), GrU, GrD)
                else
                    call this%Acc_sigma%count(.false.)
                    call debug_log_sigma_flip('post_reject', tag, dir, nt, k, i, j, sigma_val, GrU, GrD)
                endif
            else
                sigma_val = Gauge%sigma_y(i, nt)
                call debug_log_sigma_flip('pre', tag, dir, nt, k, i, j, sigma_val, GrU, GrD)
                call sigma_flip_rank2(GrU, i, j, sigma_val, r_up, .false.)
                call sigma_flip_rank2(GrD, i, j, sigma_val, r_dn, .false.)
                r_up = max(r_up, 0.d0)
                r_dn = max(r_dn, 0.d0)
                Rf = max(r_up, 1.d-300) * max(r_dn, 1.d-300)
                thr = ranf(iseed)
                accepted = (min(1.d0, Rf * exp(-dS)) > thr)
                if (accepted) then
                    call this%Acc_sigma%count(.true.)
                    call sigma_flip_rank2(GrU, i, j, sigma_val, r_up, .true.)
                    call sigma_flip_rank2(GrD, i, j, sigma_val, r_dn, .true.)
                    Gauge%sigma_y(i, nt) = - Gauge%sigma_y(i, nt)
                    call debug_log_sigma_flip('post_accept', tag, dir, nt, k, i, j, Gauge%sigma_y(i, nt), GrU, GrD)
                else
                    call this%Acc_sigma%count(.false.)
                    call debug_log_sigma_flip('post_reject', tag, dir, nt, k, i, j, sigma_val, GrU, GrD)
                endif
            endif
        enddo
        return
    end subroutine sweep_sigma_dir

    subroutine write_accept_logs(this)
        class(LocalSweep), intent(inout) :: this
        include 'mpif.h'
        real(kind=8) :: sig_all, lam_all, worm_all
        call this%Acc_sigma%ratio()
        call this%Acc_lambda%ratio()
        call this%Acc_worm%ratio()
        sig_all = 0.d0; lam_all = 0.d0; worm_all = 0.d0
        call MPI_Reduce(this%Acc_sigma%acc, sig_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        call MPI_Reduce(this%Acc_lambda%acc, lam_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        call MPI_Reduce(this%Acc_worm%acc, worm_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) then
            sig_all = sig_all / dble(ISIZE)
            lam_all = lam_all / dble(ISIZE)
            worm_all = worm_all / dble(ISIZE)
            open(unit=84, file='accept_sigma.log', status='unknown', action='write', position='append')
            write(84,*) sig_all
            call flush(84)
            close(84)
            open(unit=85, file='accept_lambda.log', status='unknown', action='write', position='append')
            write(85,*) lam_all
            call flush(85)
            close(85)
            open(unit=86, file='accept_worm.log', status='unknown', action='write', position='append')
            write(86,*) worm_all
            call flush(86)
            close(86)
        endif
        return
    end subroutine write_accept_logs

    subroutine debug_log_green(tag, nt, GrU, GrD)
        implicit none
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        complex(kind=8), intent(in) :: GrU(:, :), GrD(:, :)
        real(kind=8) :: diag_min, diag_max, diag_avg
        real(kind=8) :: diag_min_dn, diag_max_dn, diag_avg_dn
        real(kind=8) :: max_off
        integer :: i, j, n
        integer, parameter :: limit = 400
        integer, save :: counter = 0

        n = size(GrU,1)
        if (n <= 0) return

        diag_min = real(GrU(1,1)); diag_max = real(GrU(1,1)); diag_avg = 0.d0
        diag_min_dn = real(GrD(1,1)); diag_max_dn = real(GrD(1,1)); diag_avg_dn = 0.d0
        max_off = 0.d0

        do i = 1, n
            if (.not.(real(GrU(i,i)) == real(GrU(i,i))) .or. .not.(aimag(GrU(i,i)) == aimag(GrU(i,i)))) then
                call debug_log_nan('GrU', tag, nt, i, 0, GrU(i,i))
                return
            endif
            if (.not.(real(GrD(i,i)) == real(GrD(i,i))) .or. .not.(aimag(GrD(i,i)) == aimag(GrD(i,i)))) then
                call debug_log_nan('GrD', tag, nt, i, 0, GrD(i,i))
                return
            endif
            diag_min = min(diag_min, real(GrU(i,i)))
            diag_max = max(diag_max, real(GrU(i,i)))
            diag_avg = diag_avg + real(GrU(i,i))
            diag_min_dn = min(diag_min_dn, real(GrD(i,i)))
            diag_max_dn = max(diag_max_dn, real(GrD(i,i)))
            diag_avg_dn = diag_avg_dn + real(GrD(i,i))
        enddo
        diag_avg    = diag_avg    / dble(n)
        diag_avg_dn = diag_avg_dn / dble(n)

        do i = 1, n
            do j = 1, size(GrU,2)
                if (i == j) cycle
                if (.not.(real(GrU(i,j)) == real(GrU(i,j))) .or. .not.(aimag(GrU(i,j)) == aimag(GrU(i,j)))) then
                    call debug_log_nan('GrU', tag, nt, i, j, GrU(i,j))
                    return
                endif
                if (.not.(real(GrD(i,j)) == real(GrD(i,j))) .or. .not.(aimag(GrD(i,j)) == aimag(GrD(i,j)))) then
                    call debug_log_nan('GrD', tag, nt, i, j, GrD(i,j))
                    return
                endif
                max_off = max(max_off, abs(GrU(i,j)), abs(GrD(i,j)))
            enddo
        enddo

        if (debug_green_log_enabled) then
            if (counter < limit) then
                counter = counter + 1
                open(unit=209, file='debug_local_sweep.log', status='unknown', action='write', position='append')
                write(209,'(I8,1X,A12,1X,I8,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6)') &
                    counter, tag, nt, diag_min, diag_max, diag_avg, diag_min_dn, diag_max_dn, diag_avg_dn, max_off
                close(209)
            endif
        endif
        return
    end subroutine debug_log_green

    subroutine reset_debug_green()
        implicit none
        if (.not. debug_green_log_enabled) return
        open(unit=209, file='debug_local_sweep.log', status='unknown', action='write')
        close(209)
        return
    end subroutine reset_debug_green

    subroutine reset_debug_sigma_flip()
        implicit none
        debug_sigma_flip_count = 0
        if (.not. debug_sigma_flip_enabled) return
        open(unit=213, file='debug_sigma_flip.log', status='unknown', action='write')
        close(213)
        return
    end subroutine reset_debug_sigma_flip

    subroutine debug_log_sigma_flip(stage, tag, dir, nt, k, i, j, sigma, GrU, GrD)
        implicit none
        character(len=*), intent(in) :: stage, tag, dir
        integer, intent(in) :: nt, k, i, j, sigma
        complex(kind=8), intent(in) :: GrU(:, :), GrD(:, :)
        complex(kind=8) :: uu_ii, uu_jj, uu_ij, uu_ji
        complex(kind=8) :: dd_ii, dd_jj, dd_ij, dd_ji
        real(kind=8) :: max_u, max_d
        logical :: nan_u, nan_d

        if (.not. debug_sigma_flip_enabled) return
        if (nt /= debug_sigma_flip_target_nt) return
        if (debug_sigma_flip_count >= debug_sigma_flip_limit) return

        uu_ii = GrU(i, i); uu_jj = GrU(j, j)
        uu_ij = GrU(i, j); uu_ji = GrU(j, i)
        dd_ii = GrD(i, i); dd_jj = GrD(j, j)
        dd_ij = GrD(i, j); dd_ji = GrD(j, i)

        nan_u = .false.
        nan_d = .false.

        if (.not.(real(uu_ii) == real(uu_ii)) .or. .not.(aimag(uu_ii) == aimag(uu_ii))) nan_u = .true.
        if (.not.(real(uu_jj) == real(uu_jj)) .or. .not.(aimag(uu_jj) == aimag(uu_jj))) nan_u = .true.
        if (.not.(real(uu_ij) == real(uu_ij)) .or. .not.(aimag(uu_ij) == aimag(uu_ij))) nan_u = .true.
        if (.not.(real(uu_ji) == real(uu_ji)) .or. .not.(aimag(uu_ji) == aimag(uu_ji))) nan_u = .true.

        if (.not.(real(dd_ii) == real(dd_ii)) .or. .not.(aimag(dd_ii) == aimag(dd_ii))) nan_d = .true.
        if (.not.(real(dd_jj) == real(dd_jj)) .or. .not.(aimag(dd_jj) == aimag(dd_jj))) nan_d = .true.
        if (.not.(real(dd_ij) == real(dd_ij)) .or. .not.(aimag(dd_ij) == aimag(dd_ij))) nan_d = .true.
        if (.not.(real(dd_ji) == real(dd_ji)) .or. .not.(aimag(dd_ji) == aimag(dd_ji))) nan_d = .true.

        max_u = max(abs(uu_ii), abs(uu_jj))
        max_u = max(max_u, abs(uu_ij))
        max_u = max(max_u, abs(uu_ji))

        max_d = max(abs(dd_ii), abs(dd_jj))
        max_d = max(max_d, abs(dd_ij))
        max_d = max(max_d, abs(dd_ji))

        debug_sigma_flip_count = debug_sigma_flip_count + 1
        open(unit=213, file='debug_sigma_flip.log', status='unknown', action='write', position='append')
        write(213,*) '#', debug_sigma_flip_count, trim(stage), trim(tag), trim(dir), 'nt', nt, 'k', k, 'i', i, 'j', j, 'sigma', sigma
        write(213,*) 'GrU_block', uu_ii, uu_jj, uu_ij, uu_ji, 'has_nan', nan_u, 'max_abs', max_u
        write(213,*) 'GrD_block', dd_ii, dd_jj, dd_ij, dd_ji, 'has_nan', nan_d, 'max_abs', max_d
        close(213)

        if (nan_u .or. nan_d) stop 'debug_sigma_flip detected NaN'

        return
    end subroutine debug_log_sigma_flip

    subroutine debug_log_nan(spin_tag, step_tag, nt, i, j, val)
        implicit none
        character(len=*), intent(in) :: spin_tag
        character(len=*), intent(in) :: step_tag
        integer, intent(in) :: nt, i, j
        complex(kind=8), intent(in) :: val
        open(unit=212, file='debug_nan.log', status='unknown', action='write', position='append')
        write(212,'(A4,1X,A12,1X,I8,1X,I6,1X,I6,1X,ES14.6,1X,ES14.6)') &
            spin_tag, step_tag, nt, i, j, real(val), aimag(val)
        close(212)
        return
    end subroutine debug_log_nan

end module LocalSweepV1_mod

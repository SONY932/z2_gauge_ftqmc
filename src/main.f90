program DQMC_GaugeMain
    use calc_basic
    use MyLattice
    use BondList_mod
    use GaugeFields_mod
    use ProcessMatrix
    use Stabilize_mod
    use LocalSweepV1_mod
    implicit none

    integer :: nbc
    type(SquareLattice) :: Latt
    type(BondList) :: Bonds
    type(GaugeConf) :: Gauge
    type(Propagator) :: PropU, PropD
    type(WrapList) :: WrU, WrD
    type(LocalSweep) :: Sweep
    integer :: iseed
    integer(kind=8) :: clock_now, clock_rate
    logical :: conf_loaded

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)

! 读参与初始化
    call read_input()
    call Params_set()
    call write_info()
    call Lattice_make(Latt)
    call Bonds%make(Latt)
    call SYSTEM_CLOCK(COUNT=clock_now, COUNT_RATE=clock_rate)
    if (clock_rate > 0_8) then
        iseed = int(mod(clock_now + int(IRANK, kind=8) + 1_8, 2147483646_8)) + 1
    else
        iseed = 13579 + IRANK
    endif
    if (iseed <= 0) iseed = abs(iseed) + 1
    call Gauge%make(iseed, Latt)
    call Gauge_conf_in(Gauge, iseed, conf_loaded)
    if (.not. conf_loaded) then
        if (IRANK == 0) write(6,*) 'Gauge_conf_in: 未检测到 confin/seeds，采用随机初始化'
    endif

! 稳定化/传播器/预处理
    call PropU%make(); call PropD%make()
    call WrU%make();   call WrD%make()
    call Sweep%init()

! 简单 bin 循环
    do nbc = 1, Nbin
        ! AFM-ISING 每个 sweep 前都会重建段栈；这里仿照其流程在每个 bin 开始时调用 pre，
        ! 保证 Wrap_L/Wrap_R 恢复时段栈是最新的（避免读取被清空的缓存）。
        call Sweep%pre(PropU, PropD, WrU, WrD, Latt, Bonds, Gauge)
        call Sweep%sweep(PropU, PropD, WrU, WrD, Latt, Bonds, Gauge, iseed)
    enddo

! 清理
    call Gauge_conf_out(Gauge, iseed)
    call Lattice_clear(Latt)
    call Stabilize_clear()
    call MPI_FINALIZE(IERR)
end program DQMC_GaugeMain



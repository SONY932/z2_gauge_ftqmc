module calc_basic
    implicit none

! 常数
    real(kind = 8), parameter :: Zero = 1.0d-10
    real(kind = 8), parameter :: PI = acos(-1.d0)
    real(kind = 8), parameter :: upbound = 1.0d+200
! 晶格参数
    integer, parameter :: Norb = 2 ! 该模型有两个轨道，对应论文里面的λ = 1、2
    integer, parameter :: Nbond = 2 ! 晶格有两个方向，水平和竖直
    integer, parameter :: Nspin = 2 ! 自旋σ = 1、2
    integer, public :: Nlx, Nly, NlxTherm, NlyTherm ! 晶格大小和热化时的晶格大小
    integer, public :: Lq, LqTherm ! 总格点数目
    integer, public :: Ndim ! 格林函数维度
    real(kind = 8), public :: Dtau ! 虚时间隔
    real(kind = 8), public :: Beta ! 逆温度
    integer, public :: Ltrot, LtrotTherm ! 虚时长度
! 哈密顿量参数
    real(kind = 8), public :: RT ! 跃迁项的系数 t
    real(kind = 8), public :: J  ! 空间规范耦合强度（沿用 J）
    real(kind = 8), public :: h  ! 横场（用于 γ = (1/2) ln coth(Δτ h)）
    real(kind = 8), public :: xi ! Ising场耦合强度（旧参数，后续不再直接使用）
    real(kind = 8), public :: mu ! 化学势
! 扭转载边界角已不使用（保持占位为 0）
    real(kind = 8), public :: ThetaX, ThetaY
! 新增：扇区选择（末片 λ 的乘积约束）
    integer, public :: lambda_sector ! +1 或 -1
! 更新参数
    logical, public :: is_global ! 是否全局更新
    integer, public :: Nglobal ! 全局更新频率（多少步更新一次）
! 初始化状态参数（NB_field 不再使用，置 0）
    real(kind = 8), public :: NB_field
    integer, public :: absolute ! 选择初始化Ising场的方式（沿用旧接口）
! 过程控制参数
    logical, public :: is_tau ! 是否计算含时观测量
    integer, public :: Nthermal ! 计算含时观测量前热化步数（自定义）
    logical, public :: is_warm ! 是否进行热化（只使用Ising作用量更新场）
    integer, public :: Nwarm ! 热化步数（自定义）
    integer, public :: Nst ! 存储稳定矩阵的格式，Ltrot/Nwrap
    integer, public :: Nwrap ! 稳定化频率，隔 Nwrap 步存储一次
    integer, public :: Nwrap_input = 0
    integer, public :: Nbin ! 数据存储频率，隔 Nbin 步存储一次
    integer, public :: Nsweep ! 来回扫的更新次数（蒙卡步）
    integer, public :: ISIZE, IRANK, IERR ! 并行参数
! 守卫与修正选项
    logical, public :: dtau_guard
    logical, public :: zero_fix
    logical, public :: enable_reorth = .true.
    integer, public :: ortho_log_stride = 10
    integer, public :: nwrap_divider = 1
! 稳定性/全局更新附加参数
    real(kind = 8), public :: ortho_threshold = 5.5d-5
    integer, public :: worm_max_span = 4

contains
    subroutine read_input()
        include 'mpif.h'
        integer :: ios
        if (IRANK == 0) then
            open(unit = 20, file = 'paramC_sets.txt', status = 'unknown')
            read(20, *) RT, J, h, xi, mu
            read(20, *) Nlx, Nly, Ltrot, Beta
            read(20, *) NlxTherm, NlyTherm, LtrotTherm
            read(20, *) Nwrap, Nbin, Nsweep
            Nwrap_input = Nwrap
            read(20, *) is_tau, Nthermal
            read(20, *) is_warm, Nwarm
            read(20, *) is_global, Nglobal
            read(20, *) absolute
            NB_field = 0.d0
            ThetaX = 0.d0; ThetaY = 0.d0
! 新增：lambda 扇区（可选行，缺省取 +1）
            ios = 0
            lambda_sector = 1
            read(20, *, iostat=ios) lambda_sector
            if (ios /= 0) then
                lambda_sector = 1
                ios = 0
            endif
! 新增：守卫与扩展修正（可选行，缺省关闭）
            dtau_guard = .false.
            zero_fix   = .false.
            read(20, *, iostat=ios) dtau_guard, zero_fix
            if (ios /= 0) then
                dtau_guard = .false.; zero_fix = .false.
                ios = 0
            endif
! 新增：稳定化重正交与日志节流（可选行，缺省开启/1）
            enable_reorth = .true.
            ortho_log_stride = 10
            read(20, *, iostat=ios) enable_reorth, ortho_log_stride
            if (ios /= 0) then
                enable_reorth = .true.
                ortho_log_stride = 10
                ios = 0
            endif
! 新增：Nwrap 除数（可选行，缺省 1）
            nwrap_divider = 1
            read(20, *, iostat=ios) nwrap_divider
            if (ios /= 0) then
                nwrap_divider = 1
                ios = 0
            endif
! 新增：正交阈值与 worm 最大时间跨度（可选行，缺省 5.5e-5 / 16）
            ortho_threshold = 5.5d-5
            worm_max_span   = 4
            read(20, *, iostat=ios) ortho_threshold, worm_max_span
            if (ios /= 0) then
                ortho_threshold = 5.5d-5
                worm_max_span   = 4
                ios = 0
            endif
            close(20)
        endif
! 广播参数
        call MPI_BCAST(Beta, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(RT, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(J, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(h, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(xi, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(mu, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NB_field, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(ThetaX, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(ThetaY, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(lambda_sector, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(absolute, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(dtau_guard, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(zero_fix,   1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(enable_reorth, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(ortho_log_stride, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwrap_input, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(nwrap_divider, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(ortho_threshold, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(worm_max_span, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nlx, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nly, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Ltrot, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlxTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlyTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(LtrotTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwrap, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nbin, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwarm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nglobal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nthermal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nsweep, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_tau, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_warm, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_global, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        return
    end subroutine read_input

    subroutine Params_set()
        Lq = Nlx * Nly
        LqTherm = NlxTherm * NlyTherm
        Dtau = Beta / dble(Ltrot)
        Ndim = Lq
        if (Nwrap_input == 0) Nwrap_input = Nwrap
        if (nwrap_divider < 1) nwrap_divider = 1
        if (nwrap_divider > 1) then
            Nwrap = Nwrap_input / nwrap_divider
            if (Nwrap < 1) Nwrap = 1
        endif
        if (Nwrap <= 0) Nwrap = 1
        if (mod(Ltrot, Nwrap) == 0) then
            Nst = Ltrot / Nwrap
        else
            write(6, *) 'Ltrot必须为Nwrap的整数倍'; stop
        endif
        if (dtau_guard) then
            if (abs(RT) > Zero) then
                if (Dtau > 1.d0 / (12.d0 * abs(RT))) then
                    write(6,*) 'Dtau 守卫触发：Dtau =', Dtau, '，上限 =', 1.d0/(12.d0*abs(RT))
                    write(6,*) '请增大 Ltrot 或减小 Beta/t 以满足 Δτ ≤ 1/(12|t|)'; stop
                endif
            endif
        endif
        if (ortho_log_stride <= 0) ortho_log_stride = 1
        return
    end subroutine Params_set

    integer function nranf(iseed, N)
        integer, intent( inout ) :: iseed
        integer, intent( in ) :: N
        real(kind = 8), external :: ranf
        nranf = nint( ranf(iseed) * dble(N) + 0.5 )
        if (nranf < 1) nranf = 1
        if (nranf > N) nranf = N
        return
    end function nranf

    integer function npbc(nr, L)
        integer, intent( in ) :: nr, L
        npbc = nr
        if (nr < 1) npbc = nr + L
        if (nr > L) npbc = nr - L
        return
    end function npbc

    subroutine write_info()
        if (IRANK == 0) then
            open (unit = 50, file = 'info.txt', status = 'unknown', action = 'write')
            write(50, *) '=============================================='
            write(50, *) 'NatPhys2017 Z2 规范场 DQMC (v1)'
            write(50, *) '晶格 x 方向长度 Lx          :', Nlx
            write(50, *) '晶格 y 方向长度 Ly          :', Nly
            write(50, *) '跃迁系数 t                  :', RT
            write(50, *) '空间规范耦合 J              :', J
            write(50, *) '横场 h                      :', h
            write(50, *) '化学势 mu                   :', mu
            write(50, *) 'lambda_sector (∏λ)          :', lambda_sector
            write(50, *) 'dtau_guard / zero_fix        :', dtau_guard, zero_fix
            write(50, *) 'enable_reorth (stabilize)    :', enable_reorth
            write(50, *) 'ortho_log_stride             :', ortho_log_stride
            write(50, *) '观测与日志输出目录          : 当前工作目录（运行参数文件所在层级）'
            if (is_global) then
                write(50, *) '全局更新频率 Nglobal          :', Nglobal
            endif
            if (is_warm) then
                write(50, *) '热化步数 Nwarm               :', Nwarm
            endif
            write(50, *) '选择的初始化场构型方式 absolute :', absolute
            write(50, *) '逆温度 Beta                   :', Beta
            write(50, *) '虚时长度 Ltrot                :', Ltrot
            write(50, *) '=> Dtau = Beta / Ltrot      :', Dtau
            write(50, *) '稳定化频率 Nwrap (输入/实际) :', Nwrap_input, Nwrap
            write(50, *) 'Nwrap 除数 nwrap_divider      :', nwrap_divider
            write(50, *) 'ortho_threshold (relative)   :', ortho_threshold
            write(50, *) 'worm_max_span                :', worm_max_span
            write(50, *) '# Bins                       :', Nbin
            write(50, *) '来回扫的更新次数 Nsweep        :', Nsweep
            ! NB_field 与 Twist 已禁用
            write(50, *) '=============================================='
            call flush(50)
        endif
        return
    end subroutine write_info

end module calc_basic



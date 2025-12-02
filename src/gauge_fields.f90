module GaugeFields_mod
    use calc_basic
    use MyLattice
    implicit none
    include 'mpif.h'

    public :: GaugeConf, Gauge_conf_in, Gauge_conf_out, temporal_gamma
    private :: Gauge_conf_random_fill
    
    type :: GaugeConf
        integer, allocatable :: sigma_x(:,:) ! (Lq, Ltrot), bond i->i+ex at time τ ∈ {±1}
        integer, allocatable :: sigma_y(:,:) ! (Lq, Ltrot), bond i->i+ey at time τ ∈ {±1}
        integer, allocatable :: lambda(:)    ! (Lq), 末片轴规 λ_i ∈ {±1}
    contains
        procedure :: make => GaugeConf_make
        final     :: GaugeConf_clear
    end type GaugeConf

contains
    subroutine GaugeConf_make(this, iseed, Latt)
        class(GaugeConf), intent(inout) :: this
        integer, intent(inout) :: iseed
        class(SquareLattice), intent(in) :: Latt
        integer :: ii, nt
        real(kind=8), external :: ranf

        allocate(this%sigma_x(Lq, Ltrot))
        allocate(this%sigma_y(Lq, Ltrot))
        allocate(this%lambda(Lq))

! 初始化 σ：absolute=1 随机；3 全 +1；其余默认随机
        do nt = 1, Ltrot
            do ii = 1, Lq
                select case (absolute)
                case (3)
                    this%sigma_x(ii, nt) = 1
                    this%sigma_y(ii, nt) = 1
                case default
                    if (ranf(iseed) >= 0.5d0) then
                        this%sigma_x(ii, nt) = 1
                    else
                        this%sigma_x(ii, nt) = -1
                    endif
                    if (ranf(iseed) >= 0.5d0) then
                        this%sigma_y(ii, nt) = 1
                    else
                        this%sigma_y(ii, nt) = -1
                    endif
                end select
            enddo
        enddo

! 初始化 λ：全 +1，再按扇区修正一次（确保 ∏λ = lambda_sector）
        this%lambda = 1
        if (lambda_sector == -1) then
            this%lambda(1) = -1
        endif
        return
    end subroutine GaugeConf_make

    subroutine GaugeConf_clear(this)
        type(GaugeConf), intent(inout) :: this
        if (allocated(this%sigma_x)) deallocate(this%sigma_x)
        if (allocated(this%sigma_y)) deallocate(this%sigma_y)
        if (allocated(this%lambda))  deallocate(this%lambda)
        return
    end subroutine GaugeConf_clear

    subroutine Gauge_conf_in(Gauge, iseed, loaded)
! 兼容 AFM-ISING 的 confin/seeds 管理：
!  - confin.txt 第一行 = 0  => 冷启动，从 seeds.txt 读取每个进程的初始种子并生成随机场
!  - confin.txt 第一行 ≠ 0 => 热启动，依次读取各进程的种子与规范场配置
        type(GaugeConf), intent(inout) :: Gauge
        integer, intent(inout) :: iseed
        logical, intent(out), optional :: loaded
        integer :: ii, nt, dest
        integer :: ios_conf, ios_seed
        integer :: seed_from_file, seed_tmp
        logical :: has_loaded, cold_start
        logical :: conf_open, seeds_open
        integer, allocatable :: out_seed(:)
        integer, allocatable :: sx_bank(:,:), sy_bank(:,:), lambda_bank(:,:)
        integer, allocatable :: sigma_x_tmp(:,:), sigma_y_tmp(:,:)
        integer, allocatable :: lambda_tmp(:)
        integer :: status(MPI_STATUS_SIZE)

        has_loaded = .false.
        cold_start = .false.
        seed_from_file = 0
        ios_conf = -1
        conf_open = .false.
        seeds_open = .false.

        if (IRANK == 0) then
            ios_conf = 0
            open(unit=31, file='confin.txt', status='old', action='read', iostat=ios_conf)
            if (ios_conf == 0) then
                conf_open = .true.
                read(31, *, iostat=ios_conf) seed_from_file
                if (ios_conf == 0) then
                    if (seed_from_file == 0) then
                        cold_start = .true.
                        ios_seed = 0
                        open(unit=10, file='seeds.txt', status='old', action='read', iostat=ios_seed)
                        if (ios_seed == 0) then
                            seeds_open = .true.
                            if (ISIZE > 1) then
                                allocate(out_seed(ISIZE-1))
                                allocate(sx_bank(Lq*Ltrot, ISIZE-1))
                                allocate(sy_bank(Lq*Ltrot, ISIZE-1))
                                allocate(lambda_bank(Lq, ISIZE-1))
                            endif
                            allocate(sigma_x_tmp(Lq, Ltrot))
                            allocate(sigma_y_tmp(Lq, Ltrot))
                            allocate(lambda_tmp(Lq))

                            read(10, *, iostat=ios_seed) seed_tmp
                            if (ios_seed == 0) then
                                call Gauge_conf_random_fill(Gauge%sigma_x, Gauge%sigma_y, Gauge%lambda, seed_tmp)
                                iseed = seed_tmp
                                seed_from_file = iseed
                                if (ISIZE > 1) then
                                    do dest = 1, ISIZE - 1
                                        read(10, *, iostat=ios_seed) seed_tmp
                                        if (ios_seed /= 0) exit
                                        call Gauge_conf_random_fill(sigma_x_tmp, sigma_y_tmp, lambda_tmp, seed_tmp)
                                        out_seed(dest) = seed_tmp
                                        sx_bank(:, dest) = reshape(sigma_x_tmp, (/ Lq*Ltrot /))
                                        sy_bank(:, dest) = reshape(sigma_y_tmp, (/ Lq*Ltrot /))
                                        lambda_bank(:, dest) = lambda_tmp
                                    enddo
                                endif
                            endif
                        endif
                        if (ios_seed == 0) then
                            has_loaded = .true.
                        else
                            has_loaded = .false.
                        endif
                        if (seeds_open) then
                            close(10)
                            seeds_open = .false.
                        endif
                    else
                        iseed = seed_from_file
                        allocate(sigma_x_tmp(Lq, Ltrot))
                        allocate(sigma_y_tmp(Lq, Ltrot))
                        allocate(lambda_tmp(Lq))
                        if (ISIZE > 1) then
                            allocate(out_seed(ISIZE-1))
                            allocate(sx_bank(Lq*Ltrot, ISIZE-1))
                            allocate(sy_bank(Lq*Ltrot, ISIZE-1))
                            allocate(lambda_bank(Lq, ISIZE-1))
                        endif
                        do nt = 1, Ltrot
                            do ii = 1, Lq
                                read(31, *, iostat=ios_conf) Gauge%sigma_x(ii, nt)
                                if (ios_conf /= 0) exit
                            enddo
                            if (ios_conf /= 0) exit
                        enddo
                        if (ios_conf == 0) then
                            do nt = 1, Ltrot
                                do ii = 1, Lq
                                    read(31, *, iostat=ios_conf) Gauge%sigma_y(ii, nt)
                                    if (ios_conf /= 0) exit
                                enddo
                                if (ios_conf /= 0) exit
                            enddo
                        endif
                        if (ios_conf == 0) then
                            do ii = 1, Lq
                                read(31, *, iostat=ios_conf) Gauge%lambda(ii)
                                if (ios_conf /= 0) exit
                            enddo
                        endif
                        if (ios_conf == 0 .and. ISIZE > 1) then
                            do dest = 1, ISIZE - 1
                                read(31, *, iostat=ios_conf) seed_tmp
                                if (ios_conf /= 0) exit
                                do nt = 1, Ltrot
                                    do ii = 1, Lq
                                        read(31, *, iostat=ios_conf) sigma_x_tmp(ii, nt)
                                        if (ios_conf /= 0) exit
                                    enddo
                                    if (ios_conf /= 0) exit
                                enddo
                                if (ios_conf /= 0) exit
                                do nt = 1, Ltrot
                                    do ii = 1, Lq
                                        read(31, *, iostat=ios_conf) sigma_y_tmp(ii, nt)
                                        if (ios_conf /= 0) exit
                                    enddo
                                    if (ios_conf /= 0) exit
                                enddo
                                if (ios_conf /= 0) exit
                                do ii = 1, Lq
                                    read(31, *, iostat=ios_conf) lambda_tmp(ii)
                                    if (ios_conf /= 0) exit
                                enddo
                                if (ios_conf /= 0) exit
                                out_seed(dest) = seed_tmp
                                sx_bank(:, dest) = reshape(sigma_x_tmp, (/ Lq*Ltrot /))
                                sy_bank(:, dest) = reshape(sigma_y_tmp, (/ Lq*Ltrot /))
                                lambda_bank(:, dest) = lambda_tmp
                            enddo
                        endif
                        if (ios_conf == 0) has_loaded = .true.
                    endif
                endif
            endif
            if (conf_open) then
                close(31)
                conf_open = .false.
            endif
        endif

        call MPI_BCAST(has_loaded, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
        if (.not. has_loaded) then
            if (present(loaded)) loaded = .false.
            if (IRANK == 0) then
                if (allocated(out_seed))   deallocate(out_seed)
                if (allocated(sx_bank))    deallocate(sx_bank)
                if (allocated(sy_bank))    deallocate(sy_bank)
                if (allocated(lambda_bank)) deallocate(lambda_bank)
                if (allocated(sigma_x_tmp)) deallocate(sigma_x_tmp)
                if (allocated(sigma_y_tmp)) deallocate(sigma_y_tmp)
                if (allocated(lambda_tmp))  deallocate(lambda_tmp)
            endif
            return
        endif

        call MPI_BCAST(cold_start, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(seed_from_file, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) iseed = seed_from_file

        if (IRANK == 0) then
            if (ISIZE > 1) then
                do dest = 1, ISIZE - 1
                    call MPI_SEND(out_seed(dest), 1, MPI_Integer, dest, dest, MPI_COMM_WORLD, IERR)
                    call MPI_SEND(sx_bank(:, dest), Lq*Ltrot, MPI_Integer, dest, dest + 1024, MPI_COMM_WORLD, IERR)
                    call MPI_SEND(sy_bank(:, dest), Lq*Ltrot, MPI_Integer, dest, dest + 2048, MPI_COMM_WORLD, IERR)
                    call MPI_SEND(lambda_bank(:, dest), Lq, MPI_Integer, dest, dest + 3072, MPI_COMM_WORLD, IERR)
                enddo
            endif
        else
            call MPI_RECV(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, STATUS, IERR)
            call MPI_RECV(Gauge%sigma_x, Lq*Ltrot, MPI_Integer, 0, IRANK + 1024, MPI_COMM_WORLD, STATUS, IERR)
            call MPI_RECV(Gauge%sigma_y, Lq*Ltrot, MPI_Integer, 0, IRANK + 2048, MPI_COMM_WORLD, STATUS, IERR)
            call MPI_RECV(Gauge%lambda, Lq, MPI_Integer, 0, IRANK + 3072, MPI_COMM_WORLD, STATUS, IERR)
        endif

        if (present(loaded)) loaded = .true.

        if (IRANK == 0) then
            if (allocated(out_seed))    deallocate(out_seed)
            if (allocated(sx_bank))     deallocate(sx_bank)
            if (allocated(sy_bank))     deallocate(sy_bank)
            if (allocated(lambda_bank)) deallocate(lambda_bank)
            if (allocated(sigma_x_tmp)) deallocate(sigma_x_tmp)
            if (allocated(sigma_y_tmp)) deallocate(sigma_y_tmp)
            if (allocated(lambda_tmp))  deallocate(lambda_tmp)
        endif
        return
    end subroutine Gauge_conf_in

    subroutine Gauge_conf_out(Gauge, iseed)
        type(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: iseed
        integer :: ii, nt, src
        integer :: status(MPI_STATUS_SIZE)
        integer :: seed_tmp
        integer, allocatable :: sigma_x_tmp(:,:), sigma_y_tmp(:,:)
        integer, allocatable :: lambda_tmp(:)

        if (IRANK == 0) then
            open(unit=36, file='confout.txt', status='unknown')
            write(36,*) iseed
            do nt = 1, Ltrot
                do ii = 1, Lq
                    write(36,*) Gauge%sigma_x(ii, nt)
                enddo
            enddo
            do nt = 1, Ltrot
                do ii = 1, Lq
                    write(36,*) Gauge%sigma_y(ii, nt)
                enddo
            enddo
            do ii = 1, Lq
                write(36,*) Gauge%lambda(ii)
            enddo
            if (ISIZE > 1) then
                allocate(sigma_x_tmp(Lq, Ltrot))
                allocate(sigma_y_tmp(Lq, Ltrot))
                allocate(lambda_tmp(Lq))
                do src = 1, ISIZE - 1
                    call MPI_RECV(seed_tmp, 1, MPI_Integer, src, src, MPI_COMM_WORLD, STATUS, IERR)
                    call MPI_RECV(sigma_x_tmp, Lq*Ltrot, MPI_Integer, src, src + 1024, MPI_COMM_WORLD, STATUS, IERR)
                    call MPI_RECV(sigma_y_tmp, Lq*Ltrot, MPI_Integer, src, src + 2048, MPI_COMM_WORLD, STATUS, IERR)
                    call MPI_RECV(lambda_tmp, Lq, MPI_Integer, src, src + 3072, MPI_COMM_WORLD, STATUS, IERR)
                    write(36,*) seed_tmp
                    do nt = 1, Ltrot
                        do ii = 1, Lq
                            write(36,*) sigma_x_tmp(ii, nt)
                        enddo
                    enddo
                    do nt = 1, Ltrot
                        do ii = 1, Lq
                            write(36,*) sigma_y_tmp(ii, nt)
                        enddo
                    enddo
                    do ii = 1, Lq
                        write(36,*) lambda_tmp(ii)
                    enddo
                enddo
                deallocate(sigma_x_tmp)
                deallocate(sigma_y_tmp)
                deallocate(lambda_tmp)
            endif
            close(36)
        else
            call MPI_SEND(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, IERR)
            call MPI_SEND(Gauge%sigma_x, Lq*Ltrot, MPI_Integer, 0, IRANK + 1024, MPI_COMM_WORLD, IERR)
            call MPI_SEND(Gauge%sigma_y, Lq*Ltrot, MPI_Integer, 0, IRANK + 2048, MPI_COMM_WORLD, IERR)
            call MPI_SEND(Gauge%lambda, Lq, MPI_Integer, 0, IRANK + 3072, MPI_COMM_WORLD, IERR)
        endif
        return
    end subroutine Gauge_conf_out

    subroutine Gauge_conf_random_fill(sx, sy, lam, seed)
        integer, intent(out) :: sx(Lq, Ltrot)
        integer, intent(out) :: sy(Lq, Ltrot)
        integer, intent(out) :: lam(Lq)
        integer, intent(inout) :: seed
        integer :: ii, nt
        real(kind=8), external :: ranf

        do nt = 1, Ltrot
            do ii = 1, Lq
                select case (absolute)
                case (3)
                    sx(ii, nt) = 1
                    sy(ii, nt) = 1
                case default
                    if (ranf(seed) >= 0.5d0) then
                        sx(ii, nt) = 1
                    else
                        sx(ii, nt) = -1
                    endif
                    if (ranf(seed) >= 0.5d0) then
                        sy(ii, nt) = 1
                    else
                        sy(ii, nt) = -1
                    endif
                end select
            enddo
        enddo

        lam = 1
        if (lambda_sector == -1) lam(1) = -1
        return
    end subroutine Gauge_conf_random_fill

    real(kind=8) function temporal_gamma() result(gamma)
! γ = (1/2) ln coth(Δτ h)；h→0 时返回 0，避免数值发散
        real(kind=8) :: th
        if (abs(h) <= Zero) then
            gamma = 0.d0
        else
            th = tanh(Dtau * h)
            if (th < 1.d-300) th = 1.d-300
            gamma = -0.5d0 * log(th)
        endif
        return
    end function temporal_gamma

end module GaugeFields_mod



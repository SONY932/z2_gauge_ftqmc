module ObserMin_mod
    use calc_basic
    use MyLattice
    use GaugeFields_mod
    implicit none
    logical, parameter :: debug_green_diag_enabled = .false.

    type, public :: ObserMin
        real(kind=8) :: sum_flux
        real(kind=8) :: sum_flux2
        real(kind=8) :: sum_density
        real(kind=8) :: sum_wilson_22
        real(kind=8) :: sum_S_pi_pi
        integer :: n_acc
    contains
        procedure :: reset  => obs_reset
        procedure :: acc_gauge => obs_acc_gauge
        procedure :: acc_fermion => obs_acc_fermion
        procedure :: reduce_and_write => obs_reduce_and_write
    end type ObserMin

contains
    subroutine obs_reset(this)
        class(ObserMin), intent(inout) :: this
        this%sum_flux = 0.d0
        this%sum_flux2 = 0.d0
        this%sum_density = 0.d0
        this%sum_wilson_22 = 0.d0
        this%sum_S_pi_pi = 0.d0
        this%n_acc = 0
        return
    end subroutine obs_reset

    subroutine obs_acc_gauge(this, Gauge, Latt, nt)
        class(ObserMin), intent(inout) :: this
        class(GaugeConf), intent(in) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: nt
        integer :: x, y, i, ix, iy, flux_sum, flux_sq_sum, plaq
        real(kind=8) :: flux_avg, flux_sq_avg
        flux_sum = 0
        flux_sq_sum = 0
        do y = 1, Nly
            do x = 1, Nlx
                i  = Latt%inv_n_list(x, y)
                ix = Latt%inv_n_list(npbc(x+1, Nlx), y)
                iy = Latt%inv_n_list(x, npbc(y+1, Nly))
                plaq = Gauge%sigma_x(i, nt) * Gauge%sigma_y(ix, nt) * &
                        Gauge%sigma_x(iy, nt) * Gauge%sigma_y(i, nt)
                flux_sum = flux_sum + plaq
                flux_sq_sum = flux_sq_sum + plaq * plaq
            enddo
        enddo
        flux_avg = dble(flux_sum) / dble(Lq)
        flux_sq_avg = dble(flux_sq_sum) / dble(Lq)
        this%sum_flux = this%sum_flux + flux_avg
        this%sum_flux2 = this%sum_flux2 + flux_sq_avg
! Wilson 2x2 平均
        this%sum_wilson_22 = this%sum_wilson_22 + wilson_loop_avg(Gauge, Latt, nt, 2, 2)
        return
    end subroutine obs_acc_gauge

    subroutine obs_acc_fermion(this, GrU, GrD, Latt)
        class(ObserMin), intent(inout) :: this
        complex(kind=8), intent(in) :: GrU(:, :), GrD(:, :)
        class(SquareLattice), intent(in) :: Latt
        integer :: i
        real(kind=8) :: nloc
        real(kind=8) :: Spi
        real(kind=8) :: test_val
        
        ! 检查格林函数是否有效（无NaN或Inf）
        do i = 1, min(size(GrU,1), 10)  ! 只检查前10个对角元素以提高效率
            test_val = real(GrU(i,i)) + aimag(GrU(i,i)) + real(GrD(i,i)) + aimag(GrD(i,i))
            ! 如果检测到NaN或Inf，跳过本次采样
            if (test_val /= test_val .or. abs(test_val) > 1.d10) then
                return
            endif
        enddo
        
        nloc = 0.d0
        do i = 1, size(GrU,1)
!           等时格林函数存储的是 G = <c c†>，粒子数为 1 - G
            nloc = nloc + real( (1.d0 - GrU(i,i)) + (1.d0 - GrD(i,i)) )
        enddo
        this%sum_density = this%sum_density + nloc / dble(size(GrU,1))
! 结构因子 S(pi,pi)
        Spi = structure_factor_pi_pi(GrU, GrD, Latt)
        this%sum_S_pi_pi = this%sum_S_pi_pi + Spi
        this%n_acc = this%n_acc + 1

! 调试辅助：记录前若干次采样的对角/非对角行为，定位 G 是否被冻结
        if (debug_green_diag_enabled) call log_green_diag(GrU, GrD)
        return
    end subroutine obs_acc_fermion

    subroutine obs_reduce_and_write(this)
        class(ObserMin), intent(inout) :: this
        include 'mpif.h'
        real(kind=8) :: flux_all, flux2_all, dens_all, wil22_all, Spi_all
        real(kind=8) :: flux_mean, flux2_mean, dens_mean, wil22_mean, Spi_mean
        integer :: n_all, n_local
        flux_all = 0.d0; flux2_all = 0.d0; dens_all = 0.d0; wil22_all = 0.d0; Spi_all = 0.d0
        n_all = 0; n_local = this%n_acc
        call MPI_Reduce(this%sum_flux, flux_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        call MPI_Reduce(this%sum_flux2, flux2_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        call MPI_Reduce(this%sum_density, dens_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        call MPI_Reduce(this%sum_wilson_22, wil22_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        call MPI_Reduce(this%sum_S_pi_pi, Spi_all, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        call MPI_Reduce(n_local, n_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0 .and. n_all > 0) then
            flux_mean = flux_all / dble(n_all)
            flux2_mean = flux2_all / dble(n_all)
            dens_mean = dens_all / dble(n_all)
            wil22_mean = wil22_all / dble(n_all)
            Spi_mean = Spi_all / dble(n_all)
            open(unit=82, file='avg_flux_min', status='unknown', action='write', position='append')
            write(82,*) flux_mean
            call flush(82)
            close(82)
            open(unit=83, file='density_min', status='unknown', action='write', position='append')
            write(83,*) dens_mean
            call flush(83)
            close(83)
! 磁通易感 χ_B = Lq*(<B^2> - <B>^2)
            open(unit=86, file='chi_B', status='unknown', action='write', position='append')
            write(86,*) dble(Lq) * ( flux2_mean - flux_mean*flux_mean )
            call flush(86)
            close(86)
! Wilson 2x2 环
            open(unit=87, file='wilson_2x2', status='unknown', action='write', position='append')
            write(87,*) wil22_mean
            call flush(87)
            close(87)
! 结构因子 S(pi,pi)
            open(unit=88, file='S_pipi', status='unknown', action='write', position='append')
            write(88,*) Spi_mean
            call flush(88)
            close(88)
        endif
        call this%reset()
        return
    end subroutine obs_reduce_and_write

    function wilson_loop_avg(G, Latt, nt, Lx, Ly) result(avg)
        class(GaugeConf), intent(in) :: G
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: nt, Lx, Ly
        real(kind=8) :: avg
        integer :: x, y
        integer :: xp, yp, dx, dy
        integer :: i_src
        integer :: prod
        avg = 0.d0
        do y = 1, Nly
            do x = 1, Nlx
                prod = 1
! bottom edge along +x
                xp = x
                do dx = 0, Lx-1
                    i_src = Latt%inv_n_list(xp, y)
                    prod = prod * G%sigma_x(i_src, nt)
                    xp = npbc(xp+1, Nlx)
                enddo
! right edge along +y
                yp = y
                do dy = 0, Ly-1
                    i_src = Latt%inv_n_list(xp, yp)
                    prod = prod * G%sigma_y(i_src, nt)
                    yp = npbc(yp+1, Nly)
                enddo
! top edge along -x (Z2: inverse == itself)
                do dx = 0, Lx-1
                    i_src = Latt%inv_n_list(npbc(xp-1, Nlx), yp)
                    prod = prod * G%sigma_x(i_src, nt)
                    xp = npbc(xp-1, Nlx)
                enddo
! left edge along -y
                do dy = 0, Ly-1
                    i_src = Latt%inv_n_list(x, npbc(yp-1, Nly))
                    prod = prod * G%sigma_y(i_src, nt)
                    yp = npbc(yp-1, Nly)
                enddo
                avg = avg + dble(prod)
            enddo
        enddo
        avg = avg / dble(Lq)
        return
    end function wilson_loop_avg

    pure function structure_factor_pi_pi(GrU, GrD, Latt) result(Sq)
        complex(kind=8), intent(in) :: GrU(:, :), GrD(:, :)
        class(SquareLattice), intent(in) :: Latt
        real(kind=8) :: Sq
        integer :: i, j
        integer :: xi, yi, xj, yj, phase
        real(kind=8) :: cij
        Sq = 0.d0
        do i = 1, size(GrU,1)
            xi = Latt%n_list(i,1); yi = Latt%n_list(i,2)
            do j = 1, size(GrU,2)
                xj = Latt%n_list(j,1); yj = Latt%n_list(j,2)
                if (mod((xi - xj + yi - yj), 2) == 0) then
                    phase = 1
                else
                    phase = -1
                endif
                cij = - real( GrU(i,j)*GrU(j,i) + GrD(i,j)*GrD(j,i) )
                Sq = Sq + dble(phase) * cij
            enddo
        enddo
! 归一化到 1/Lq^2，对应标准结构因子定义
        Sq = Sq / ( dble(size(GrU,1)) * dble(size(GrU,1)) )
        return
    end function structure_factor_pi_pi

    subroutine log_green_diag(GrU, GrD)
        implicit none
        complex(kind=8), intent(in) :: GrU(:, :), GrD(:, :)
        integer, parameter :: limit = 200
        integer, save :: counter = 0
        integer :: i, j, ndim_local
        real(kind=8) :: diag_min, diag_max, diag_sum
        real(kind=8) :: diag_min_dn, diag_max_dn, diag_sum_dn
        real(kind=8) :: max_off
        if (.not. debug_green_diag_enabled) return
        ndim_local = size(GrU,1)
        if (counter >= limit) return
        if (ndim_local < 1) return

        diag_min     = real(GrU(1,1)); diag_max     = real(GrU(1,1)); diag_sum     = 0.d0
        diag_min_dn  = real(GrD(1,1)); diag_max_dn  = real(GrD(1,1)); diag_sum_dn  = 0.d0
        max_off = 0.d0

        do i = 1, ndim_local
            diag_min    = min(diag_min,    real(GrU(i,i)))
            diag_max    = max(diag_max,    real(GrU(i,i)))
            diag_sum    = diag_sum + real(GrU(i,i))
            diag_min_dn = min(diag_min_dn, real(GrD(i,i)))
            diag_max_dn = max(diag_max_dn, real(GrD(i,i)))
            diag_sum_dn = diag_sum_dn + real(GrD(i,i))
        enddo

        do i = 1, ndim_local
            do j = 1, size(GrU,2)
                if (i == j) cycle
                max_off = max(max_off, abs(GrU(i,j)), abs(GrD(i,j)))
            enddo
        enddo

        counter = counter + 1
        open(unit=199, file='debug_G_diag.log', status='unknown', action='write', position='append')
        write(199,'(I8,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6)') &
             counter, diag_min, diag_max, diag_sum/dble(ndim_local), &
             diag_min_dn, diag_max_dn, diag_sum_dn/dble(ndim_local), max_off
        close(199)
        return
    end subroutine log_green_diag

end module ObserMin_mod


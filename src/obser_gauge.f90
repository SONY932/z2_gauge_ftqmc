module ObserGauge_mod
    use calc_basic
    use MyLattice
    use GaugeFields_mod
    use GaugeAction_mod
    implicit none

    type, public :: ObserGauge
        real(kind=8) :: avg_flux
    contains
        procedure :: reset => ObserGauge_reset
        procedure :: accumulate => ObserGauge_acc
        procedure :: writeout => ObserGauge_write
    end type ObserGauge

contains
    subroutine ObserGauge_reset(this)
        class(ObserGauge), intent(inout) :: this
        this%avg_flux = 0.d0
        return
    end subroutine ObserGauge_reset

    subroutine ObserGauge_acc(this, Gauge, Latt, nt)
        class(ObserGauge), intent(inout) :: this
        class(GaugeConf), intent(in) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: nt
        integer :: x, y
        integer :: flux_sum
        flux_sum = 0
        do y = 1, Nly
            do x = 1, Nlx
                flux_sum = flux_sum + plaq_product(Gauge, Latt, x, y, nt)
            enddo
        enddo
        this%avg_flux = this%avg_flux + dble(flux_sum) / dble(Lq)
        return
    end subroutine ObserGauge_acc

    subroutine ObserGauge_write(this)
        class(ObserGauge), intent(in) :: this
        open(unit=81, file='avg_flux', status='unknown', action='write', position='append')
        write(81,*) this%avg_flux
        close(81)
        return
    end subroutine ObserGauge_write

end module ObserGauge_mod



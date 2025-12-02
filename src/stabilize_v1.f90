module StabilizeV1_mod
    use calc_basic
    implicit none

    type, public :: WrapCounter
        integer :: step
    contains
        procedure :: reset => Wrap_reset
        procedure :: tick  => Wrap_tick
        procedure :: should_wrap => Wrap_should
    end type WrapCounter

contains
    subroutine Wrap_reset(this)
        class(WrapCounter), intent(inout) :: this
        this%step = 0
        return
    end subroutine Wrap_reset

    subroutine Wrap_tick(this)
        class(WrapCounter), intent(inout) :: this
        this%step = this%step + 1
        return
    end subroutine Wrap_tick

    logical function Wrap_should(this) result(flag)
        class(WrapCounter), intent(in) :: this
        if (Nwrap <= 0) then
            flag = .false.
        else
            flag = (mod(this%step, Nwrap) == 0)
        endif
        return
    end function Wrap_should

end module StabilizeV1_mod



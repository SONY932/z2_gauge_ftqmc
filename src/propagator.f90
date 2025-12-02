module Propagator_mod
    use calc_basic
    implicit none

    type, public :: Propagator
        complex(kind=8), allocatable :: Gr(:, :)
    contains
        procedure :: make => Prop_make
        procedure :: reset_half => Prop_reset_half
    end type Propagator

contains
    subroutine Prop_make(this)
        class(Propagator), intent(inout) :: this
        allocate(this%Gr(Ndim, Ndim))
        call Prop_reset_half(this)
        return
    end subroutine Prop_make

    subroutine Prop_reset_half(this)
        class(Propagator), intent(inout) :: this
        integer :: i
        this%Gr = dcmplx(0.d0, 0.d0)
        do i = 1, Ndim
            this%Gr(i, i) = dcmplx(0.5d0, 0.d0)
        enddo
        return
    end subroutine Prop_reset_half

end module Propagator_mod



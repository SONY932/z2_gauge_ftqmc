module OperatorGauge_mod
    use calc_basic
    implicit none

    type, public :: OperatorGauge
        real(kind=8) :: cval   ! cosh(scale*Δτ t)
        real(kind=8) :: sval   ! sinh(scale*Δτ t)
    contains
        procedure :: set => opG_set
        procedure :: params => opG_params
        procedure :: mmult_R_2x2 => opG_mmult_R_2x2
        procedure :: mmult_L_2x2 => opG_mmult_L_2x2
    end type OperatorGauge

contains
    subroutine opG_set(this, scale)
        class(OperatorGauge), intent(inout) :: this
        real(kind=8), intent(in), optional :: scale
        real(kind=8) :: sc
        sc = 1.d0
        if (present(scale)) sc = scale
        this%cval = cosh(sc * Dtau * RT)
        this%sval = sinh(sc * Dtau * RT)
        return
    end subroutine opG_set

    subroutine opG_params(this, sigma, c, s)
        class(OperatorGauge), intent(in) :: this
        integer, intent(in) :: sigma  ! ±1
        real(kind=8), intent(out) :: c, s
        c = this%cval
        s = dble(sigma) * this%sval
        return
    end subroutine opG_params

    subroutine opG_mmult_R_2x2(this, Mat, i, j, sigma, nflag)
! In:  Mat; Out: O(σ) * Mat for nflag=+1; O(σ)^{-1} * Mat for nflag=-1
        class(OperatorGauge), intent(in) :: this
        complex(kind=8), intent(inout) :: Mat(:, :)
        integer, intent(in) :: i, j, sigma, nflag
        complex(kind=8) :: row_i(size(Mat,2)), row_j(size(Mat,2))
        real(kind=8) :: c, s
        call this%params(sigma, c, s)
        if (nflag == -1) s = -s
        row_i = Mat(i, :)
        row_j = Mat(j, :)
        Mat(i, :) = dcmplx(c,0.d0) * row_i + dcmplx(s,0.d0) * row_j
        Mat(j, :) = dcmplx(s,0.d0) * row_i + dcmplx(c,0.d0) * row_j
        return
    end subroutine opG_mmult_R_2x2

    subroutine opG_mmult_L_2x2(this, Mat, i, j, sigma, nflag)
! In:  Mat; Out: Mat * O(σ) for nflag=+1; Mat * O(σ)^{-1} for nflag=-1
        class(OperatorGauge), intent(in) :: this
        complex(kind=8), intent(inout) :: Mat(:, :)
        integer, intent(in) :: i, j, sigma, nflag
        complex(kind=8) :: col_i(size(Mat,1)), col_j(size(Mat,1))
        real(kind=8) :: c, s
        call this%params(sigma, c, s)
        if (nflag == -1) s = -s
        col_i = Mat(:, i)
        col_j = Mat(:, j)
        Mat(:, i) = dcmplx(c,0.d0) * col_i + dcmplx(s,0.d0) * col_j
        Mat(:, j) = dcmplx(s,0.d0) * col_i + dcmplx(c,0.d0) * col_j
        return
    end subroutine opG_mmult_L_2x2

end module OperatorGauge_mod


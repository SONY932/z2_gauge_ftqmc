module OpMu_mod
    use calc_basic
    implicit none

contains
    subroutine opMu_mmult_R(Mat, nflag)
! In: Mat; Out: exp((+nflag)*Δτ μ I) * Mat
! 与 apply_group_* 的左右乘保持一致，正向/逆向调用提供成对的 ±μ 半步
        complex(kind=8), intent(inout) :: Mat(:, :)
        integer, intent(in) :: nflag
        complex(kind=8) :: s
        integer :: i
        if (nflag == 1) then
            s = dcmplx(exp( Dtau * mu), 0.d0)
        elseif (nflag == -1) then
            s = dcmplx(exp(-Dtau * mu), 0.d0)
        else
            write(6,*) 'opMu_mmult_R: illegal nflag=', nflag; stop
        endif
        do i = 1, size(Mat,1)
            Mat(i, :) = s * Mat(i, :)
        enddo
        return
    end subroutine opMu_mmult_R

    subroutine opMu_mmult_L(Mat, nflag)
! In: Mat; Out: Mat * exp((+nflag)*Δτ μ I)
! 采用与右乘相同的指数因子，保证左右传播器互为厄米共轭
        complex(kind=8), intent(inout) :: Mat(:, :)
        integer, intent(in) :: nflag
        complex(kind=8) :: s
        integer :: j
        if (nflag == 1) then
            s = dcmplx(exp( Dtau * mu), 0.d0)
        elseif (nflag == -1) then
            s = dcmplx(exp(-Dtau * mu), 0.d0)
        else
            write(6,*) 'opMu_mmult_L: illegal nflag=', nflag; stop
        endif
        do j = 1, size(Mat,2)
            Mat(:, j) = Mat(:, j) * s
        enddo
        return
    end subroutine opMu_mmult_L

end module OpMu_mod


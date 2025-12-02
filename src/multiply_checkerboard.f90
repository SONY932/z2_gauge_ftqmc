module MultiplyChecker_mod
    use calc_basic
    use MyLattice
    use BondList_mod
    use GaugeFields_mod
    use OperatorGauge_mod
    implicit none
    logical, parameter :: debug_operator_enabled = .false.
    integer, save :: debug_operator_write_count = 0
    integer, save :: debug_operator_call_count = 0
    integer, save :: debug_operator_detail_count = 0
    integer, parameter :: debug_operator_detail_target_nt = 6
    integer, parameter :: debug_operator_detail_limit = 64
contains
    subroutine apply_group_R(Mat, Latt, Bonds, Gauge, nt, group_no, dir, nflag, scale)
! Mat <- Mat * Π_{b∈group(group_no,dir)} O_b(σ_b,τ)^{±1}
        complex(kind=8), intent(inout) :: Mat(:, :)
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt, group_no, nflag
        character(len=*), intent(in) :: dir ! 'x' or 'y'
        real(kind=8), intent(in), optional :: scale
        type(OperatorGauge) :: Op
        integer :: k, bidx, i, j, sigma
        if (present(scale)) then
            call Op%set(scale)
        else
            call Op%set()
        endif
        select case (dir)
        case ('x')
            if (group_no == 1) then
                do k = 1, size(Bonds%group_ex_even)
                    bidx = Bonds%group_ex_even(k)
                    i = Bonds%ex_src(bidx); j = Bonds%ex_dst(bidx)
                    sigma = Gauge%sigma_x(i, nt)
                    call Op%mmult_R_2x2(Mat, i, j, sigma, nflag)
                enddo
            else
                do k = 1, size(Bonds%group_ex_odd)
                    bidx = Bonds%group_ex_odd(k)
                    i = Bonds%ex_src(bidx); j = Bonds%ex_dst(bidx)
                    sigma = Gauge%sigma_x(i, nt)
                    call Op%mmult_R_2x2(Mat, i, j, sigma, nflag)
                enddo
            endif
        case ('y')
            if (group_no == 1) then
                do k = 1, size(Bonds%group_ey_even)
                    bidx = Bonds%group_ey_even(k)
                    i = Bonds%ey_src(bidx); j = Bonds%ey_dst(bidx)
                    sigma = Gauge%sigma_y(i, nt)
                    call Op%mmult_R_2x2(Mat, i, j, sigma, nflag)
                enddo
            else
                do k = 1, size(Bonds%group_ey_odd)
                    bidx = Bonds%group_ey_odd(k)
                    i = Bonds%ey_src(bidx); j = Bonds%ey_dst(bidx)
                    sigma = Gauge%sigma_y(i, nt)
                    call Op%mmult_R_2x2(Mat, i, j, sigma, nflag)
                enddo
            endif
        case default
            write(6,*) 'apply_group_R: illegal dir=', dir; stop
        end select
        return
    end subroutine apply_group_R

    subroutine apply_group_L(Mat, Latt, Bonds, Gauge, nt, group_no, dir, nflag, scale)
! Mat <- Π_{b∈group(group_no,dir)} O_b(σ_b,τ)^{±1} * Mat
        complex(kind=8), intent(inout) :: Mat(:, :)
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt, group_no, nflag
        character(len=*), intent(in) :: dir ! 'x' or 'y'
        real(kind=8), intent(in), optional :: scale
        type(OperatorGauge) :: Op
        integer :: k, bidx, i, j, sigma
        real(kind=8) :: c, s
        if (present(scale)) then
            call Op%set(scale)
        else
            call Op%set()
        endif
        select case (dir)
        case ('x')
            if (group_no == 1) then
                do k = 1, size(Bonds%group_ex_even)
                    bidx = Bonds%group_ex_even(k)
                    i = Bonds%ex_src(bidx); j = Bonds%ex_dst(bidx)
                    sigma = Gauge%sigma_x(i, nt)
                    call check_sigma_input('apply_group_L', dir, nt, group_no, k, i, j, sigma, size(Mat,1))
                    if (nflag == -1) then
                        call Op%params(sigma, c, s)
                        s = -s
                        call debug_log_operator_detail('pre', 'L_x_even', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                    endif
                    call Op%mmult_L_2x2(Mat, i, j, sigma, nflag)
                    if (nflag == -1) then
                        call debug_log_operator_detail('post', 'L_x_even', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                        call debug_log_operator('L_x_even', nt, group_no, k, i, j, sigma, Mat)
                    endif
                enddo
            else
                do k = 1, size(Bonds%group_ex_odd)
                    bidx = Bonds%group_ex_odd(k)
                    i = Bonds%ex_src(bidx); j = Bonds%ex_dst(bidx)
                    sigma = Gauge%sigma_x(i, nt)
                    call check_sigma_input('apply_group_L', dir, nt, group_no, k, i, j, sigma, size(Mat,1))
                    if (nflag == -1) then
                        call Op%params(sigma, c, s)
                        s = -s
                        call debug_log_operator_detail('pre', 'L_x_odd', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                    endif
                    call Op%mmult_L_2x2(Mat, i, j, sigma, nflag)
                    if (nflag == -1) then
                        call debug_log_operator_detail('post', 'L_x_odd', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                        call debug_log_operator('L_x_odd', nt, group_no, k, i, j, sigma, Mat)
                    endif
                enddo
            endif
        case ('y')
            if (group_no == 1) then
                do k = 1, size(Bonds%group_ey_even)
                    bidx = Bonds%group_ey_even(k)
                    i = Bonds%ey_src(bidx); j = Bonds%ey_dst(bidx)
                    sigma = Gauge%sigma_y(i, nt)
                    call check_sigma_input('apply_group_L', dir, nt, group_no, k, i, j, sigma, size(Mat,1))
                    if (nflag == -1) then
                        call Op%params(sigma, c, s)
                        s = -s
                        call debug_log_operator_detail('pre', 'L_y_even', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                    endif
                    call Op%mmult_L_2x2(Mat, i, j, sigma, nflag)
                    if (nflag == -1) then
                        call debug_log_operator_detail('post', 'L_y_even', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                        call debug_log_operator('L_y_even', nt, group_no, k, i, j, sigma, Mat)
                    endif
                enddo
            else
                do k = 1, size(Bonds%group_ey_odd)
                    bidx = Bonds%group_ey_odd(k)
                    i = Bonds%ey_src(bidx); j = Bonds%ey_dst(bidx)
                    sigma = Gauge%sigma_y(i, nt)
                    call check_sigma_input('apply_group_L', dir, nt, group_no, k, i, j, sigma, size(Mat,1))
                    if (nflag == -1) then
                        call Op%params(sigma, c, s)
                        s = -s
                        call debug_log_operator_detail('pre', 'L_y_odd', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                    endif
                    call Op%mmult_L_2x2(Mat, i, j, sigma, nflag)
                    if (nflag == -1) then
                        call debug_log_operator_detail('post', 'L_y_odd', nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
                        call debug_log_operator('L_y_odd', nt, group_no, k, i, j, sigma, Mat)
                    endif
                enddo
            endif
        case default
            write(6,*) 'apply_group_L: illegal dir=', dir; stop
        end select
        return
    end subroutine apply_group_L

    subroutine check_sigma_input(tag, dir, nt, group_no, k, i, j, sigma, nmax)
        character(len=*), intent(in) :: tag
        character(len=*), intent(in) :: dir
        integer, intent(in) :: nt, group_no, k, i, j, sigma, nmax
        if (i < 1 .or. i > nmax .or. j < 1 .or. j > nmax) then
            write(6,*) tag, ': index out of range ->', dir, nt, group_no, k, i, j, ' nmax=', nmax
            stop 'apply_group index'
        endif
        if (sigma /= 1 .and. sigma /= -1) then
            write(6,*) tag, ': sigma invalid ->', dir, nt, group_no, k, i, j, sigma
            stop 'apply_group sigma'
        endif
        return
    end subroutine check_sigma_input

    subroutine reset_debug_operator()
        implicit none
        debug_operator_write_count = 0
        debug_operator_call_count = 0
        debug_operator_detail_count = 0
        if (.not. debug_operator_enabled) return
        open(unit=211, file='debug_operator.log', status='unknown', action='write')
        close(211)
        open(unit=212, file='debug_operator_detail.log', status='unknown', action='write')
        close(212)
        return
    end subroutine reset_debug_operator

    subroutine debug_log_operator(tag, nt, group_no, k, i, j, sigma, Mat)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt, group_no, k, i, j, sigma
        complex(kind=8), intent(in) :: Mat(:, :)
        integer :: n, r, c
        real(kind=8) :: max_abs
        logical :: has_nan
        if (.not. debug_operator_enabled) return
        debug_operator_call_count = debug_operator_call_count + 1
        n = size(Mat,1)
        if (n > 0) then
            has_nan = .false.
            max_abs = 0.d0
            do r = 1, n
                do c = 1, size(Mat,2)
                    if (.not.(real(Mat(r,c)) == real(Mat(r,c)))) has_nan = .true.
                    if (.not.(aimag(Mat(r,c)) == aimag(Mat(r,c)))) has_nan = .true.
                    max_abs = max(max_abs, abs(Mat(r,c)))
                enddo
            enddo
        else
            has_nan = .true.
            max_abs = 0.d0
        endif
        if (has_nan) then
            open(unit=211, file='debug_operator.log', status='unknown', action='write', position='append')
            write(211,'(I8,1X,A20,1X,I6,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,A3,1X,ES14.6)') &
                debug_operator_call_count, tag, nt, group_no, k, i, j, sigma, 'NaN', max_abs
            close(211)
            stop 'debug_operator detected NaN'
        endif
        debug_operator_write_count = debug_operator_write_count + 1
        open(unit=211, file='debug_operator.log', status='unknown', action='write', position='append')
        write(211,'(I8,1X,A20,1X,I6,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,ES14.6)') &
            debug_operator_call_count, tag, nt, group_no, k, i, j, sigma, max_abs
        close(211)
        return
    end subroutine debug_log_operator

    subroutine debug_log_operator_detail(stage, tag, nt, group_no, k, i, j, sigma, nflag, c, s, Mat)
        character(len=*), intent(in) :: stage, tag
        integer, intent(in) :: nt, group_no, k, i, j, sigma, nflag
        real(kind=8), intent(in) :: c, s
        complex(kind=8), intent(in) :: Mat(:, :)
        complex(kind=8) :: diag_i, diag_j, off_ij, off_ji
        real(kind=8) :: row_i_max, row_j_max, mat_max
        if (.not. debug_operator_enabled) return
        if (nt /= debug_operator_detail_target_nt) return
        if (debug_operator_detail_count >= debug_operator_detail_limit) return
        debug_operator_detail_count = debug_operator_detail_count + 1
        diag_i = Mat(i, i)
        diag_j = Mat(j, j)
        off_ij = Mat(i, j)
        off_ji = Mat(j, i)
        row_i_max = maxval(abs(Mat(i, :)))
        row_j_max = maxval(abs(Mat(j, :)))
        mat_max = maxval(abs(Mat))
        open(unit=212, file='debug_operator_detail.log', status='unknown', action='write', position='append')
        write(212,*) '#', debug_operator_detail_count, trim(stage), trim(tag), 'nt', nt, 'group', group_no, 'k', k, 'i', i, 'j', j, &
            'sigma', sigma, 'nflag', nflag
        write(212,*) 'c', c, 's', s, 'row_i_max', row_i_max, 'row_j_max', row_j_max, 'mat_max', mat_max
        write(212,*) 'Mat(i,i)', diag_i, 'Mat(j,j)', diag_j, 'Mat(i,j)', off_ij, 'Mat(j,i)', off_ji
        close(212)
        return
    end subroutine debug_log_operator_detail

end module MultiplyChecker_mod



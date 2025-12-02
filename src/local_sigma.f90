module LocalSigma_mod
    use calc_basic
    use MyLattice
    use BondList_mod
    use GaugeFields_mod
    use MultiplyChecker_mod
    implicit none
    logical, parameter :: debug_sigma_log_enabled = .false.

contains
    subroutine apply_lambda_left(Gr, Gauge)
! 在末片左乘 P[λ]
        complex(kind=8), intent(inout) :: Gr(:, :)
        class(GaugeConf), intent(in) :: Gauge
        integer :: i
        do i = 1, size(Gr,1)
            Gr(i, :) = cmplx(dble(Gauge%lambda(i)), 0.d0, kind=8) * Gr(i, :)
        enddo
        return
    end subroutine apply_lambda_left

    subroutine apply_lambda_right(Gr, Gauge)
! 在末片右乘 P[λ]（对称化λ投影）
        complex(kind=8), intent(inout) :: Gr(:, :)
        class(GaugeConf), intent(in) :: Gauge
        integer :: i
        do i = 1, size(Gr,2)
            Gr(:, i) = Gr(:, i) * cmplx(dble(Gauge%lambda(i)), 0.d0, kind=8)
        enddo
        return
    end subroutine apply_lambda_right

    subroutine apply_lambda_both(Gr, Gauge)
! 先左乘再右乘 λ，对等时格林函数做一次扇区投影（再次调用即可恢复）
        complex(kind=8), intent(inout) :: Gr(:, :)
        class(GaugeConf), intent(in) :: Gauge
        call apply_lambda_left(Gr, Gauge)
        call apply_lambda_right(Gr, Gauge)
        return
    end subroutine apply_lambda_both

    subroutine apply_trotter_layer_R(Gr, Latt, Bonds, Gauge, nt, nflag_in)
! 对等时 Green 矩阵执行一次"右乘"的对称二阶 Trotter 层
! 关键修复：顺序应该与_L相反（对称性要求）
        complex(kind=8), intent(inout) :: Gr(:, :)
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        integer, intent(in), optional :: nflag_in
        integer :: nflag
        real(kind=8), parameter :: half = 0.5d0
        nflag = 1
        if (present(nflag_in)) nflag = nflag_in
        ! 对称Trotter的序列是回文，正反序相同
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
        call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        return
    end subroutine apply_trotter_layer_R

    subroutine apply_half_trotter_L(Gr, Latt, Bonds, Gauge, nt, nflag_in, reverse)
! 左乘半步Trotter：exp(±ε/2 K_gauge)，nflag_in控制符号
        complex(kind=8), intent(inout) :: Gr(:, :)
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        integer, intent(in), optional :: nflag_in
        logical, intent(in), optional :: reverse
        integer :: nflag
        real(kind=8), parameter :: half = 0.5d0
        logical :: rev
        nflag = 1
        if (present(nflag_in)) nflag = nflag_in
        rev = .false.
        if (present(reverse)) rev = reverse
        if (.not. rev) then
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
        else
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
        endif
        return
    end subroutine apply_half_trotter_L

    subroutine apply_half_trotter_R(Gr, Latt, Bonds, Gauge, nt, nflag_in, reverse)
! 右乘半步Trotter：exp(±ε/2 K_gauge)，nflag_in控制符号
        complex(kind=8), intent(inout) :: Gr(:, :)
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        integer, intent(in), optional :: nflag_in
        logical, intent(in), optional :: reverse
        integer :: nflag
        real(kind=8), parameter :: half = 0.5d0
        logical :: rev
        nflag = 1
        if (present(nflag_in)) nflag = nflag_in
        rev = .false.
        if (present(reverse)) rev = reverse
        if (.not. rev) then
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
        else
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
        endif
        return
    end subroutine apply_half_trotter_R

    subroutine apply_trotter_layer_L(Gr, Latt, Bonds, Gauge, nt, nflag_in)
! 对等时 Green 矩阵执行一次"右乘"的对称二阶 Trotter 层
! 注意：由于函数名混乱，apply_group_L 实际执行右乘！
! 为了产生与 apply_trotter_layer_R 相同的 B_gauge，需要使用**反向**顺序
! apply_trotter_layer_R 使用：2x, 1x, 2y, 1y, 1y, 2y, 1x, 2x（左乘顺序）
! 对于右乘，后应用的算符在左边，所以应该反向应用：2x, 1x, 2y, 1y, 1y, 2y, 1x, 2x
! 这样产生 Mat * (O_x2 * O_x1 * ... * O_x2) = Mat * B_gauge
        complex(kind=8), intent(inout) :: Gr(:, :)
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: nt
        integer, intent(in), optional :: nflag_in
        integer :: nflag
        real(kind=8), parameter :: half = 0.5d0
        nflag = 1
        if (present(nflag_in)) nflag = nflag_in
! 反向顺序以匹配 apply_trotter_layer_R 产生的 B_gauge
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
        call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        if (debug_sigma_log_enabled) call debug_log_sigma('after_apply_trotter_L', nt, Gr)
        return
    end subroutine apply_trotter_layer_L

    subroutine debug_log_sigma(tag, nt, Mat)
        implicit none
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        complex(kind=8), intent(in) :: Mat(:, :)
        integer, parameter :: limit = 400
        integer, save :: counter = 0
        integer :: i, j, n
        real(kind=8) :: diag_min, diag_max, diag_avg, max_off
        logical :: has_nan

        if (.not. debug_sigma_log_enabled) return
        n = size(Mat,1)
        if (counter >= limit .or. n <= 0) return

        has_nan = .false.
        diag_min = 0.d0; diag_max = 0.d0; diag_avg = 0.d0
        max_off = 0.d0

        do i = 1, n
            if (.not.(real(Mat(i,i)) == real(Mat(i,i)))) has_nan = .true.
            if (.not.(aimag(Mat(i,i)) == aimag(Mat(i,i)))) has_nan = .true.
        enddo

        if (.not. has_nan) then
            diag_min = real(Mat(1,1))
            diag_max = real(Mat(1,1))
            diag_avg = 0.d0
            do i = 1, n
                diag_min = min(diag_min, real(Mat(i,i)))
                diag_max = max(diag_max, real(Mat(i,i)))
                diag_avg = diag_avg + real(Mat(i,i))
            enddo
            diag_avg = diag_avg / dble(n)
        endif

        do i = 1, n
            do j = 1, size(Mat,2)
                if (i == j) cycle
                if (.not.(real(Mat(i,j)) == real(Mat(i,j)))) has_nan = .true.
                if (.not.(aimag(Mat(i,j)) == aimag(Mat(i,j)))) has_nan = .true.
                max_off = max(max_off, abs(Mat(i,j)))
            enddo
        enddo

        counter = counter + 1
        open(unit=210, file='debug_sigma.log', status='unknown', action='write', position='append')
        if (has_nan) then
            write(210,'(I8,1X,A24,1X,I8,1X,A3,1X,ES14.6)') counter, tag, nt, 'NaN', max_off
        else
            write(210,'(I8,1X,A24,1X,I8,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6)') &
                counter, tag, nt, diag_min, diag_max, diag_avg, max_off
        endif
        close(210)
        return
    end subroutine debug_log_sigma

end module LocalSigma_mod

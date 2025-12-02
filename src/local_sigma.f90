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
! 当 nflag=-1 时（逆向），顺序需要反转以得到正确的 B^{-1}
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
        if (nflag == 1) then
! 正向 G * B：标准顺序
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        else
! 逆向 G * B^{-1}：反转顺序（B^{-1} = O_n^{-1} * ... * O_1^{-1}）
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        endif
        return
    end subroutine apply_trotter_layer_R

    subroutine apply_half_trotter_L(Gr, Latt, Bonds, Gauge, nt, nflag_in, reverse)
! 左乘半步Trotter：exp(±ε/2 K_gauge)，nflag_in控制符号
! 顺序与 apply_trotter_layer_L 一致：
!   前半(rev=false): x-even, x-odd, y-even, y-odd
!   后半(rev=true):  y-odd, y-even, x-odd, x-even
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
! 前半：x-even(1), x-odd(2), y-even(1), y-odd(2)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
        else
! 后半：y-odd(2), y-even(1), x-odd(2), x-even(1)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
        endif
        return
    end subroutine apply_half_trotter_L

    subroutine apply_half_trotter_R(Gr, Latt, Bonds, Gauge, nt, nflag_in, reverse)
! 右乘半步Trotter：exp(±ε/2 K_gauge)，nflag_in控制符号
! 顺序与 apply_trotter_layer_R 一致：
!   前半(rev=false): x-odd, x-even, y-odd, y-even
!   后半(rev=true):  y-even, y-odd, x-even, x-odd
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
! 前半：x-odd(2), x-even(1), y-odd(2), y-even(1)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
        else
! 后半：y-even(1), y-odd(2), x-even(1), x-odd(2)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_R(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        endif
        return
    end subroutine apply_half_trotter_R

    subroutine apply_trotter_layer_L(Gr, Latt, Bonds, Gauge, nt, nflag_in)
! 对等时 Green 矩阵执行一次"左乘"的对称二阶 Trotter 层
! 当 nflag=-1 时（逆向），顺序需要反转以得到正确的 B^{-1}
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
        if (nflag == 1) then
! 正向 B * G：使用与 apply_trotter_layer_R 反转的顺序（因为左乘从右到左累积）
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        else
! 逆向 B^{-1} * G：使用与正向反转的顺序（先应用最后一个因子的逆）
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'y', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 1, 'x', nflag, half)
            call apply_group_L(Gr, Latt, Bonds, Gauge, nt, 2, 'x', nflag, half)
        endif
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

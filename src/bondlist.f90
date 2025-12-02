module BondList_mod
    use calc_basic
    use MyLattice
    implicit none

    type, public :: BondList
        integer :: num_bonds_ex, num_bonds_ey
        integer, allocatable :: ex_src(:), ex_dst(:)
        integer, allocatable :: ey_src(:), ey_dst(:)
        integer, allocatable :: group_ex_even(:), group_ex_odd(:)
        integer, allocatable :: group_ey_even(:), group_ey_odd(:)
    contains
        procedure :: make  => BondList_make
        final     :: BondList_clear
    end type BondList

contains
    subroutine BondList_make(this, Latt)
        class(BondList), intent(inout) :: this
        class(SquareLattice), intent(in) :: Latt
        integer :: ii, ix, iy, bcount
        integer :: even_count, odd_count

! x-方向：每个格点 i 的 (i -> i+ex)
        this%num_bonds_ex = Lq
        allocate(this%ex_src(this%num_bonds_ex), this%ex_dst(this%num_bonds_ex))
        bcount = 0
        do ii = 1, Lq
            bcount = bcount + 1
            this%ex_src(bcount) = ii
            this%ex_dst(bcount) = Latt%L_bonds(ii, 1)
        enddo
! y-方向：每个格点 i 的 (i -> i+ey)
        this%num_bonds_ey = Lq
        allocate(this%ey_src(this%num_bonds_ey), this%ey_dst(this%num_bonds_ey))
        bcount = 0
        do ii = 1, Lq
            bcount = bcount + 1
            this%ey_src(bcount) = ii
            this%ey_dst(bcount) = Latt%L_bonds(ii, 2)
        enddo

! 四组棋盘分解：
! ex-偶列/奇列：按源点 x 坐标的奇偶划分
        even_count = 0; odd_count = 0
        do ii = 1, this%num_bonds_ex
            ix = Latt%n_list(this%ex_src(ii), 1)
            if (mod(ix, 2) == 0) then
                even_count = even_count + 1
            else
                odd_count = odd_count + 1
            endif
        enddo
        allocate(this%group_ex_even(even_count), this%group_ex_odd(odd_count))
        even_count = 0; odd_count = 0
        do ii = 1, this%num_bonds_ex
            ix = Latt%n_list(this%ex_src(ii), 1)
            if (mod(ix, 2) == 0) then
                even_count = even_count + 1
                this%group_ex_even(even_count) = ii
            else
                odd_count = odd_count + 1
                this%group_ex_odd(odd_count) = ii
            endif
        enddo

! ey-偶行/奇行：按源点 y 坐标的奇偶划分
        even_count = 0; odd_count = 0
        do ii = 1, this%num_bonds_ey
            iy = Latt%n_list(this%ey_src(ii), 2)
            if (mod(iy, 2) == 0) then
                even_count = even_count + 1
            else
                odd_count = odd_count + 1
            endif
        enddo
        allocate(this%group_ey_even(even_count), this%group_ey_odd(odd_count))
        even_count = 0; odd_count = 0
        do ii = 1, this%num_bonds_ey
            iy = Latt%n_list(this%ey_src(ii), 2)
            if (mod(iy, 2) == 0) then
                even_count = even_count + 1
                this%group_ey_even(even_count) = ii
            else
                odd_count = odd_count + 1
                this%group_ey_odd(odd_count) = ii
            endif
        enddo
        return
    end subroutine BondList_make

    subroutine BondList_clear(this)
        type(BondList), intent(inout) :: this
        if (allocated(this%ex_src)) deallocate(this%ex_src)
        if (allocated(this%ex_dst)) deallocate(this%ex_dst)
        if (allocated(this%ey_src)) deallocate(this%ey_src)
        if (allocated(this%ey_dst)) deallocate(this%ey_dst)
        if (allocated(this%group_ex_even)) deallocate(this%group_ex_even)
        if (allocated(this%group_ex_odd))  deallocate(this%group_ex_odd)
        if (allocated(this%group_ey_even)) deallocate(this%group_ey_even)
        if (allocated(this%group_ey_odd))  deallocate(this%group_ey_odd)
        return
    end subroutine BondList_clear

end module BondList_mod



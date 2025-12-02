module MyLattice ! definition on space geometry
    use calc_basic
    implicit none
    
    type, public :: SquareLattice
        integer, dimension(:,:), allocatable :: n_list, inv_n_list, t_list, inv_t_list
        integer, dimension(:,:), allocatable :: L_bonds, nn_bonds, LT_bonds, imj
        real(kind=8), dimension(:,:), allocatable :: xk_v, aimj_v, k_dot_r
        real(kind=8) :: a1_v(2), a2_v(2), b1_v(2), b2_v(2)
    contains
        final :: Lattice_clear
    end type SquareLattice
    
    complex(kind=8), dimension(:,:), public, allocatable, save :: ZKRON
    
    interface Fourier_R_to_K
        module procedure FFT_RtoK_1, FFT_RtoK_2, FFT_RtoK_3,  FFT_RtoK_4
    end interface
    
contains
    subroutine Lattice_make(Latt)
        class(SquareLattice), intent(inout) :: Latt
        integer :: nc, nx, ny
        integer :: ii, jj, ix, jx, iy, jy, nt, iit, imjx, imjy
! 定义晶格坐标到格点编号的映射
        allocate(Latt%n_list(Lq, 1:2), Latt%inv_n_list(Nlx, Nly))
        nc = 0
        do ny = 1, Nly
            do nx = 1, Nlx
                nc = nc + 1
                Latt%n_list(nc, 1) = nx
                Latt%n_list(nc, 2) = ny
                Latt%inv_n_list(nx, ny) = nc
            enddo
        enddo
! 定义虚时间和格点编号到总 ising 场采样空间的编号的映射
        allocate(Latt%t_list(Lq*Ltrot, 2), Latt%inv_t_list(Lq, Ltrot))
        nc = 0
        do nt = 1, Ltrot
            do ii = 1, Lq
                nc = nc + 1
                Latt%t_list(nc, 1) = ii
                Latt%t_list(nc, 2) = nt
                Latt%inv_t_list(ii, nt) = nc
            enddo
        enddo
! 定义两个格点差到格点编号的映射（利用平移对称性）
        allocate(Latt%imj(Lq, Lq))
        do jj = 1, Lq
            do ii =1, Lq
                ix = Latt%n_list(ii, 1)
                iy = Latt%n_list(ii, 2)
                jx = Latt%n_list(jj, 1)
                jy = Latt%n_list(jj, 2)
                imjx = npbc(ix - jx, Nlx)
                imjy = npbc(iy - jy, Nly)
                Latt%imj(ii, jj) = Latt%inv_n_list(imjx, imjy)
            enddo
        enddo
        
! 定义晶格基矢和倒格矢
        Latt%a1_v(1) = 1.d0; Latt%a1_v(2) = 0.d0
        Latt%a2_v(1) = 0.d0; Latt%a2_v(2) = 1.d0
        Latt%b1_v(1) = 2.d0 * PI; Latt%b1_v(2) = 0.d0
        Latt%b2_v(1) = 0.d0; Latt%b2_v(2) = 2.d0 * PI
! 定义傅里叶变换的数组（预计算）
        allocate(Latt%xk_v(Lq, 2), Latt%aimj_v(Lq, 2), Latt%k_dot_r(Lq, Lq))
        do ii = 1, Lq
            ix = Latt%n_list(ii, 1)
            iy = Latt%n_list(ii, 2)
            Latt%aimj_v(ii, 1) = dble(ix) * Latt%a1_v(1) + dble(iy) * Latt%a2_v(1)
            Latt%aimj_v(ii, 2) = dble(ix) * Latt%a1_v(2) + dble(iy) * Latt%a2_v(2)
            Latt%xk_v(ii, 1) = dble(ix - 1) * Latt%b1_v(1)/dble(Nlx) + dble(iy - 1) * Latt%b2_v(1)/dble(Nly)
            Latt%xk_v(ii, 2) = dble(ix - 1) * Latt%b1_v(2)/dble(Nlx) + dble(iy - 1) * Latt%b2_v(2)/dble(Nly)
        enddo
        do jj = 1, Lq
            do ii = 1, Lq
                Latt%k_dot_r(ii, jj) = Latt%xk_v(ii, 1)*Latt%aimj_v(jj, 1) + Latt%xk_v(ii, 2)*Latt%aimj_v(jj, 2)
            enddo
        enddo
! 定义单位矩阵，用来计算 sm 更新中的 Delta函数
        allocate(ZKRON(Ndim, Ndim))
        ZKRON = dcmplx(0.d0, 0.d0)
        do ii = 1, Ndim
            ZKRON(ii, ii) = dcmplx(1.d0, 0.d0)
        enddo !define delta function

        allocate(Latt%L_Bonds(Lq, 0:4), Latt%nn_bonds(Lq, 0:4), Latt%LT_bonds(Lq*Ltrot, 0:6))
!define the nearest neighbor bonds
        do iy = 1, Nly
            do ix = 1, Nlx
                ii = Latt%inv_n_list(ix, iy)
                Latt%L_bonds(ii, 0) = ii
                Latt%L_bonds(ii, 1) = Latt%inv_n_list( npbc(ix+1, Nlx), iy )
                Latt%L_bonds(ii, 2) = Latt%inv_n_list( ix, npbc(iy+1, Nly) )
                Latt%L_bonds(ii, 3) = Latt%inv_n_list( npbc(ix-1, Nlx), iy )
                Latt%L_bonds(ii, 4) = Latt%inv_n_list( ix, npbc(iy-1, Nly) )
            enddo
        enddo
! define the second nearest neighbor bonds
        do iy = 1, Nly
            do ix = 1, Nlx
                ii = Latt%inv_n_list(ix, iy)
                Latt%nn_bonds(ii, 0) = ii
                Latt%nn_bonds(ii, 1) = Latt%inv_n_list( npbc(ix+1, Nlx), npbc(iy+1, Nly) )
                Latt%nn_bonds(ii, 2) = Latt%inv_n_list( npbc(ix-1, Nlx), npbc(iy+1, Nly) )
                Latt%nn_bonds(ii, 3) = Latt%inv_n_list( npbc(ix-1, Nlx), npbc(iy-1, Nly) )
                Latt%nn_bonds(ii, 4) = Latt%inv_n_list( npbc(ix+1, Nlx), npbc(iy-1, Nly) )
            enddo
        enddo
! define the nearest neighbors on space-time
        do nt = 1, Ltrot
            do ii = 1, Lq
                iit = Latt%inv_t_list(ii, nt)
                Latt%LT_bonds(iit, 0) = iit
                Latt%LT_bonds(iit, 1) = Latt%inv_t_list(Latt%L_bonds(ii, 1), nt)
                Latt%LT_bonds(iit, 2) = Latt%inv_t_list(Latt%L_bonds(ii, 2), nt)
                Latt%LT_bonds(iit, 3) = Latt%inv_t_list(Latt%L_bonds(ii, 3), nt)
                Latt%LT_bonds(iit, 4) = Latt%inv_t_list(Latt%L_bonds(ii, 4), nt)
                Latt%LT_bonds(iit, 5) = Latt%inv_t_list(ii, npbc(nt+1, Ltrot))
                Latt%LT_bonds(iit, 6) = Latt%inv_t_list(ii, npbc(nt-1, Ltrot))
            enddo
        enddo
	    return
    end subroutine Lattice_make
    
    subroutine Lattice_clear(this)
        type(SquareLattice), intent(inout) :: this
        deallocate(this%n_list, this%inv_n_list, this%t_list, this%inv_t_list)
        deallocate(this%L_bonds, this%nn_bonds, this%LT_bonds, this%imj)
        deallocate(this%xk_v, this%aimj_v, this%k_dot_r)
        deallocate(ZKRON)
        return
    end subroutine Lattice_clear
    
    subroutine FFT_RtoK_1(gr, gk, Latt)
        complex(kind=8), dimension(Lq),      intent(in)       :: gr
        complex(kind=8), dimension(Lq),      intent(out)     :: gk
        type(SquareLattice),                            intent(in)       :: Latt
        integer :: imj, nk
        
        gk = dcmplx(0.d0, 0.d0)
        do imj = 1, Lq
            do nk = 1, Lq
                gk(nk) = gk(nk) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj)
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_1
    
    subroutine FFT_RtoK_2(gr, gk, Latt)
        complex(kind=8), dimension(:,:),    intent(in)    :: gr
        complex(kind=8), dimension(:,:),    intent(out)  :: gk
        type(SquareLattice),                         intent(in)    :: Latt
        integer :: imj, nk, nf, NN2
        
        NN2 = size(gr, 2)
        if ((size(gr, 1) .ne. Lq) .or. (size(gk, 1) .ne. Lq)) then
            write(6,*), 'incorrect matrix size in FFT_RtoK_2'; stop
        endif
        gk = dcmplx(0.d0, 0.d0)
        do nf = 1, NN2
            do imj = 1, Lq
                do nk = 1, Lq
                    gk(nk, nf) = gk(nk, nf) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj, nf)
                enddo
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_2
    
    subroutine FFT_RtoK_3(gr, gk, Latt)
        complex(kind=8), dimension(:,:,:),      intent(in)          :: gr
        complex(kind=8), dimension(:,:,:),      intent(out)        :: gk
        type(SquareLattice),                             intent(in)          :: Latt
        integer :: imj, nk, nf2, nf3, NN2, NN3
        
        NN2 = size(gr, 2); NN3 = size(gr, 3)
        if ((size(gr, 1) .ne. Lq) .or. (size(gk, 1) .ne. Lq)) then
            write(6,*), 'incorrect matrix size in FFT_RtoK_3'; stop
        endif
        gk = dcmplx(0.d0, 0.d0)
        do nf3 = 1, NN3
            do nf2 = 1, NN2
                do imj = 1, Lq
                    do nk = 1, Lq
                        gk(nk, nf2, nf3) = gk(nk, nf2, nf3) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj, nf2, nf3)
                    enddo
                enddo
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_3
    
    subroutine FFT_RtoK_4(gr, gk, Latt)
        complex(kind=8), dimension(:,:,:,:),      intent(in)          :: gr
        complex(kind=8), dimension(:,:,:,:),      intent(out)        :: gk
        type(SquareLattice),                               intent(in)          :: Latt
        integer :: imj, nk, nf2, nf3, nf4, NN2, NN3, NN4
        
        NN2 = size(gr, 2); NN3 = size(gr, 3); NN4 = size(gr, 4)
        if ((size(gr, 1) .ne. Lq) .or. (size(gk, 1) .ne. Lq)) then
            write(6,*), 'incorrect matrix size in FFT_RtoK_4'; stop
        endif
        gk = dcmplx(0.d0, 0.d0)
        do nf4 = 1, NN4
            do nf3 = 1, NN3
                do nf2 = 1, NN2
                    do imj = 1, Lq
                        do nk = 1, Lq
                            gk(nk, nf2, nf3, nf4) = gk(nk, nf2, nf3, nf4) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj, nf2, nf3, nf4)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_4
end module MyLattice



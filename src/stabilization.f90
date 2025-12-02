module Stabilize_mod
    use ProcessMatrix
    use MyMats
    use calc_basic
    implicit none
    
    private
    public :: Wrap_pre, Wrap_L, Wrap_R, Wrap_tau, Stabilize_init, Stabilize_clear, reset_debug_wrap_detail
    
    complex(kind=8) :: Z_one
    complex(kind=8), dimension(:,:), allocatable :: matUDV
    complex(kind=8), dimension(:), allocatable :: TAU, DUP
    integer, dimension(:), allocatable :: IPVT
    complex(kind=8), dimension(:,:), allocatable :: invUUR, invUUL, invVUR, invVUL, temp
    complex(kind=8), dimension(:), allocatable :: bigD
    complex(kind=8), dimension(:,:), allocatable :: mat_big, mat_left, mat_right, Gr_tot, bigU, bigV, invbigU, invbigV
    type(PropGreen), allocatable :: Gr_tmp
    logical, parameter :: wrap_debug_enabled = .false.
    logical, parameter :: wrap_lightlog_enabled = .true.
    integer, save :: debug_wrap_detail_count = 0
    integer, parameter :: debug_wrap_detail_limit = 2000
    integer, save :: debug_stab_green_count = 0
    integer, parameter :: debug_stab_green_limit = 200
    integer, save :: wrap_lightlog_count = 0
    integer, parameter :: wrap_lightlog_limit = 3000
    integer, save :: stab_stage_log_count = 0
    integer, parameter :: stab_stage_log_limit = 3000
    
contains
    subroutine Stabilize_init()
        if (allocated(TAU)) call Stabilize_clear()
        Z_one = dcmplx(1.d0, 0.d0)
        allocate(TAU(Ndim), IPVT(Ndim))
        allocate(DUP(Ndim))
        allocate(matUDV(Ndim, Ndim))
        allocate(invUUR(Ndim, Ndim), invUUL(Ndim, Ndim), invVUR(Ndim, Ndim), invVUL(Ndim, Ndim), temp(Ndim, Ndim))
        allocate(mat_big(2*Ndim, 2*Ndim), mat_left(2*Ndim, 2*Ndim), mat_right(2*Ndim, 2*Ndim))
        allocate(Gr_tot(2*Ndim, 2*Ndim), bigU(2*Ndim, 2*Ndim), bigV(2*Ndim, 2*Ndim), invbigU(2*Ndim, 2*Ndim), invbigV(2*Ndim, 2*Ndim))
        allocate(bigD(2*Ndim))
        allocate(Gr_tmp)
        call Gr_tmp%make()
        return
    end subroutine Stabilize_init
    
    subroutine Stabilize_clear()
        if (allocated(matUDV))  deallocate(matUDV)
        if (allocated(TAU))     deallocate(TAU)
        if (allocated(IPVT))    deallocate(IPVT)
        if (allocated(DUP))     deallocate(DUP)
        if (allocated(invUUR))  deallocate(invUUR)
        if (allocated(invUUL))  deallocate(invUUL)
        if (allocated(invVUR))  deallocate(invVUR)
        if (allocated(invVUL))  deallocate(invVUL)
        if (allocated(temp))    deallocate(temp)
        if (allocated(mat_big))    deallocate(mat_big)
        if (allocated(mat_left))   deallocate(mat_left)
        if (allocated(mat_right))  deallocate(mat_right)
        if (allocated(Gr_tot))     deallocate(Gr_tot)
        if (allocated(bigU))       deallocate(bigU)
        if (allocated(bigV))       deallocate(bigV)
        if (allocated(invbigU))    deallocate(invbigU)
        if (allocated(invbigV))    deallocate(invbigV)
        if (allocated(bigD))       deallocate(bigD)
        if (allocated(Gr_tmp))     deallocate(Gr_tmp)
        return
    end subroutine Stabilize_clear
    
    subroutine QDRP_decompose(Mat, D, IPVT, TAU, WORK, Lwork)
! Arguments: 
        complex(kind=8), dimension(:,:), intent(inout) :: Mat
        complex(kind=8), dimension(:), intent(inout) :: D
        integer, dimension(Ndim), intent(inout) :: IPVT
        complex(kind=8), dimension(Ndim), intent(inout) :: TAU
        complex(kind=8), dimension(:), intent(inout), allocatable :: WORK
        integer, intent(inout) :: Lwork
! Local: 
        real(kind=8), dimension(2*Ndim) :: RWORK
        complex(kind=8) :: WQ
        integer :: info, i, j
        real(kind=8) :: X, eps
! Query optimal amount of memory
        call ZGEQP3(Ndim, Ndim, Mat, Ndim, IPVT, TAU, WQ, -1, RWORK, info)
        Lwork = nint(real(WQ)); allocate(WORK(Lwork))
! QR decomposition of Mat with full column pivoting, Mat * P = Q * R
        call ZGEQP3(Ndim, Ndim, Mat, Ndim, IPVT, TAU, WORK, Lwork, RWORK, info)
! separate off D
        eps = 1d-300
        do i = 1, Ndim
            X = max(abs(Mat(i, i)), eps); D(i) = dcmplx(X, 0.d0) ! protected diagonal entry
            do j = i, Ndim
                Mat(i, j) = Mat(i, j) / X
            enddo
        enddo
        return
    end subroutine QDRP_decompose
    
    subroutine stab_UR(Prop)
        class(Propagator), intent(inout) :: Prop
        integer :: info, i, j, Lwork, Lreq
        complex(kind=8) :: WQ2
        complex(kind=8), dimension(:), allocatable :: WORK
! QR(TMP * U * D) * V
! U*D
        matUDV = Prop%UUR
        do i = 1,Ndim
            matUDV(:, i) = matUDV(:, i) * Prop%DUR(i)
        enddo
! QR(TMP * U * D)
        IPVT = 0
! output in matUDV and DUR
        call QDRP_decompose(matUDV, Prop%DUR, IPVT, TAU, WORK, Lwork)
! Permute V. Since we multiply with V from the right we have to permute the rows of V.
! A V = A P P^-1 V = Q R P^-1 V; apply P^-1 to rows of V (backward permutation)
        call ZLAPMR(.true., Ndim, Ndim, Prop%VUR, Ndim, IPVT) ! lapack 3.3
! V = R * V, output in VUR; character U: only use upper triangular part of matUDV
        call ZTRMM('L', 'U', 'N', 'N', Ndim, Ndim, Z_one, matUDV, Ndim, Prop%VUR, Ndim)
! Generate explicitly U in the previously abused storage of U, output in matUDV
! workspace query for ZUNGQR
        call ZUNGQR(Ndim, Ndim, Ndim, matUDV, Ndim, TAU, WQ2, -1, info)
        Lreq = nint(real(WQ2))
        if (.not. allocated(WORK)) then
            allocate(WORK(Lreq))
        elseif (size(WORK) < Lreq) then
            deallocate(WORK)
            allocate(WORK(Lreq))
        endif
        call ZUNGQR(Ndim, Ndim, Ndim, matUDV, Ndim, TAU, WORK, Lreq, info)
        deallocate(WORK)
        Prop%UUR = matUDV
        return
    end subroutine stab_UR
    
    subroutine  stab_UL(Prop)
        class(Propagator), intent(inout) :: Prop
        integer :: info, i, j, Lwork, Lreq
        complex(kind=8) :: WQ2
        complex(kind=8), dimension(:), allocatable :: WORK
! QR(TMP^\dagger * U^dagger * D) * V^dagger
! U*D
        matUDV = dconjg(transpose(Prop%UUL))
        do i = 1, Ndim
            matUDV(:, i) = matUDV(:, i) * Prop%DUL(i)
        enddo
! QR(TMP^\dagger * U^dagger * D)
        IPVT = 0 
! output in matUDV and DUL
        call QDRP_decompose(matUDV, Prop%DUL, IPVT, TAU, WORK, Lwork)
! Permute V. Since we multiply with V^dagger from the right we have to permute the columns of V.
! A V^dagger = A P P^-1 V^dagger = Q R P^-1 V^dagger; apply P^-1 to columns of V (backward permutation)
        call ZLAPMT(.true., Ndim, Ndim, Prop%VUL, Ndim, IPVT)
! V = V * R^dagger; output in VUL; character U: only use upper triangular part of matUDV
        call ZTRMM('R', 'U', 'C', 'N', Ndim, Ndim, Z_one, matUDV, Ndim, Prop%VUL, Ndim)
! Generate explicitly U in the previously abused storage of U, output in matUDV
! workspace query for ZUNGQR
        call ZUNGQR(Ndim, Ndim, Ndim, matUDV, Ndim, TAU, WQ2, -1, info)
        Lreq = nint(real(WQ2))
        if (.not. allocated(WORK)) then
            allocate(WORK(Lreq))
        elseif (size(WORK) < Lreq) then
            deallocate(WORK)
            allocate(WORK(Lreq))
        endif
        call ZUNGQR(Ndim, Ndim, Ndim, matUDV, Ndim, TAU, WORK, Lreq, info)
        deallocate(WORK)
        Prop%UUL = dconjg(transpose(matUDV))
        return
    end subroutine stab_UL
    
    subroutine stab_green(Gr, Prop, nt)
! Arguments: 
        complex(kind=8), dimension(Ndim, Ndim), intent(out) :: Gr
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: nt
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: VRVL, invULUR
        complex(kind=8), dimension(:), allocatable :: WORK
        integer :: nl, nr, Lwork, info, Lreq
        complex(kind=8) :: WQ2
! coefficient of zgemm: alpha * op(A) * op(B) + beta * op(C); here alpha=Z_one=1, beta=0
! VRVL = VR*VL
        call mmult(VRVL, Prop%VUR, Prop%VUL)
! invULUR = UR^dagger UL^dagger = (UL*UR)^-1; character C: conjugate transpose of UUR and UUL
        call ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_one, Prop%UUR, Ndim, Prop%UUL, Ndim, dcmplx(0.d0, 0.d0), invULUR, Ndim)
! compute: matUDV = (UL*UR)^-1 + DR VR VL DL
        temp = dcmplx(0.d0, 0.d0)
        do nr = 1, Ndim ! temp(1:Ndim, nr) = VRVL(1:Ndim, nr) * DUL(nr)
            call zaxpy(Ndim, Prop%DUL(nr), VRVL(1, nr), 1, temp(1, nr), 1)
        enddo
        VRVL = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim ! VRVL(nl, 1:Ndim) = DUR(nl) * temp(nl, 1:Ndim)
            call zaxpy(Ndim, Prop%DUR(nl), temp(nl, 1), Ndim, VRVL(nl, 1), Ndim)
        enddo
        matUDV = invULUR + VRVL
        if (wrap_debug_enabled) call debug_stab_green_matrix('stab_green_before', nt, matUDV)
        if (nt <= Ltrot/2) then
! 左半区段稳定性更差，直接使用大矩阵公式确保数值稳定
            call stab_green_big(Prop)
            Gr = Gr_tmp%Gr00
            return
        else
! if ntau >Ltrot/2, decompose matUDV^dagger
            matUDV = dconjg(transpose(matUDV))
! matUDV * P  = U D V
            IPVT = 0
            call QDRP_decompose(matUDV, DUP, IPVT, TAU, WORK, Lwork)
            if (wrap_debug_enabled) call debug_stab_green_diag('stab_green_after', nt, DUP)
! matUDV^dagger * P = U D V
! G=UL^dagger *  U * D^-1 * V * P^-1 * UR^dagger; multiply from left to right
            temp = dconjg(transpose(Prop%UUL))
! UL^dagger * U
! workspace query for ZUNMQR (right, no-transpose)
            call ZUNMQR('R', 'N', Ndim, Ndim, Ndim, matUDV, Ndim, TAU, temp, Ndim, WQ2, -1, info)
            Lreq = nint(real(WQ2))
            if (.not. allocated(WORK)) then
                allocate(WORK(Lreq))
            elseif (size(WORK) < Lreq) then
                deallocate(WORK)
                allocate(WORK(Lreq))
            endif
            call ZUNMQR('R', 'N', Ndim, Ndim, Ndim, matUDV, Ndim, TAU, temp, Ndim, WORK, Lreq, info)
! (UL^dagger * U) * D^-1
            do nr = 1, Ndim
                temp(:, nr) = temp(:, nr)/ DUP(nr)
            enddo
! compute (UL^dagger * U * D^-1) * V by solving X * V^dagger = UL^dagger * U * D^-1 * V
            call ZTRSM('R', 'U', 'C', 'N', Ndim, Ndim, Z_one, matUDV(1, 1), Ndim, temp(1, 1), Ndim)
! apply inverse permutation matrix: (UL^dagger * U * D^-1 * V) * P^-1; rearrange columns
            call ZLAPMT(.false., Ndim, Ndim, temp(1, 1), Ndim, IPVT(1))
! (UL^dagger * U * D^-1 * V * P^-1) * UR^dagger = G
! output Gr
            call ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_one, temp(1, 1), Ndim, Prop%UUR(1, 1), Ndim, dcmplx(0.d0, 0.d0), Gr(1, 1), Ndim)
        endif
        deallocate(WORK)
        return
    end subroutine stab_green
        
    subroutine stab_green_big(Prop)
! Arguments:
        class(Propagator), intent(in) :: Prop
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: ULUR, VRVL, invULUR, invVRVL
        complex(kind=8) :: det
        integer :: nl, nr
        
        call mmult(ULUR, Prop%UUL, Prop%UUR)
        call mmult(VRVL, Prop%VUR, Prop%VUL)
        call inv(ULUR, invULUR, det)
        call inv(VRVL, invVRVL, det)
        mat_big = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim
            do nr = 1, Ndim
                mat_big(nl, nr) = invVRVL(nl, nr)
                mat_big(nl + Ndim, nr + Ndim) = invULUR(nl, nr)
            enddo
        enddo
        do nr = 1, Ndim
            mat_big(nr, nr + Ndim) = Prop%DUL(nr)
            mat_big(nr + Ndim, nr) = - Prop%DUR(nr)
        enddo
        call udv(mat_big, bigU, bigD, bigV, 0)
        call inv(bigU, invbigU, det)
        call inv(bigV, invbigV, det)
        call inv(Prop%UUR, invUUR, det)
        call inv(Prop%UUL, invUUL, det)
        call inv(Prop%VUR, invVUR, det)
        call inv(Prop%VUL, invVUL, det)
        mat_big = dcmplx(0.d0, 0.d0); mat_left = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim
            do nr = 1, Ndim
                mat_big(nl, nr) = invVUR(nl, nr)
                mat_big(nl + Ndim, nr + Ndim) =invUUL(nl, nr)
            enddo
        enddo
        call mmult(mat_left, mat_big, invbigV)
        mat_big = dcmplx(0.d0, 0.d0); mat_right = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim
            do nr = 1, Ndim
                mat_big(nl, nr) = invVUL(nl, nr)
                mat_big(nl + Ndim, nr + Ndim) =invUUR(nl, nr)
            enddo
        enddo
        call mmult(mat_right, invbigU, mat_big)
        do nr = 1, 2*Ndim
            do nl = 1, 2*Ndim
                mat_left(nl, nr) = mat_left(nl, nr) / bigD(nr)
            enddo
        enddo
        call mmult(Gr_tot, mat_left, mat_right)
! output time-sliced Green function
        do nl = 1, Ndim
            do nr = 1, Ndim
                Gr_tmp%Gr00(nl, nr) = Gr_tot(nl, nr)
                Gr_tmp%Gr0t(nl, nr) = Gr_tot(nl, nr + Ndim)
                Gr_tmp%Grt0(nl, nr) = Gr_tot(nl + Ndim, nr)
                Gr_tmp%Grtt(nl, nr) = Gr_tot(nl + Ndim, nr + Ndim)
            enddo
        enddo
        return
    end subroutine stab_green_big
    
    real(kind=8) function compare_mat(Gr, Gr2) result(dif)
        complex(kind=8), dimension(Ndim, Ndim), intent(in) :: Gr, Gr2
        dif = maxval(abs(Gr - Gr2))
        return
    end function compare_mat
    
    subroutine Wrap_pre(Prop, WrList, nt)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr
        integer :: nt_st
        if (wrap_debug_enabled) call debug_wrap_detail('Wrap_pre_entry', nt, Prop%Gr)
        if (mod(nt, Nwrap) .ne. 0) then
            write(6,*) "incorrect preortho time slice, NT = ", nt; stop
        endif
        nt_st = nt/Nwrap
        if (nt .ne. 0) call stab_UR(Prop)
        call ensure_matrix_finite('Wrap_pre:Prop%UUR', nt, Prop%UUR)
        WrList%URlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUR(1:Ndim, 1:Ndim)
        WrList%VRlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUR(1:Ndim, 1:Ndim)
        WrList%DRlist(1:Ndim, nt_st) = Prop%DUR(1:Ndim)
        call ensure_matrix_finite('Wrap_pre:URlist', nt, WrList%URlist(:, :, nt_st))
        if (nt == Ltrot) then
            Gr = dcmplx(0.d0, 0.d0)
            call stab_green(Gr, Prop, nt)
            Prop%Gr = Gr
            call ensure_matrix_finite('Wrap_pre:Prop%Gr', nt, Prop%Gr)
            call wrap_lightlog('Wrap_pre', nt, minval(abs(Prop%UUR)), maxval(abs(Prop%UUR)), &
                 minval(abs(WrList%URlist(:, :, nt_st))), maxval(abs(WrList%URlist(:, :, nt_st))), &
                 minval(abs(Prop%Gr)), maxval(abs(Prop%Gr)))
        else
            call wrap_lightlog('Wrap_pre', nt, minval(abs(Prop%UUR)), maxval(abs(Prop%UUR)), &
                 minval(abs(WrList%URlist(:, :, nt_st))), maxval(abs(WrList%URlist(:, :, nt_st))), &
                 0.d0, 0.d0)
        endif
        if (wrap_debug_enabled) call debug_wrap_detail('Wrap_pre_after', nt, Prop%Gr)
        return
    end subroutine Wrap_pre
    
    subroutine Wrap_L(Prop, WrList, nt, flag)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        character(len=*), optional, intent(in) :: flag
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr
        integer :: nt_st
        real(kind=8) :: dif
        if (wrap_debug_enabled) call debug_wrap_detail('Wrap_L_entry', nt, Prop%Gr)
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0) then
            write(6,*) "incorrect ortholeft time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        call ensure_matrix_finite('Wrap_L:Prop%Gr_entry', nt, Prop%Gr)
        Gr = dcmplx(0.d0, 0.d0)
        Prop%UUR(1:Ndim, 1:Ndim) = WrList%URlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUR(1:Ndim, 1:Ndim) = WrList%VRlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUR(1:Ndim) = WrList%DRlist(1:Ndim, nt_st)
        call ensure_matrix_finite('Wrap_L:UR_retrieved', nt, Prop%UUR)
        if (nt == 0) then ! clear URlist
            WrList%URlist = dcmplx(0.d0, 0.d0)
            WrList%VRlist = dcmplx(0.d0, 0.d0)
            WrList%DRlist = dcmplx(0.d0, 0.d0)
        endif
        if (nt .ne. Ltrot) then
            call stab_UL(Prop)
            call stab_green(Gr, Prop, nt)
            dif = compare_mat(Gr, Prop%Gr)
            if (dif > Prop%Xmaxm) Prop%Xmaxm = dif
            if (dif .ge. 5.5d-5) write(6,*) nt, dif, "left ortho unstable in RANK ", IRANK
            if (present(flag)) Prop%Xmeanm = Prop%Xmeanm + dif
            Prop%Gr = Gr
        endif
        call ensure_matrix_finite('Wrap_L:Prop%Gr', nt, Prop%Gr)
        WrList%ULlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUL(1:Ndim, 1:Ndim)
        WrList%VLlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUL(1:Ndim, 1:Ndim)
        WrList%DLlist(1:Ndim, nt_st) = Prop%DUL(1:Ndim)
        call ensure_matrix_finite('Wrap_L:ULlist', nt, WrList%ULlist(:, :, nt_st))
        call wrap_lightlog('Wrap_L', nt, minval(abs(Prop%UUL)), maxval(abs(Prop%UUL)), &
             minval(abs(WrList%ULlist(:, :, nt_st))), maxval(abs(WrList%ULlist(:, :, nt_st))), &
             minval(abs(Prop%Gr)), maxval(abs(Prop%Gr)))
        if (wrap_debug_enabled) call debug_wrap_detail('Wrap_L_after', nt, Prop%Gr)
        return
    end subroutine Wrap_L
    
    subroutine Wrap_R(Prop, WrList, nt, flag)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        character(len=*), optional, intent(in) :: flag
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr
        integer :: nt_st
        real(kind=8) :: dif
        if (wrap_debug_enabled) call debug_wrap_detail('Wrap_R_entry', nt, Prop%Gr)
        call ensure_matrix_finite('Wrap_R:Prop%Gr_entry', nt, Prop%Gr)
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0) then
            write(6,*) "incorrect orthoright time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        Gr = dcmplx(0.d0, 0.d0)
        Prop%UUL(1:Ndim, 1:Ndim) = WrList%ULlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUL(1:Ndim, 1:Ndim) = WrList%VLlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUL(1:Ndim) = WrList%DLlist(1:Ndim, nt_st)
        call ensure_matrix_finite('Wrap_R:UL_retrieved', nt, Prop%UUL)
        if (nt .ne. 0) then
            call stab_stage_log('Wrap_R_before_stabUR', nt, minval(abs(Prop%UUL)), maxval(abs(Prop%UUL)))
            call stab_UR(Prop)
            if (wrap_debug_enabled) call debug_wrap_detail('Wrap_R_post_stabUR', nt, Prop%Gr)
            call stab_stage_log('Wrap_R_after_stabUR', nt, minval(abs(Prop%UUR)), maxval(abs(Prop%UUR)))
            call stab_green(Gr, Prop, nt)
            call stab_stage_log('Wrap_R_after_stabGreen', nt, minval(abs(Gr)), maxval(abs(Gr)))
            dif = compare_mat(Gr, Prop%Gr)
            if (dif > Prop%Xmaxm) Prop%Xmaxm = dif
            if (dif .ge. 5.5d-5) write(6,*) nt, dif, "right ortho unstable in RANK ", IRANK
            if (present(flag)) Prop%Xmeanm = Prop%Xmeanm + dif
            Prop%Gr = Gr
        endif
        call ensure_matrix_finite('Wrap_R:Prop%Gr', nt, Prop%Gr)
        WrList%URlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUR(1:Ndim, 1:Ndim)
        WrList%VRlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUR(1:Ndim, 1:Ndim)
        WrList%DRlist(1:Ndim, nt_st) = Prop%DUR(1:Ndim)
        call ensure_matrix_finite('Wrap_R:URlist', nt, WrList%URlist(:, :, nt_st))
        call wrap_lightlog('Wrap_R', nt, minval(abs(Prop%UUR)), maxval(abs(Prop%UUR)), &
             minval(abs(WrList%URlist(:, :, nt_st))), maxval(abs(WrList%URlist(:, :, nt_st))), &
             minval(abs(Prop%Gr)), maxval(abs(Prop%Gr)))
        if (wrap_debug_enabled) call debug_wrap_detail('Wrap_R_after', nt, Prop%Gr)
        return
    end subroutine Wrap_R

    subroutine Wrap_tau(Prop, PropGr, WrList, nt)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(PropGreen), intent(inout) :: PropGr
        class(WrapList), intent(in) :: WrList
        integer, intent(in) :: nt
! Local: 
        integer :: nt_st
        real(kind=8) :: dif
        if (mod(nt, Nwrap) .ne. 0 .or. nt == 0) then
            write(6,*) "incorrect orthobig time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        Prop%UUL(1:Ndim, 1:Ndim) = WrList%ULlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUL(1:Ndim, 1:Ndim) = WrList%VLlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUL(1:Ndim) = WrList%DLlist(1:Ndim, nt_st)
        call stab_UR(Prop)
        call stab_green_big(Prop)
! stabilization test
        dif = compare_mat(Gr_tmp%Gr0t, PropGr%Gr0t)
        if (dif > PropGr%Xmaxm(1)) PropGr%Xmaxm(1) = dif
        if (dif .ge. 5.5d-5) write(6,*) nt, dif, "GR0T ortho unstable in RANK ", IRANK
        PropGr%Xmeanm(1) = PropGr%Xmeanm(1) + dif
        dif = compare_mat(Gr_tmp%Grt0, PropGr%Grt0)
        if (dif > PropGr%Xmaxm(2)) PropGr%Xmaxm(2) = dif
        if (dif .ge. 5.5d-5) write(6,*) nt, dif, "GRT0 ortho unstable in RANK ", IRANK
        PropGr%Xmeanm(2) = PropGr%Xmeanm(2) + dif
        dif = compare_mat(Gr_tmp%Grtt, PropGr%Grtt)
        if (dif > PropGr%Xmaxm(3)) PropGr%Xmaxm(3) = dif
        if (dif .ge. 5.5d-5) write(6,*) nt, dif, "GRTT ortho unstable in RANK ", IRANK
        PropGr%Xmeanm(3) = PropGr%Xmeanm(3) + dif
        PropGr%Gr00 = Gr_tmp%Gr00
        PropGr%Gr0t = Gr_tmp%Gr0t
        PropGr%Grt0 = Gr_tmp%Grt0
        PropGr%Grtt = Gr_tmp%Grtt
        return
    end subroutine Wrap_tau
    
    subroutine reset_debug_wrap_detail()
        if (.not. wrap_debug_enabled) return
        debug_wrap_detail_count = 0
        debug_stab_green_count = 0
        open(unit=214, file='debug_wrap_detail.log', status='unknown', action='write')
        close(214)
        open(unit=215, file='debug_stab_green.log', status='unknown', action='write')
        close(215)
        return
    end subroutine reset_debug_wrap_detail

    subroutine debug_wrap_detail(tag, nt, Mat)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        complex(kind=8), intent(in) :: Mat(:, :)
        logical :: has_nan
        real(kind=8) :: max_abs
        integer :: r, c
        if (.not. wrap_debug_enabled) return
        if (debug_wrap_detail_count >= debug_wrap_detail_limit) return
        has_nan = .false.
        max_abs = 0.d0
        do r = 1, size(Mat,1)
            do c = 1, size(Mat,2)
                if (.not.(real(Mat(r,c)) == real(Mat(r,c))) .or. .not.(aimag(Mat(r,c)) == aimag(Mat(r,c)))) then
                    has_nan = .true.
                endif
                max_abs = max(max_abs, abs(Mat(r,c)))
            enddo
        enddo
        debug_wrap_detail_count = debug_wrap_detail_count + 1
        open(unit=214, file='debug_wrap_detail.log', status='unknown', action='write', position='append')
        write(214,*) '#', debug_wrap_detail_count, trim(tag), 'nt', nt, 'has_nan', has_nan, 'max_abs', max_abs
        close(214)
        if (has_nan) stop 'debug_wrap_detail detected NaN'
        return
    end subroutine debug_wrap_detail

    subroutine debug_stab_green_matrix(tag, nt, Mat)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        complex(kind=8), intent(in) :: Mat(:, :)
        logical :: has_nan
        real(kind=8) :: max_abs
        integer :: r, c
        if (.not. wrap_debug_enabled) return
        if (debug_stab_green_count >= debug_stab_green_limit) return
        has_nan = .false.
        max_abs = 0.d0
        do r = 1, size(Mat,1)
            do c = 1, size(Mat,2)
                if (.not.(real(Mat(r,c)) == real(Mat(r,c))) .or. .not.(aimag(Mat(r,c)) == aimag(Mat(r,c)))) has_nan = .true.
                max_abs = max(max_abs, abs(Mat(r,c)))
            enddo
        enddo
        debug_stab_green_count = debug_stab_green_count + 1
        open(unit=215, file='debug_stab_green.log', status='unknown', action='write', position='append')
        write(215,*) '#', debug_stab_green_count, trim(tag), 'nt', nt, 'has_nan', has_nan, 'max_abs', max_abs
        close(215)
        if (has_nan) stop 'debug_stab_green matrix NaN'
        return
    end subroutine debug_stab_green_matrix

    subroutine debug_stab_green_diag(tag, nt, diag)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        complex(kind=8), intent(in) :: diag(:)
        logical :: has_nan
        real(kind=8) :: min_abs, max_abs
        integer :: i
        if (.not. wrap_debug_enabled) return
        if (debug_stab_green_count >= debug_stab_green_limit) return
        has_nan = .false.
        min_abs = huge(1.d0)
        max_abs = 0.d0
        do i = 1, size(diag)
            if (.not.(real(diag(i)) == real(diag(i))) .or. .not.(aimag(diag(i)) == aimag(diag(i)))) has_nan = .true.
            min_abs = min(min_abs, abs(diag(i)))
            max_abs = max(max_abs, abs(diag(i)))
        enddo
        debug_stab_green_count = debug_stab_green_count + 1
        open(unit=215, file='debug_stab_green.log', status='unknown', action='write', position='append')
        write(215,*) '#', debug_stab_green_count, trim(tag), 'nt', nt, 'has_nan', has_nan, 'min_abs', min_abs, 'max_abs', max_abs
        close(215)
        if (has_nan) stop 'debug_stab_green diag NaN'
        return
    end subroutine debug_stab_green_diag

    subroutine wrap_lightlog(tag, nt, min_a, max_a, min_b, max_b, min_c, max_c)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        real(kind=8), intent(in) :: min_a, max_a, min_b, max_b, min_c, max_c
        if (.not. wrap_lightlog_enabled) return
        if (wrap_lightlog_count >= wrap_lightlog_limit) return
        wrap_lightlog_count = wrap_lightlog_count + 1
        open(unit=216, file='wrap_light_debug.log', status='unknown', action='write', position='append')
        write(216,'(A16,1X,I8,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6)') &
            trim(tag), nt, min_a, max_a, min_b, max_b, min_c, max_c
        close(216)
        return
    end subroutine wrap_lightlog

    subroutine ensure_matrix_finite(tag, nt, Mat)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        complex(kind=8), intent(in) :: Mat(:, :)
        integer :: r, c
        logical :: bad
        real(kind=8) :: max_abs, min_abs, val
        bad = .false.
        max_abs = 0.d0
        min_abs = huge(0.d0)
        do r = 1, size(Mat,1)
            do c = 1, size(Mat,2)
                val = abs(Mat(r,c))
                max_abs = max(max_abs, val)
                min_abs = min(min_abs, val)
                if (.not.(real(Mat(r,c)) == real(Mat(r,c)))) bad = .true.
                if (.not.(aimag(Mat(r,c)) == aimag(Mat(r,c)))) bad = .true.
                if (bad) exit
            enddo
            if (bad) exit
        enddo
        if (min_abs == huge(0.d0)) min_abs = 0.d0
        if (bad) then
            open(unit=217, file='wrap_nan.log', status='unknown', action='write', position='append')
            write(217,'(A20,1X,I8,1X,ES14.6,1X,ES14.6)') trim(tag), nt, min_abs, max_abs
            close(217)
            stop 'ensure_matrix_finite detected NaN'
        endif
        return
    end subroutine ensure_matrix_finite

    subroutine stab_stage_log(tag, nt, minv, maxv)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nt
        real(kind=8), intent(in) :: minv, maxv
        if (stab_stage_log_count >= stab_stage_log_limit) return
        stab_stage_log_count = stab_stage_log_count + 1
        open(unit=218, file='wrap_stab_debug.log', status='unknown', action='write', position='append')
        write(218,'(A20,1X,I8,1X,ES14.6,1X,ES14.6)') trim(tag), nt, minv, maxv
        close(218)
        return
    end subroutine stab_stage_log

end module Stabilize_mod

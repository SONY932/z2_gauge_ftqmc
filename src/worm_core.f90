module WormCore_mod
    use calc_basic
    use WormProposal_mod
    use MyLattice
    use GaugeFields_mod
    implicit none

contains
    ! ===== 规范作用ΔS累积（三部分：空间+时间+边界）=====
    
    subroutine accumulate_gauge_dS_timeseg(p, bdir, i_src, tau1, tau2, Gauge, Latt)
        ! 为时间段 [tau1, tau2] 累积规范作用量变化
        class(WormProposal), intent(inout) :: p
        character(len=*), intent(in) :: bdir
        integer, intent(in) :: i_src, tau1, tau2
        class(GaugeConf), intent(in) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        integer :: t, tau1m1, tau2p1
        integer :: x, y, xp1, xm1, yp1, ym1
        integer :: i_ll, i_lr, i_ul, i_ur
        integer :: p1, p2
        real(kind=8) :: gam
        logical :: flip0, flipM1
        
        gam = temporal_gamma()
        x = Latt%n_list(i_src, 1)
        y = Latt%n_list(i_src, 2)
        xp1 = npbc(x+1, Nlx); xm1 = npbc(x-1, Nlx)
        yp1 = npbc(y+1, Nly); ym1 = npbc(y-1, Nly)
        
        ! 1) 时间NN项：只有两端边界
        tau1m1 = npbc(tau1-1, Ltrot)
        tau2p1 = npbc(tau2+1, Ltrot)
        
        if (bdir == 'x') then
            p%dS_time = p%dS_time - 2.d0 * gam * dble( &
                Gauge%sigma_x(i_src, tau1) * Gauge%sigma_x(i_src, tau1m1) + &
                Gauge%sigma_x(i_src, tau2) * Gauge%sigma_x(i_src, tau2p1) )
        else
            p%dS_time = p%dS_time - 2.d0 * gam * dble( &
                Gauge%sigma_y(i_src, tau1) * Gauge%sigma_y(i_src, tau1m1) + &
                Gauge%sigma_y(i_src, tau2) * Gauge%sigma_y(i_src, tau2p1) )
        endif
        
        ! 2) 空间plaquette项：逐时片，每片2个plaquette
        do t = tau1, tau2
            if (bdir == 'x') then
                ! σx 影响两个plaquette
                i_ll = Latt%inv_n_list(x, y)
                i_lr = Latt%inv_n_list(xp1, y)
                i_ul = Latt%inv_n_list(x, yp1)
                p1 = p%shadow_sigma_x(Gauge, i_ll, t) * p%shadow_sigma_y(Gauge, i_lr, t) * &
                     p%shadow_sigma_x(Gauge, i_ul, t) * p%shadow_sigma_y(Gauge, i_ll, t)
                
                i_ll = Latt%inv_n_list(x, ym1)
                i_lr = Latt%inv_n_list(xp1, ym1)
                i_ul = Latt%inv_n_list(x, y)
                p2 = p%shadow_sigma_x(Gauge, i_ll, t) * p%shadow_sigma_y(Gauge, i_lr, t) * &
                     p%shadow_sigma_x(Gauge, i_ul, t) * p%shadow_sigma_y(Gauge, i_ll, t)
                
                p%dS_space = p%dS_space + 2.d0 * J * dble(p1 + p2)
            else
                ! σy 影响两个plaquette
                i_ll = Latt%inv_n_list(x, y)
                i_lr = Latt%inv_n_list(xp1, y)
                i_ul = Latt%inv_n_list(x, yp1)
                i_ur = Latt%inv_n_list(xp1, y)
                p1 = p%shadow_sigma_x(Gauge, i_ll, t) * p%shadow_sigma_y(Gauge, i_ll, t) * &
                     p%shadow_sigma_x(Gauge, i_ul, t) * p%shadow_sigma_y(Gauge, i_ur, t)
                
                i_ll = Latt%inv_n_list(xm1, y)
                i_lr = Latt%inv_n_list(x, y)
                i_ul = Latt%inv_n_list(xm1, yp1)
                i_ur = Latt%inv_n_list(xm1, y)
                p2 = p%shadow_sigma_x(Gauge, i_ll, t) * p%shadow_sigma_y(Gauge, i_ll, t) * &
                     p%shadow_sigma_x(Gauge, i_ul, t) * p%shadow_sigma_y(Gauge, i_ur, t)
                
                p%dS_space = p%dS_space + 2.d0 * J * dble(p1 + p2)
            endif
        enddo
        
        ! 3) Gauss边界项（λ相关）：只在跨越τ=1或τ=Ltrot时非零
        flip0 = in_segment(1, tau1, tau2, Ltrot)
        flipM1 = in_segment(Ltrot, tau1, tau2, Ltrot)
        
        if ((flip0 .and. .not. flipM1) .or. (.not. flip0 .and. flipM1)) then
            ! XOR为真：需要加边界项
            call accumulate_boundary_lambda(p, i_src, Gauge, Latt)
        endif
        
        return
    end subroutine accumulate_gauge_dS_timeseg
    
    logical function in_segment(t, tau1, tau2, Ltrot) result(inside)
        integer, intent(in) :: t, tau1, tau2, Ltrot
        if (tau1 <= tau2) then
            inside = (t >= tau1 .and. t <= tau2)
        else
            ! 跨越边界的情况：[tau1..Ltrot, 1..tau2]
            inside = (t >= tau1 .or. t <= tau2)
        endif
        return
    end function in_segment
    
    subroutine accumulate_boundary_lambda(p, i_src, Gauge, Latt)
        use GaugeAction_mod
        class(WormProposal), intent(inout) :: p
        integer, intent(in) :: i_src
        class(GaugeConf), intent(in) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        integer :: s_0, s_M1
        real(kind=8) :: gam
        
        gam = temporal_gamma()
        ! 端点的star乘积（用影子视图）
        s_0 = star_product_shadow(p, Gauge, Latt, i_src, 1)
        s_M1 = star_product_shadow(p, Gauge, Latt, i_src, Ltrot)
        
        p%dS_bdry = p%dS_bdry - 2.d0 * gam * dble(s_0 * Gauge%lambda(i_src) * s_M1)
        return
    end subroutine accumulate_boundary_lambda
    
    integer function star_product_shadow(p, Gauge, Latt, i_site, nt) result(s)
        ! 计算格点的star乘积（使用影子σ）
        class(WormProposal), intent(in) :: p
        class(GaugeConf), intent(in) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: i_site, nt
        integer :: x, y, xm1, ym1, i_left, i_down
        
        x = Latt%n_list(i_site, 1)
        y = Latt%n_list(i_site, 2)
        xm1 = npbc(x-1, Nlx)
        ym1 = npbc(y-1, Nly)
        i_left = Latt%inv_n_list(xm1, y)
        i_down = Latt%inv_n_list(x, ym1)
        
        s = p%shadow_sigma_x(Gauge, i_site, nt) * p%shadow_sigma_y(Gauge, i_site, nt) * &
            p%shadow_sigma_x(Gauge, i_left, nt)  * p%shadow_sigma_y(Gauge, i_down, nt)
        return
    end function star_product_shadow
    
    ! ===== rank-2低秩更新（稳定、可复用）=====
    
    subroutine build_rank2_UV_for_bond_flip(U, V, i, j, sigma_old)
        ! 构造 U,V 用于 O'_b O_b^{-1} 的 rank-2 变换
        ! X = O(σ_new) O(σ_old)^{-1} - I
        complex(kind=8), intent(out) :: U(:,:), V(:,:)  ! (Ndim, 2)
        integer, intent(in) :: i, j, sigma_old
        real(kind=8) :: cval, sval
        
        cval = cosh(Dtau * RT)
        sval = sinh(Dtau * RT)
        
        ! X - I 的矩阵元（只在{i,j}子空间非零）
        ! X - I = [[2s^2, -2cs*σ], [-2cs*σ, 2s^2]]
        U = dcmplx(0.d0, 0.d0)
        V = dcmplx(0.d0, 0.d0)
        
        ! Sherman-Morrison-Woodbury: G' = G - G*U*(I2 + V^T*(I-G)*U)^{-1}*V^T*(I-G)
        ! U的两列在i,j位置
        U(i, 1) = dcmplx(1.d0, 0.d0)
        U(j, 2) = dcmplx(1.d0, 0.d0)
        
        ! V^T = (X-I) = A (参见sigma_rank2.f90)
        V(i, 1) = dcmplx(2.d0*sval*sval, 0.d0)
        V(j, 1) = dcmplx(-2.d0*cval*sval*dble(sigma_old), 0.d0)
        V(i, 2) = dcmplx(-2.d0*cval*sval*dble(sigma_old), 0.d0)
        V(j, 2) = dcmplx(2.d0*sval*sval, 0.d0)
        
        return
    end subroutine build_rank2_UV_for_bond_flip
    
    subroutine rank2_det_and_update(G, U, V, logdet2)
        ! 计算 det(I2 + V^T (I-G) U) 并更新 G
        ! Sherman-Morrison-Woodbury: G' = G - GU(I2 + V^T(I-G)U)^{-1}V^T(I-G)
        complex(kind=8), intent(inout) :: G(:,:)
        complex(kind=8), intent(in) :: U(:,:), V(:,:)
        real(kind=8), intent(out) :: logdet2
        complex(kind=8) :: S(2,2), K(2,2)
        complex(kind=8), allocatable :: VT_ImG(:,:), X(:,:), GU(:,:)
        complex(kind=8) :: detW
        integer :: info, ipiv(2), jj, Ns
        
        Ns = size(G, 1)
        allocate(VT_ImG(2, Ns), X(2, Ns), GU(Ns, 2))
        
        ! K = V^T (I - G) U = V^T U - V^T G U
        K = matmul(transpose(V), U) - matmul(transpose(V), matmul(G, U))
        
        ! S = I2 + K
        S = dcmplx(0.d0, 0.d0)
        S(1,1) = dcmplx(1.d0, 0.d0)
        S(2,2) = dcmplx(1.d0, 0.d0)
        S = S + K
        
        ! logdet2（稳定公式）
        detW = S(1,1)*S(2,2) - S(1,2)*S(2,1)
        logdet2 = log(max(real(detW), 1.d-300))
        
        ! 构造 V^T (I - G) = V^T - V^T * G
        ! 注意：V是(Ns,2)，V^T是(2,Ns)，I-G是(Ns,Ns)
        ! V^T * (I-G) = V^T - V^T*G，结果是(2,Ns)
        VT_ImG = transpose(V) - matmul(transpose(V), G)
        
        ! 解 S * X = V^T (I - G)
        X = VT_ImG
        call zgesv(2, Ns, S, 2, ipiv, X, 2, info)
        if (info /= 0) then
            write(6,*) 'rank2_det_and_update: zgesv failed, info=', info
            logdet2 = -1d30
            deallocate(VT_ImG, X, GU)
            return
        endif
        
        ! 计算 GU
        GU = matmul(G, U)
        
        ! 更新 G ← G - (GU) * X
        G = G - matmul(GU, X)
        
        deallocate(VT_ImG, X, GU)
        return
    end subroutine rank2_det_and_update

    subroutine propagate_shadow_G(p, Gauge, Latt, Bonds, t)
        ! 用影子σ传播格林函数一个时间片（对称Trotter）
        use BondList_mod
        use OpMu_mod
        use MultiplyChecker_mod
        class(WormProposal), intent(inout) :: p
        class(GaugeConf), intent(in) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        integer, intent(in) :: t
        integer :: kk, bidx, ii, jj, sigma
        real(kind=8), parameter :: half = 0.5d0
        
        ! 对称Trotter传播：B = exp(+ε/2 K) exp(-εμ) exp(+ε/2 K)
        ! Gr ← B^{-1} * Gr * B（等时格林函数）
        
        ! 左乘 exp(-εμ)
        call opMu_mmult_L(p%GUtmp, +1)
        call opMu_mmult_L(p%GDtmp, +1)
        
        ! 左乘 exp(+ε/2 K)（使用影子σ）
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 1, 'x', +1, half)
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 2, 'x', +1, half)
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 1, 'y', +1, half)
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 2, 'y', +1, half)
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 2, 'y', +1, half)
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 1, 'y', +1, half)
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 2, 'x', +1, half)
        call apply_shadow_group_L(p, Gauge, Latt, Bonds, t, 1, 'x', +1, half)
        
        ! 右乘 exp(-εμ)^{-1} = exp(+εμ)（注意nflag=+1）
        call opMu_mmult_R(p%GUtmp, -1)
        call opMu_mmult_R(p%GDtmp, -1)
        
        return
    end subroutine propagate_shadow_G
    
    subroutine apply_shadow_group_L(p, Gauge, Latt, Bonds, nt, group_no, dir, nflag, scale)
        ! 左乘一组算符（使用影子σ）
        use OperatorGauge_mod
        use BondList_mod
        class(WormProposal), intent(inout) :: p
        class(GaugeConf), intent(in) :: Gauge
        class(SquareLattice), intent(in) :: Latt
        class(BondList), intent(in) :: Bonds
        integer, intent(in) :: nt, group_no, nflag
        character(len=*), intent(in) :: dir
        real(kind=8), intent(in) :: scale
        type(OperatorGauge) :: Op
        integer :: kk, bidx, ii, jj, sigma
        
        call Op%set(scale)
        
        select case (dir)
        case ('x')
            if (group_no == 1) then
                do kk = 1, size(Bonds%group_ex_even)
                    bidx = Bonds%group_ex_even(kk)
                    ii = Bonds%ex_src(bidx)
                    jj = Bonds%ex_dst(bidx)
                    sigma = p%shadow_sigma_x(Gauge, ii, nt)
                    call Op%mmult_L_2x2(p%GUtmp, ii, jj, sigma, nflag)
                    call Op%mmult_L_2x2(p%GDtmp, ii, jj, sigma, nflag)
                enddo
            else
                do kk = 1, size(Bonds%group_ex_odd)
                    bidx = Bonds%group_ex_odd(kk)
                    ii = Bonds%ex_src(bidx)
                    jj = Bonds%ex_dst(bidx)
                    sigma = p%shadow_sigma_x(Gauge, ii, nt)
                    call Op%mmult_L_2x2(p%GUtmp, ii, jj, sigma, nflag)
                    call Op%mmult_L_2x2(p%GDtmp, ii, jj, sigma, nflag)
                enddo
            endif
        case ('y')
            if (group_no == 1) then
                do kk = 1, size(Bonds%group_ey_even)
                    bidx = Bonds%group_ey_even(kk)
                    ii = Bonds%ey_src(bidx)
                    jj = Bonds%ey_dst(bidx)
                    sigma = p%shadow_sigma_y(Gauge, ii, nt)
                    call Op%mmult_L_2x2(p%GUtmp, ii, jj, sigma, nflag)
                    call Op%mmult_L_2x2(p%GDtmp, ii, jj, sigma, nflag)
                enddo
            else
                do kk = 1, size(Bonds%group_ey_odd)
                    bidx = Bonds%group_ey_odd(kk)
                    ii = Bonds%ey_src(bidx)
                    jj = Bonds%ey_dst(bidx)
                    sigma = p%shadow_sigma_y(Gauge, ii, nt)
                    call Op%mmult_L_2x2(p%GUtmp, ii, jj, sigma, nflag)
                    call Op%mmult_L_2x2(p%GDtmp, ii, jj, sigma, nflag)
                enddo
            endif
        end select
        return
    end subroutine apply_shadow_group_L

end module WormCore_mod


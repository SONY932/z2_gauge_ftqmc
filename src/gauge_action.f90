module GaugeAction_mod
    use calc_basic
    use MyLattice
    use GaugeFields_mod
    implicit none

    public :: delta_S_sigma_flip, delta_S_plaquette_flip, star_product, delta_S_lambda_boundary

contains
    integer function star_product(G, Latt, i_site, nt) result(s)
! 计算格点 i_site 在时片 nt 的"星"乘积：s_r^(ℓ) = ∏_{b∈+r} σ^z_b
! 即该格点4个出射边的乘积：σx(→x), σy(→y), σx(←x), σy(←y)
        class(GaugeConf), intent(in) :: G
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: i_site, nt
        integer :: x, y, xm1, ym1, i_left, i_down
        x = Latt%n_list(i_site, 1)
        y = Latt%n_list(i_site, 2)
        xm1 = npbc(x-1, Nlx)
        ym1 = npbc(y-1, Nly)
        i_left = Latt%inv_n_list(xm1, y)
        i_down = Latt%inv_n_list(x, ym1)
! 4条边：→x, →y, ←x, ←y
        s = G%sigma_x(i_site, nt) * G%sigma_y(i_site, nt) * &
            G%sigma_x(i_left, nt)  * G%sigma_y(i_down, nt)
        return
    end function star_product

    integer function plaq_product(G, Latt, x, y, nt) result(p)
! 约定：以 (x,y) 为左下角的 plaquette：
! 边为 σx(x,y), σy(x+1,y), σx(x,y+1), σy(x,y)
        class(GaugeConf), intent(in) :: G
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: x, y, nt
        integer :: i_ll, i_lr, i_ul
        integer :: xp1, yp1
        xp1 = npbc(x+1, Nlx)
        yp1 = npbc(y+1, Nly)
        i_ll = Latt%inv_n_list(x,   y)
        i_lr = Latt%inv_n_list(xp1, y)
        i_ul = Latt%inv_n_list(x,   yp1)
        p = G%sigma_x(i_ll, nt) * G%sigma_y(i_lr, nt) * G%sigma_x(i_ul, nt) * G%sigma_y(i_ll, nt)
! 注意：上式使用同一时片 nt 的空间边
        return
    end function plaq_product

    real(kind=8) function delta_S_sigma_flip(G, Latt, dir, i_src, nt) result(dS)
! 计算空间 plaquette + 时间耦合（γ）的 ΔS = S_new - S_old
! 输入：
!   dir  ∈ {'x','y'} 表示翻转的是 σx(i,nt) 或 σy(i,nt)
!   i_src: 源格点编号（σx/σy 的起点）
        class(GaugeConf), intent(in) :: G
        class(SquareLattice), intent(in) :: Latt
        character(len=*), intent(in) :: dir
        integer, intent(in) :: i_src, nt
        integer :: x, y
        integer :: p1, p2
        integer :: ym1, xm1
        integer :: ntp1, ntm1
        integer :: s_now, s_p1, s_m1
        real(kind=8) :: gam
        dS = 0.d0
        x = Latt%n_list(i_src, 1)
        y = Latt%n_list(i_src, 2)
        if (dir == 'x') then
! 空间项：翻转 σx 影响两个 plaquette，每个都改变符号
! 哈密顿量 H = -K Σ ∏σ，作用量 S = -εK Σ ∏σ
! 翻转一个σ使∏σ→-∏σ，所以 ΔS = -εK(-∏σ) - (-εK·∏σ) = +2εK·∏σ
! 影响2个plaquette：ΔS_plaq = +2K(p1 + p2)
            p1 = plaq_product(G, Latt, x, y, nt)
            ym1 = npbc(y-1, Nly)
            p2 = plaq_product(G, Latt, x, ym1, nt)
            dS = +2.d0 * J * dble(p1 + p2)
! 时间项：作用量 S = +γ Σ σ(t)σ(t+1)，翻转σ(t)→-σ(t)
! ΔS_time = +γ[-σ(t)(σ(t-1)+σ(t+1))] - γ[σ(t)(σ(t-1)+σ(t+1))]
!         = -2γ σ(t)(σ(t-1)+σ(t+1))
            gam = temporal_gamma()
            ntp1 = npbc(nt+1, Ltrot); ntm1 = npbc(nt-1, Ltrot)
            s_now = G%sigma_x(i_src, nt)
            s_p1  = G%sigma_x(i_src, ntp1)
            s_m1  = G%sigma_x(i_src, ntm1)
            dS = dS - 2.d0 * gam * dble( s_now * (s_p1 + s_m1) )
        elseif (dir == 'y') then
! 空间项：翻转 σy 影响两个 plaquette，ΔS = +2K(p1 + p2)
            p1 = plaq_product(G, Latt, x, y, nt)
            xm1 = npbc(x-1, Nlx)
            p2 = plaq_product(G, Latt, xm1, y, nt)
            dS = +2.d0 * J * dble(p1 + p2)
! 时间项：类似 x 方向
            gam = temporal_gamma()
            ntp1 = npbc(nt+1, Ltrot); ntm1 = npbc(nt-1, Ltrot)
            s_now = G%sigma_y(i_src, nt)
            s_p1  = G%sigma_y(i_src, ntp1)
            s_m1  = G%sigma_y(i_src, ntm1)
            dS = dS - 2.d0 * gam * dble( s_now * (s_p1 + s_m1) )
        else
            write(6,*) 'delta_S_sigma_flip: illegal dir=', dir; stop
        endif
        if (nt == 1 .or. nt == Ltrot) then
            if (dir == 'x') then
                dS = dS + delta_S_lambda_boundary(G, Latt, (/ i_src, Latt%L_bonds(i_src, 1) /))
            else
                dS = dS + delta_S_lambda_boundary(G, Latt, (/ i_src, Latt%L_bonds(i_src, 2) /))
            endif
        endif
        return
    end function delta_S_sigma_flip
    
    real(kind=8) function delta_S_plaquette_flip(G, Latt, x, y, nt) result(dS)
! 计算同时翻转一个plaquette的4条边时的作用量变化
! 空间项：中心plaquette不变，只有4个相邻plaquettes改变符号（上、下、左、右）
! 时间项：4条边各贡献时间耦合
        class(GaugeConf), intent(in) :: G
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: x, y, nt
        integer :: xp1, xm1, yp1, ym1
        integer :: ntp1, ntm1
        integer :: i_ll, i_lr, i_ul, i_ur
        integer :: p
        real(kind=8) :: gam
        
        dS = 0.d0
        xp1 = npbc(x+1, Nlx); xm1 = npbc(x-1, Nlx)
        yp1 = npbc(y+1, Nly); ym1 = npbc(y-1, Nly)
        i_ll = Latt%inv_n_list(x, y)
        i_lr = Latt%inv_n_list(xp1, y)
        i_ul = Latt%inv_n_list(x, yp1)
        i_ur = Latt%inv_n_list(xp1, yp1)
        
! 空间项：中心plaquette不变（4个-1相乘=1）
! 只有4个相邻plaquettes会改变符号（每个只有1条边被翻转）：
! - P(x, y-1)：使用 σx(x,y)（下边的上边）
! - P(x-1, y)：使用 σy(x,y)（左边的右边）
! - P(x+1, y)：使用 σy(x+1,y)（右边的左边）
! - P(x, y+1)：使用 σx(x,y+1)（上边的下边）
! ΔS = +2K Σ p_相邻
        ! 下
        p = plaq_product(G, Latt, x, ym1, nt)
        dS = dS + 2.d0 * J * dble(p)
        ! 左
        p = plaq_product(G, Latt, xm1, y, nt)
        dS = dS + 2.d0 * J * dble(p)
        ! 右
        p = plaq_product(G, Latt, xp1, y, nt)
        dS = dS + 2.d0 * J * dble(p)
        ! 上
        p = plaq_product(G, Latt, x, yp1, nt)
        dS = dS + 2.d0 * J * dble(p)
        
! 时间项：4条边的时间耦合
        gam = temporal_gamma()
        ntp1 = npbc(nt+1, Ltrot); ntm1 = npbc(nt-1, Ltrot)
        
        ! 下边 sigma_x(x,y)
        dS = dS - 2.d0 * gam * dble(G%sigma_x(i_ll, nt) * (G%sigma_x(i_ll, ntm1) + G%sigma_x(i_ll, ntp1)))
        ! 右边 sigma_y(x+1,y)
        dS = dS - 2.d0 * gam * dble(G%sigma_y(i_lr, nt) * (G%sigma_y(i_lr, ntm1) + G%sigma_y(i_lr, ntp1)))
        ! 上边 sigma_x(x,y+1)
        dS = dS - 2.d0 * gam * dble(G%sigma_x(i_ul, nt) * (G%sigma_x(i_ul, ntm1) + G%sigma_x(i_ul, ntp1)))
        ! 左边 sigma_y(x,y)
        dS = dS - 2.d0 * gam * dble(G%sigma_y(i_ll, nt) * (G%sigma_y(i_ll, ntm1) + G%sigma_y(i_ll, ntp1)))
        
        return
    end function delta_S_plaquette_flip

    real(kind=8) function delta_S_lambda_boundary(G, Latt, site_list) result(dS)
! 当边界时片 (τ=1 或 τ=Ltrot) 的 star 乘积翻转时，对 Gauss 边界项的 ΔS
        class(GaugeConf), intent(in) :: G
        class(SquareLattice), intent(in) :: Latt
        integer, intent(in) :: site_list(:)
        integer :: unique_sites(size(site_list))
        integer :: nuniq, idx, site
        integer :: s1, sL
        real(kind=8) :: gam

        gam = temporal_gamma()
        dS = 0.d0
        nuniq = 0
        do idx = 1, size(site_list)
            site = site_list(idx)
            if (site <= 0) cycle
            if (nuniq > 0) then
                if (any(unique_sites(1:nuniq) == site)) cycle
            endif
            nuniq = nuniq + 1
            unique_sites(nuniq) = site
        enddo
        if (nuniq == 0) return
        do idx = 1, nuniq
            site = unique_sites(idx)
            s1 = star_product(G, Latt, site, 1)
            sL = star_product(G, Latt, site, Ltrot)
            dS = dS - 2.d0 * gam * dble(G%lambda(site) * s1 * sL)
        enddo
        return
    end function delta_S_lambda_boundary

end module GaugeAction_mod


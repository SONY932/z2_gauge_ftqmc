program test_trotter
! 测试 apply_half_trotter_R 和 apply_half_trotter_L 是否产生相同的抽象 B 矩阵
    use calc_basic
    use MyLattice
    use BondList_mod
    use GaugeFields_mod
    use LocalSigma_mod
    implicit none
    
    type(SquareLattice) :: Latt
    type(BondList) :: Bonds
    type(GaugeConf) :: Gauge
    complex(kind=8), allocatable :: I_mat(:,:), I_R(:,:), I_L(:,:)
    complex(kind=8), allocatable :: B_R(:,:), B_L(:,:), diff(:,:)
    integer :: ii, nt_test
    real(kind=8) :: max_diff
    integer :: iseed
    
! 初始化参数
    Nlx = 4; Nly = 4; Lq = Nlx * Nly; Ndim = Norb * Lq
    Ltrot = 10; Dtau = 0.1d0; RT = 1.0d0; mu = 0.0d0
    iseed = 12345
    
! 分配数组
    allocate(I_mat(Ndim, Ndim), I_R(Ndim, Ndim), I_L(Ndim, Ndim))
    allocate(B_R(Ndim, Ndim), B_L(Ndim, Ndim), diff(Ndim, Ndim))
    
! 初始化晶格和键
    call Lattice_make(Latt)
    call Bonds%make(Latt)
    call Gauge%make(iseed, Latt)
    
! 设置所有 sigma = +1（简化测试）
    Gauge%sigma_x = 1
    Gauge%sigma_y = 1
    
! 创建单位矩阵
    I_mat = dcmplx(0.d0, 0.d0)
    do ii = 1, Ndim
        I_mat(ii, ii) = dcmplx(1.d0, 0.d0)
    enddo
    
    nt_test = 1
    
! 测试1：apply_half_trotter_R(false) 左乘 vs apply_half_trotter_L(false) 右乘
    I_R = I_mat; I_L = I_mat
    call apply_half_trotter_R(I_R, Latt, Bonds, Gauge, nt_test, reverse=.false.)  ! I_R <- B_R * I = B_R
    call apply_half_trotter_L(I_L, Latt, Bonds, Gauge, nt_test, reverse=.false.)  ! I_L <- I * B_L = B_L
    
    diff = I_R - I_L
    max_diff = maxval(abs(diff))
    write(*,*) 'Test 1: half_trotter_R(false) vs half_trotter_L(false)'
    write(*,*) '  Max difference: ', max_diff
    
! 测试2：apply_half_trotter_R(true) 左乘 vs apply_half_trotter_L(true) 右乘
    I_R = I_mat; I_L = I_mat
    call apply_half_trotter_R(I_R, Latt, Bonds, Gauge, nt_test, reverse=.true.)
    call apply_half_trotter_L(I_L, Latt, Bonds, Gauge, nt_test, reverse=.true.)
    
    diff = I_R - I_L
    max_diff = maxval(abs(diff))
    write(*,*) 'Test 2: half_trotter_R(true) vs half_trotter_L(true)'
    write(*,*) '  Max difference: ', max_diff

! 测试3：完整 Trotter 层比较
! B_R = apply_trotter_layer_R 左乘
! B_L = apply_half_trotter_L(false) + apply_half_trotter_L(true) 右乘
    I_R = I_mat
    call apply_trotter_layer_R(I_R, Latt, Bonds, Gauge, nt_test)  ! I_R <- B_full * I
    
    I_L = I_mat
    call apply_half_trotter_L(I_L, Latt, Bonds, Gauge, nt_test, reverse=.false.)
    call apply_half_trotter_L(I_L, Latt, Bonds, Gauge, nt_test, reverse=.true.)  ! I_L <- I * B_half1 * B_half2
    
    diff = I_R - I_L
    max_diff = maxval(abs(diff))
    write(*,*) 'Test 3: trotter_layer_R vs half_L(false)+half_L(true)'
    write(*,*) '  Max difference: ', max_diff

! 测试4：完整 Trotter 层比较 (两个半步组合)
    I_R = I_mat
    call apply_half_trotter_R(I_R, Latt, Bonds, Gauge, nt_test, reverse=.false.)
    call apply_half_trotter_R(I_R, Latt, Bonds, Gauge, nt_test, reverse=.true.)
    
    I_L = I_mat
    call apply_half_trotter_L(I_L, Latt, Bonds, Gauge, nt_test, reverse=.false.)
    call apply_half_trotter_L(I_L, Latt, Bonds, Gauge, nt_test, reverse=.true.)
    
    diff = I_R - I_L
    max_diff = maxval(abs(diff))
    write(*,*) 'Test 4: half_R(f)+half_R(t) vs half_L(f)+half_L(t)'
    write(*,*) '  Max difference: ', max_diff

    write(*,*) 'Tests completed.'
    
    deallocate(I_mat, I_R, I_L, B_R, B_L, diff)
    
end program test_trotter

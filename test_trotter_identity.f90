program test_trotter_identity
    use calc_basic
    use MyLattice
    use BondList_mod
    use GaugeFields_mod
    use LocalSigma_mod
    implicit none
    
    type(SquareLattice) :: Latt
    type(BondList) :: Bonds
    type(GaugeConf) :: Gauge
    complex(kind=8), allocatable :: Mat_L(:,:), Mat_R(:,:), Identity(:,:), B_gauge(:,:)
    integer :: nt, i
    real(kind=8) :: diff
    
    ! 初始化
    call Latt%init(6, 6)
    call Bonds%init(Latt)
    call Gauge%init(1234)
    
    ! 分配矩阵
    allocate(Mat_L(Ndim, Ndim), Mat_R(Ndim, Ndim), Identity(Ndim, Ndim), B_gauge(Ndim, Ndim))
    
    ! 创建单位矩阵
    Identity = dcmplx(0.d0, 0.d0)
    do i = 1, Ndim
        Identity(i,i) = dcmplx(1.d0, 0.d0)
    enddo
    
    ! 测试：对单位矩阵分别应用 apply_trotter_layer_L 和 apply_trotter_layer_R
    nt = 5
    Mat_L = Identity
    Mat_R = Identity
    
    call apply_trotter_layer_L(Mat_L, Latt, Bonds, Gauge, nt, nflag_in=+1)
    call apply_trotter_layer_R(Mat_R, Latt, Bonds, Gauge, nt, nflag_in=+1)
    
    ! Mat_L = B_gauge * I = B_gauge
    ! Mat_R = I * B_gauge = B_gauge
    ! 它们应该相等
    diff = maxval(abs(Mat_L - Mat_R))
    write(6,*) "Test 1: B_gauge from L and R on Identity"
    write(6,*) "  Max difference:", diff
    
    B_gauge = Mat_R  ! 保存 B_gauge
    
    ! 测试2：验证 B_gauge^{-1} * B_gauge = I
    Mat_L = B_gauge
    call apply_trotter_layer_L(Mat_L, Latt, Bonds, Gauge, nt, nflag_in=-1)  ! B^{-1} * B = ?
    
    diff = maxval(abs(Mat_L - Identity))
    write(6,*) "Test 2: B_gauge^{-1} * B_gauge = I?"
    write(6,*) "  Max difference from I:", diff
    
    ! 测试3：验证 B_gauge * B_gauge^{-1} = I
    Mat_R = B_gauge
    call apply_trotter_layer_R(Mat_R, Latt, Bonds, Gauge, nt, nflag_in=-1)  ! B * B^{-1} = ?
    
    diff = maxval(abs(Mat_R - Identity))
    write(6,*) "Test 3: B_gauge * B_gauge^{-1} = I?"
    write(6,*) "  Max difference from I:", diff
    
    ! 测试4：验证左乘和右乘的一致性（使用非单位矩阵）
    ! 设 M 为随机矩阵，验证 (B * M)^T = M^T * B^T
    ! 即 apply_trotter_layer_L(M) 的结果应该与 transpose(apply_trotter_layer_R(transpose(M))) 一致（如果 B 是对称的）
    write(6,*) "Test 4: Consistency between L and R"
    write(6,*) "  (skipped - requires more setup)"
    
    deallocate(Mat_L, Mat_R, Identity, B_gauge)
    
end program test_trotter_identity

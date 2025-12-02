module WormProposal_mod
    use calc_basic
    implicit none

    type, public :: WormProposal
        ! 提案上下文：影子累计，禁止提案阶段改真场
        logical :: active = .false.
        integer :: bond_idx           ! 空间键索引
        integer :: tau1, tau2         ! 闭区间 [tau1, tau2]
        real(kind=8) :: logRf         ! 单自旋 log(det-ratio) 累计
        real(kind=8) :: dS_space      ! 空间 plaquette 增量
        real(kind=8) :: dS_time       ! 时间 NN 增量
        real(kind=8) :: dS_bdry       ! Gauss 边界项增量
        complex(kind=8), allocatable :: GUtmp(:,:), GDtmp(:,:)  ! 影子 G
        ! 影子场：记录XOR翻转次数（奇数次=翻转）
        integer, allocatable :: flip_count_x(:,:)  ! (Lq, Ltrot)
        integer, allocatable :: flip_count_y(:,:)  ! (Lq, Ltrot)
    contains
        procedure :: begin => worm_begin
        procedure :: commit => worm_commit
        procedure :: discard => worm_discard
        procedure :: shadow_sigma_x => get_shadow_sigma_x
        procedure :: shadow_sigma_y => get_shadow_sigma_y
        procedure :: record_flip => shadow_flip_record
    end type WormProposal

contains
    subroutine worm_begin(this, GrU, GrD)
        class(WormProposal), intent(inout) :: this
        complex(kind=8), intent(in) :: GrU(:,:), GrD(:,:)
        
        this%active = .true.
        if (.not. allocated(this%GUtmp)) then
            allocate(this%GUtmp(Ndim, Ndim))
            allocate(this%GDtmp(Ndim, Ndim))
            allocate(this%flip_count_x(Lq, Ltrot))
            allocate(this%flip_count_y(Lq, Ltrot))
        endif
        
        this%GUtmp = GrU
        this%GDtmp = GrD
        this%logRf = 0.d0
        this%dS_space = 0.d0
        this%dS_time = 0.d0
        this%dS_bdry = 0.d0
        this%flip_count_x = 0
        this%flip_count_y = 0
        return
    end subroutine worm_begin

    subroutine worm_commit(this, GrU, GrD, Gauge)
        use GaugeFields_mod
        class(WormProposal), intent(inout) :: this
        complex(kind=8), intent(inout) :: GrU(:,:), GrD(:,:)
        class(GaugeConf), intent(inout) :: Gauge
        integer :: i, t
        
        if (.not. this%active) return
        
        ! 一次性提交格林函数
        GrU = this%GUtmp
        GrD = this%GDtmp
        
        ! 一次性提交场翻转（XOR语义：奇数次翻转）
        do t = 1, Ltrot
            do i = 1, Lq
                if (mod(this%flip_count_x(i, t), 2) == 1) then
                    Gauge%sigma_x(i, t) = -Gauge%sigma_x(i, t)
                endif
                if (mod(this%flip_count_y(i, t), 2) == 1) then
                    Gauge%sigma_y(i, t) = -Gauge%sigma_y(i, t)
                endif
            enddo
        enddo
        
        this%active = .false.
        return
    end subroutine worm_commit

    subroutine worm_discard(this)
        class(WormProposal), intent(inout) :: this
        ! 什么都不做；影子状态将被丢弃
        this%active = .false.
        return
    end subroutine worm_discard

    integer function get_shadow_sigma_x(this, Gauge, i, t) result(sigma)
        use GaugeFields_mod
        class(WormProposal), intent(in) :: this
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: i, t
        
        sigma = Gauge%sigma_x(i, t)
        if (mod(this%flip_count_x(i, t), 2) == 1) sigma = -sigma
        return
    end function get_shadow_sigma_x

    integer function get_shadow_sigma_y(this, Gauge, i, t) result(sigma)
        use GaugeFields_mod
        class(WormProposal), intent(in) :: this
        class(GaugeConf), intent(in) :: Gauge
        integer, intent(in) :: i, t
        
        sigma = Gauge%sigma_y(i, t)
        if (mod(this%flip_count_y(i, t), 2) == 1) sigma = -sigma
        return
    end function get_shadow_sigma_y

    subroutine shadow_flip_record(this, i, t, dir)
        class(WormProposal), intent(inout) :: this
        integer, intent(in) :: i, t
        character(len=*), intent(in) :: dir
        
        if (dir == 'x') then
            this%flip_count_x(i, t) = this%flip_count_x(i, t) + 1
        elseif (dir == 'y') then
            this%flip_count_y(i, t) = this%flip_count_y(i, t) + 1
        else
            write(6,*) 'shadow_flip_record: illegal dir=', dir
            stop
        endif
        return
    end subroutine shadow_flip_record

end module WormProposal_mod


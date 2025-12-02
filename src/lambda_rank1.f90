module LambdaRank1_mod
    implicit none

contains
    subroutine lambda_flip_rank1(Gr, i, det_ratio, do_update)
! 末片 λ_i 翻转的 rank-1 快更新
! r1 = 2 G_ii - 1，单通道 R_f = |r1|^2；
! G <- G - (-2)/(2G_ii - 1) * (G e_i) * (e_i^T - G[i,:])
        complex(kind=8), intent(inout) :: Gr(:, :)
        integer, intent(in) :: i
        real(kind=8), intent(out) :: det_ratio
        logical, intent(in), optional :: do_update
        logical :: upd
        complex(kind=8) :: r1, factor
        complex(kind=8), dimension(size(Gr,1)) :: Gei
        complex(kind=8), dimension(size(Gr,2)) :: eiT_minus_Grow
        integer :: n

        upd = .true.
        if (present(do_update)) upd = do_update

        r1 = dcmplx(2.d0,0.d0) * Gr(i,i) - dcmplx(1.d0,0.d0)
        det_ratio = real(r1 * dconjg(r1))
        if (det_ratio < 0.d0) det_ratio = 0.d0
        if (.not. upd) return

        if (abs(r1) < 1.d-300) then
            det_ratio = 0.d0  ! 确保设置为0
            return
        endif
        factor = dcmplx(-2.d0,0.d0) / r1
        Gei = Gr(:, i)
        eiT_minus_Grow = dcmplx(0.d0,0.d0)
        do n = 1, size(Gr,2)
            if (n == i) then
                eiT_minus_Grow(n) = dcmplx(1.d0,0.d0) - Gr(i, n)
            else
                eiT_minus_Grow(n) = - Gr(i, n)
            endif
        enddo
        Gr = Gr - factor * matmul( reshape(Gei, [size(Gei,1),1]), reshape(eiT_minus_Grow, [1,size(eiT_minus_Grow,1)]) )
        return
    end subroutine lambda_flip_rank1

end module LambdaRank1_mod



module diag_mod
  use parameters, only: dp
  implicit none

contains

  subroutine print_diag_iz(iz, sumJ_baseline, sumJ_stream)
    integer, intent(in) :: iz
    real(dp), intent(in) :: sumJ_baseline, sumJ_stream
    write(*,'(5x,a,i4,1x,a,2(1x,f12.6))') 'DIAG_iz: iz=', iz, ' sumJ_baseline, sumJ_stream =', sumJ_baseline, sumJ_stream
  end subroutine print_diag_iz

  subroutine print_diag_accum(ia,ja, idim, jdim, delta_block, Gz12, Gz21, tmp1, zstep)
    use parameters, only: dp
    implicit none
    integer, intent(in) :: ia, ja, idim, jdim
    complex(dp), intent(in) :: delta_block(idim,idim)
    complex(dp), intent(in) :: Gz12(idim,jdim)
    complex(dp), intent(in) :: Gz21(jdim,idim)
    complex(dp), intent(in) :: tmp1(idim,idim)
    complex(dp), intent(in) :: zstep
    integer :: ii, jj

    write(*,'(5x,a,i4,a,i4)') 'DIAG_accum: ia=',ia,' ja=',ja
    write(*,'(5x,a)') ' DIAG_accum: delta block (re,im):'
    do ii = 1, idim
      do jj = 1, idim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: delta(',ii,',',jj,')=', real(delta_block(ii,jj)), aimag(delta_block(ii,jj))
      end do
    end do

    write(*,'(5x,a)') ' DIAG_accum: Gz(ia,ja,*,*,spin=2) (re,im):'
    do ii = 1, idim
      do jj = 1, jdim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: Gz(ia,ja)(',ii,',',jj,')=', real(Gz12(ii,jj)), aimag(Gz12(ii,jj))
      end do
    end do

    write(*,'(5x,a)') ' DIAG_accum: Gz(ja,ia,*,*,spin=1) (re,im):'
    do ii = 1, jdim
      do jj = 1, idim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: Gz(ja,ia)(',ii,',',jj,')=', real(Gz21(ii,jj)), aimag(Gz21(ii,jj))
      end do
    end do

    write(*,'(5x,a)') ' DIAG_accum: tmp1 (re,im):'
    do ii = 1, idim
      do jj = 1, idim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: tmp1(',ii,',',jj,')=', real(tmp1(ii,jj)), aimag(tmp1(ii,jj))
      end do
    end do

    write(*,'(5x,a,1x,2(f12.6))') ' DIAG_accum: DIMAG(tmp1*zstep) sample =', DIMAG(tmp1(1,1)*zstep), DIMAG(tmp1(idim,idim)*zstep)

  end subroutine print_diag_accum

end module diag_mod

  subroutine run_diag_iz(diag_iz_max, z, nz, nnnbrs, nblocks, gdim, H, parent, taunew, block_start, block_dim, delta)
    use parameters, only: dp, tpi
    use general, only: nspin
    use green_mod
    implicit none
    integer, intent(in) :: diag_iz_max, nz, nnnbrs, nblocks, gdim
    complex(dp), intent(in) :: z(:)
    complex(dp), intent(in) :: H(:,:,:,:)
    integer, intent(in) :: parent(nnnbrs), block_start(nblocks), block_dim(nblocks)
    real(dp), intent(in) :: taunew(3,nnnbrs)
    complex(dp), intent(in) :: delta(:,:)

    complex(dp), allocatable :: Gtest(:,:,:,:,:,:)
    complex(dp), allocatable :: Gz_tmp(:,:,:,:,:)
    complex(dp), allocatable :: tmp_loc(:,:)
    complex(dp) :: zstep
    real(dp) :: sumJ_baseline, sumJ_stream
    integer :: iz, ia2, ja2, istart, jstart, idim, jdim

    if (diag_iz_max .le. 0) return

    allocate(Gtest(1,nnnbrs,nnnbrs,gdim,gdim,nspin))
    allocate(Gz_tmp(nnnbrs,nnnbrs,gdim,gdim,nspin))

    do iz = 1, min(diag_iz_max, size(z)-1)
      zstep = z(iz+1) - z(iz)
      allocate(tmp_loc(gdim,gdim))

      call compute_g(1,nnnbrs,nblocks,gdim,Gtest,H,(/z(iz)/),parent,taunew,block_start,block_dim)
      call compute_g_onez(nnnbrs,nblocks,gdim,Gz_tmp,H,z(iz),parent,taunew,block_start,block_dim)

      sumJ_baseline = 0.0_dp
      sumJ_stream = 0.0_dp

      do ia2 = 1, nnnbrs
        istart = block_start(parent(ia2))
        idim = block_dim(parent(ia2))
        do ja2 = ia2+1, nnnbrs
          jstart = block_start(parent(ja2))
          jdim = block_dim(parent(ja2))
          if (idim .ne. jdim) cycle
          tmp_loc = cmplx(0.0,0.0,dp)
          tmp_loc(1:idim,1:idim) = MATMUL( MATMUL(delta(istart:istart+idim-1,istart:istart+idim-1), Gtest(1,ia2,ja2,1:idim,1:jdim,2)), &
                                            MATMUL(delta(jstart:jstart+jdim-1,jstart:jstart+jdim-1), Gtest(1,ja2,ia2,1:jdim,1:idim,1)) )
          sumJ_baseline = sumJ_baseline + sum( DIMAG(tmp_loc(1:idim,1:idim)* zstep ) )
          tmp_loc(1:idim,1:idim) = MATMUL( MATMUL(delta(istart:istart+idim-1,istart:istart+idim-1), Gz_tmp(ia2,ja2,1:idim,1:jdim,2)), &
                                            MATMUL(delta(jstart:jstart+jdim-1,jstart:jstart+jdim-1), Gz_tmp(ja2,ia2,1:jdim,1:idim,1)) )
          sumJ_stream = sumJ_stream + sum( DIMAG(tmp_loc(1:idim,1:idim)* zstep ) )
        end do
      end do

      call print_diag_iz(iz, sumJ_baseline, sumJ_stream)
      deallocate(tmp_loc)
    end do

    deallocate(Gtest)
    deallocate(Gz_tmp)

  end subroutine run_diag_iz

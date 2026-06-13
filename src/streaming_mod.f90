module streaming_mod
  use parameters, only: dp
  use general, only: nspin
  use green_mod
  use exchange_utils
  implicit none
contains

  subroutine compute_streaming(nz, z, nnnbrs, nblocks, maxbd, H, parent, taunew, block_start, block_dim, delta, occ, Jorb, Jexc, dbg_ia, dbg_ja)
    implicit none
    integer, intent(in) :: nz, nnnbrs, nblocks, maxbd
    complex(dp), intent(in) :: z(:)
    complex(dp), intent(in) :: H(:,:,:,:)
    integer, intent(in) :: parent(:), block_start(:), block_dim(:)
    real(dp), intent(in) :: taunew(3, :)
    complex(dp), intent(in) :: delta(:,:)
    real(dp), intent(inout) :: occ(:,:,:)
    complex(dp), intent(inout) :: Jorb(:,:,:,:)
    complex(dp), intent(inout) :: Jexc(:,:)
    integer, intent(in) :: dbg_ia, dbg_ja

    complex(dp), allocatable :: Gz(:,:,:,:,:)
    integer :: iz, ia, ja, istart, jstart, iend, jend, idim, jdim
    complex(dp) :: zstep
    integer, allocatable :: idim_idx(:)
    integer :: i

    allocate(Gz(nnnbrs,nnnbrs,maxbd,maxbd,nspin))

    allocate(idim_idx(nnnbrs))
    do i = 1, nnnbrs
      idim_idx(i) = block_dim(parent(i))
    end do

    do iz = 1, nz
      zstep = z(iz+1) - z(iz)
      call compute_g_onez(nnnbrs,nblocks,maxbd,Gz,H,z(iz),parent,taunew,block_start,block_dim)

      do ia = 1, nnnbrs
        istart = block_start(parent(ia))
        idim = block_dim(parent(ia))
        iend = istart + idim - 1
        do ja = ia+1, nnnbrs
          jstart = block_start(parent(ja))
          jdim = block_dim(parent(ja))
          jend = jstart + jdim - 1

          if (idim .ne. jdim) then
            write(*,*) ia,ja,idim,jdim
            stop 'Not equal subblocks size'
          end if

          call accumulate_pair(ia, ja, idim, jdim, delta(istart:iend,istart:iend), delta(jstart:jend,jstart:jend), &
     &                        Gz(ia,ja,1:idim,1:jdim,2), Gz(ja,ia,1:jdim,1:idim,1), zstep, &
     &                        Jorb(ia,ja,1:idim,1:idim), Jexc(ia,ja), dbg_ia, dbg_ja)

        end do
      end do

      call accumulate_occupations(nnnbrs, idim_idx, nspin, Gz, zstep, occ)

    end do

    deallocate(Gz)
    deallocate(idim_idx)

  end subroutine compute_streaming

end module streaming_mod

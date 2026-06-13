module green_mod
  use parameters, only : dp, tpi
  use general, only : hdim, nkp, nspin, xk, wk, atoms, blocks, efermi
  implicit none
contains

  subroutine inverse_complex_matrix(dim,a)
    use parameters, only : dp
    implicit none
    integer :: dim, info
    integer, allocatable :: ipiv(:)
    complex(dp), intent(inout) :: a(dim,dim)
    complex(dp), allocatable :: work(:)

    allocate(ipiv(dim))
    allocate(work(dim))

    call ZGETRF(dim,dim,a,dim,ipiv,info)
    if(info /= 0) stop "inverse_complex_matrix Error in ZGETRF"
    call ZGETRI(dim,a,dim,ipiv,work,dim,info)
    if(info /= 0) stop "inverse_complex_matrix Error in ZGETRI"

    if (allocated(ipiv)) deallocate(ipiv)
    if (allocated(work)) deallocate(work)

  end subroutine inverse_complex_matrix

  SUBROUTINE compute_g_onez(natoms,nblocks,gdim,Gz,H,z,parent,taunew,block_start,block_dim, &
                            debug_print, dbg_ia, dbg_ja, dbg_i, dbg_j, dbg_ispin)

    implicit none

    INTEGER, INTENT(IN) :: natoms, nblocks, gdim, parent(natoms), block_start(nblocks), block_dim(nblocks)
    COMPLEX(DP), INTENT(IN) :: H(:,:,:,:), z
    REAL(DP), INTENT(IN) :: taunew(3,natoms)
    COMPLEX(DP), INTENT(OUT) :: Gz(:,:,:,:,:)  ! assumed-shape for safety

    INTEGER :: i,j,k, ia, ja, ik, ispin, istart, jstart
    LOGICAL, INTENT(IN), OPTIONAL :: debug_print
    INTEGER, INTENT(IN), OPTIONAL :: dbg_ia, dbg_ja, dbg_i, dbg_j, dbg_ispin
    COMPLEX(DP) :: kphase
    COMPLEX(DP), ALLOCATABLE :: Gloc(:,:)

    Gz = cmplx(0.0,0.0,dp)

    allocate(Gloc(hdim,hdim))

    do ispin = 1, nspin
      do ik = 1, nkp
        CALL compute_gloc(Gloc,H(:,:,ik,ispin),z)

        if (present(debug_print) .and. debug_print) then
          if (present(dbg_ia) .and. present(dbg_ja) .and. present(dbg_i) .and. present(dbg_j) .and. present(dbg_ispin)) then
            if (dbg_ispin == ispin) then
              ! compute indices for debug element
              istart = block_start(parent(dbg_ia))-1
              jstart = block_start(parent(dbg_ja))-1
                kphase = cdexp( 1.d0*DCMPLX(0.d0,1.d0)*tpi*&
                  DOT_PRODUCT( xk(:,ik), &
                    ((taunew(:,dbg_ia)-atoms(blocks(parent(dbg_ia))%atom)%pos)-(taunew(:,dbg_ja)-atoms(blocks(parent(dbg_ja))%atom)%pos) ) ))
              write(*,'(a,i6,a,i6,a,i4,a,i4,a,i4)') 'DIAG_onez: ik=',ik,' ia=',dbg_ia,' ja=',dbg_ja,' i=',dbg_i,' j=',dbg_j,' spin=',ispin
              write(*,'(5x,a,2(1x,2(f12.6)))') ' DIAG_onez: Gloc_re,Gloc_im, kphase_re,kphase_im =', real(Gloc(istart+dbg_i,jstart+dbg_j)), aimag(Gloc(istart+dbg_i,jstart+dbg_j)), real(kphase), aimag(kphase)
              write(*,'(5x,a,1x,f12.6)') ' DIAG_onez: wk =', wk(ik)
              write(*,'(5x,a,2(1x,2(f12.6)))') ' DIAG_onez: contrib_re,contrib_im =', real(wk(ik)*Gloc(istart+dbg_i,jstart+dbg_j)*kphase), aimag(wk(ik)*Gloc(istart+dbg_i,jstart+dbg_j)*kphase)
            end if
          end if
        end if

        ! Serial version for debugging (disable OpenMP here)
        DO ia = 1, natoms
          DO ja = 1, natoms
            ! exp(i*k*(Ri-Rj))
                  kphase = cdexp( 1.d0*DCMPLX(0.d0,1.d0)*tpi*&
                  DOT_PRODUCT( xk(:,ik), &
                    ((taunew(:,ia)-atoms(blocks(parent(ia))%atom)%pos)-(taunew(:,ja)-atoms(blocks(parent(ja))%atom)%pos) ) ))

            istart = blocks(parent(ia))%start - 1
            jstart = blocks(parent(ja))%start - 1

            DO i = 1, blocks(parent(ia))%dim
              DO j = 1, blocks(parent(ja))%dim
                Gz(ia,ja,i,j,ispin) = Gz(ia,ja,i,j,ispin) + wk(ik)*Gloc(istart+i,jstart+j)*kphase
              END DO
            END DO

          END DO
        END DO

      END DO
    end do

    deallocate(Gloc)

  END SUBROUTINE compute_g_onez

  SUBROUTINE compute_g(nz,natoms,nblocks,gdim,G,H,z,parent,taunew,block_start,block_dim, &
                       debug_print, dbg_ia, dbg_ja, dbg_i, dbg_j, dbg_ispin)

    implicit none

    INTEGER, INTENT(IN) :: nz, natoms, nblocks, gdim, parent(natoms), block_start(nblocks), block_dim(nblocks)
    COMPLEX(DP), INTENT(IN) :: H(:,:,:,:), z(:)
    REAL(DP), INTENT(IN) :: taunew(3,natoms)
    COMPLEX(DP), INTENT(OUT) :: G(:,:,:,:,:,:)  ! assumed-shape

    INTEGER :: i,j,k, ia, ja, ik, ispin, iz, istart, jstart
    LOGICAL, INTENT(IN), OPTIONAL :: debug_print
    INTEGER, INTENT(IN), OPTIONAL :: dbg_ia, dbg_ja, dbg_i, dbg_j, dbg_ispin
    COMPLEX(DP) :: kphase
    COMPLEX(DP), ALLOCATABLE :: Gloc(:,:)

    G = cmplx(0.0,0.0,dp)

    allocate(Gloc(hdim,hdim))
   
    !$OMP PARALLEL PRIVATE(Gloc,kphase,istart,jstart)
    do ispin = 1, nspin
      DO ik = 1, nkp

        !$OMP DO
        DO iz = 1, nz

          CALL compute_gloc(Gloc,H(:,:,ik,ispin),z(iz))

          if (present(debug_print) .and. debug_print) then
            if (present(dbg_ia) .and. present(dbg_ja) .and. present(dbg_i) .and. present(dbg_j) .and. present(dbg_ispin)) then
              if (dbg_ispin == ispin) then
                istart = blocks(parent(dbg_ia))%start-1
                jstart = blocks(parent(dbg_ja))%start-1
                kphase = cdexp( 1.d0*DCMPLX(0.d0,1.d0)*tpi*&
                    DOT_PRODUCT( xk(:,ik), &
                      ((taunew(:,dbg_ia)-atoms(blocks(parent(dbg_ia))%atom)%pos)-(taunew(:,dbg_ja)-atoms(blocks(parent(dbg_ja))%atom)%pos) ) ))
                write(*,'(a,i6,a,i6,a,i4,a,i4,a,i4,a,i4)') 'DIAG_g: iz=',iz,' ik=',ik,' ia=',dbg_ia,' ja=',dbg_ja,' i=',dbg_i,' j=',dbg_j,' spin=',ispin
                write(*,'(5x,a,2(1x,2(f12.6)))') ' DIAG_g: Gloc_re,Gloc_im, kphase_re,kphase_im =', real(Gloc(istart+dbg_i,jstart+dbg_j)), aimag(Gloc(istart+dbg_i,jstart+dbg_j)), real(kphase), aimag(kphase)
                write(*,'(5x,a,1x,f12.6)') ' DIAG_g: wk =', wk(ik)
                write(*,'(5x,a,2(1x,2(f12.6)))') ' DIAG_g: contrib_re,contrib_im =', real(wk(ik)*Gloc(istart+dbg_i,jstart+dbg_j)*kphase), aimag(wk(ik)*Gloc(istart+dbg_i,jstart+dbg_j)*kphase)
              end if
            end if
          end if

          DO ia = 1, natoms
            DO ja = 1, natoms
              ! exp(i*k*(Ri-Rj))
              kphase = cdexp( 1.d0*DCMPLX(0.d0,1.d0)*tpi*&
                  DOT_PRODUCT( xk(:,ik), &
                    ((taunew(:,ia)-atoms(blocks(parent(ia))%atom)%pos)-(taunew(:,ja)-atoms(blocks(parent(ja))%atom)%pos) ) ))

              istart = blocks(parent(ia))%start-1
              jstart = blocks(parent(ja))%start-1
                            
              DO i = 1, blocks(parent(ia))%dim
                DO j = 1, blocks(parent(ja))%dim
                    G(iz,ia,ja,i,j,ispin) = G(iz,ia,ja,i,j,ispin) + wk(ik)*Gloc(istart+i,jstart+j)*kphase
                END DO
              END DO

            END DO
          END DO

        END DO
        !$OMP END DO

      END DO
    end do
    !$OMP END PARALLEL 

    deallocate(Gloc)

  END SUBROUTINE compute_g

  subroutine compute_gloc(Gloc,H,z)
    implicit none
    complex(dp), intent(in) :: H(:,:), z
    complex(dp), intent(out) :: Gloc(:,:)
    complex(dp), allocatable ::  tmp(:,:)
    integer :: i, j

    Gloc = cmplx(0.0,0.0,dp)

    allocate(tmp(hdim,hdim))
    tmp(:,:) = -1.d0*H(:,:)

    do i=1, hdim
      tmp(i,i) = z + tmp(i,i) + cmplx(efermi,0.d0)
    end do 

    call inverse_complex_matrix(hdim,tmp)
    
    Gloc(:,:) = tmp(:,:)

    deallocate(tmp)

  end subroutine compute_gloc

end module green_mod

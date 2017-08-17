      SUBROUTINE grid_collect_ipe(wwg,zzg,uug,vvg,ttg,rqg,n2g,
     &                   global_lats_a,lonsperlat, lats_nodes_a, kdt)
!!
!! Revision history:
!  2007           Henry Juang, original code
!  2008           Jun Wang  modified buff for write grid component
!  Nov 23 2009    Sarah Lu, comment out 4D tracer
!  Sep 08 2010    Jun Wang  change gfsio to nemsio
!  Dec 16 2010    Jun Wang  change to nemsio library
!  Feb 20 2011    Hann-Ming Henry Juang add code for NDSL
!  Sep 24 2014    S Moorthi - some cleanup and optimization
!  Feb 04 2015    S. Moorthi - threading and optimization
!  May 17 2016    Weiyu Yang - modified from grid_collect.f for WAM-IPE
!                              coupling outputs.
!

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def

      implicit none

      integer, dimension(latg) :: global_lats_a, lonsperlat
      integer, dimension(nodes_comp) :: lats_nodes_a
!
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs) ::  uug,vvg
     &,                                                          ttg,wwg
     &,                                                          zzg,n2g
      real(kind=kind_grid), dimension(lonf,lats_node_a,levh) ::  rqg
      real(kind=kind_io4), dimension(:, :, :), allocatable :: 
     &                                              buff_mult_pieceg_ipe
!
      real(kind=kind_io8), dimension(lonf,lats_node_a) :: buffo, buffi
      integer, dimension(lonf,lats_node_a)             :: kmsk
      integer i, j, k, ngrids_gg_ipe, kdt
!
      ngrids_gg_ipe = 6*levs+levh

      if(.not. allocated(buff_mult_pieceg_ipe)) then
         allocate(buff_mult_pieceg_ipe(lonf,lats_node_a,ngrids_gg_ipe))
      endif
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = wwg(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg_ipe(1,1,k) )
!      write(0,*)'in grid collect, buff_wwg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,k))
      enddo
!!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = zzg(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg_ipe(1,1,levs+k) )
!      write(0,*)'in grid collect, buff_zzg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = uug(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg_ipe(1,1,2*levs+k) )
!      write(0,*)'in grid collect, buff_uug=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,2*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,2*levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = vvg(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg_ipe(1,1,3*levs+k) )
!      write(0,*)'in grid collect, buff_vvg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,3*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,3*levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = ttg(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg_ipe(1,1,4*levs+k) )
!      write(0,*)'in grid collect, buff_ttg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,4*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,4*levs+k))
      enddo
!
      if (levh > 0) then
        do k=1,levh
!$omp parallel do private(i,j)
          do j=1,lats_node_a
            do i=1,lonf
              buffi(i,j) = rqg(i,j,k)
!             if (abs(buffi(i,j)) < 1.0e-15) buffi(i,j) = 0.0
            enddo
          enddo
          CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                     buff_mult_pieceg_ipe(1,1,5*levs+k) )

!      write(0,*)'in grid collect, buff_rqg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+k))
        enddo
      endif
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = n2g(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg_ipe(1,1,5*levs+levh+k) )
!      write(0,*)'in grid collect, buff_n2g=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+levh+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+levh+k))
      enddo

      CALL atmgg_move_ipe(buff_mult_pieceg_ipe,ngrids_gg_ipe, kdt,
     &                    global_lats_a, lats_nodes_a)

      return
      end

      subroutine atmgg_move_ipe(buff_mult_pieceg_ipe,ngrids_gg_ipe, kdt,
     &                          global_lats_a, lats_nodes_a)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_write_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!
      integer ngrids_gg_ipe
      real(kind=kind_io4), dimension(lonf,lats_node_a,ngrids_gg_ipe)
     &                                           :: buff_mult_pieceg_ipe
      real(kind=kind_io4), dimension(lonf,lats_node_a_max,ngrids_gg_ipe)
     &                                           :: grid_node
      real(kind=kind_io4),dimension(:,:,:),  allocatable :: buff_final
      real(kind=kind_io4),dimension(:,:,:,:),allocatable :: grid_nodes
      integer, dimension(latg)       :: global_lats_a
      integer, dimension(nodes_comp) :: lats_nodes_a
      integer ioproc, kdt, lat, ipt_lats
      integer j,k,i,ierr, node
      integer lenrec
!
      ioproc = nodes_comp - 1
      DO k = 1, ngrids_gg_ipe
        DO j = 1, lats_node_a
          DO i = 1, lonf
            grid_node(i, j, k) = buff_mult_pieceg_ipe(i, j, k)
          END DO
        END DO
      END DO
!!
      if(me == ioproc) then
        if(.not. allocated(buff_final)) then
           allocate(buff_final(lonf, latg, ngrids_gg_ipe))
        endif
        if(.not. allocated(grid_nodes)) then
           allocate(grid_nodes(lonf, lats_node_a_max, ngrids_gg_ipe, 
     &                         nodes_comp))
        endif
      else
        if(.not. allocated(grid_nodes)) then
           allocate(grid_nodes(1, 1, 1, 1))
        endif
      endif
!
      if(nodes_comp>1) then
        lenrec = lonf * lats_node_a_max * ngrids_gg_ipe
!
        call mpi_gather( grid_node , lenrec, mpi_real4,
     x                 grid_nodes, lenrec, mpi_real4,
     x                 ioproc, MPI_COMM_ALL, ierr)
      else
        grid_nodes(:,:,:,1)=grid_node(:,:,:)
      endif

      IF(me == ioproc) THEN
        DO k = 1, ngrids_gg_ipe
          ipt_lats = 1
          DO node = 1, nodes_comp
            DO j = 1, lats_nodes_a(node)
              lat = global_lats_a(ipt_lats-1+j)
              DO i = 1, lonf
                buff_final(i, lat, k) = grid_nodes(i, j, k, node)
              END DO
            END DO
            ipt_lats = ipt_lats+lats_nodes_a(node)
          END DO
        END DO
      END IF

      call mpi_barrier(mpi_comm_all,ierr)
      deallocate(grid_nodes)

! Write out the wwg, zzg, uug, vvg, ttg, rqg, n2g full grid fields to
! disk.
!--------------------------------------------------------------------
! buff_final contains wwg, zzg, uug, vvg, ttg, rqg, n2g.
!-------------------------------------------------------
      if(me == ioproc) then
        write(178) kdt, lonf, latg, ngrids_gg_ipe
        print*, 'kdt, lonf, latg, ngrids_gg_ipe=',kdt, lonf, latg,
     &           ngrids_gg_ipe
        write(178) buff_final
        deallocate(buff_final)
      end if
!!
      return
      end

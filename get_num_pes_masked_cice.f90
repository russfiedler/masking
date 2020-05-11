program get_num_pes_masked_cice

!NOTE: The distribution here is ROUNDROBIN. The number of processes won't change but the distribution does for SECTROBIN
!
!

! Figure out how many blocks with ocean need to be used
! just make sure that you run this with ocean_grid.nc here

!Output distribution of tiles
! Also let the user know what the smallest number of ocean points occurring within all tiles is.

use iso_fortran_env
use netcdf
implicit none

integer, parameter :: imt=3600, jmt=2700
real, dimension(imt,jmt) :: ht
real                     :: ht_mv
integer(int16), allocatable, dimension(:,:)     :: proc_rr  ! processor allocated according to round robin 
integer(int16), allocatable, dimension(:,:)     :: num_rr   ! number of ocean points on processor 

integer :: ncid, id_v, ierr

integer :: bx,by, i,j, ib, jb
integer :: nbx, nby

integer :: num, num_good, num_pes

ierr = nf90_open('ocean_grid.nc', nf90_nowrite, ncid)
if ( ierr /= nf90_noerr ) then
   write(*,*) 'Oi, dickhead, where is ocean_grid.nc?'
   stop
endif
ierr = nf90_inq_varid(ncid, 'ht', id_v)
ierr = nf90_get_var(ncid, id_v, ht)
ierr = nf90_get_att(ncid, id_v, 'missing_value', ht_mv)
ierr = nf90_close(ncid)

do
   write(*,*) 'Enter block sizes'
   read(*,*) bx,by
   if (bx < 0 .or. by < 0 ) exit
   if ( mod(imt,bx)/=0 .or. mod(jmt,by)/=0 ) then
      write(*,*) 'Must be divisors of ', imt, jmt
      cycle
   endif

   write(*,*) 'Enter number of processors to be used.'
   read(*,*) num_pes 

   num_good = 0
   nbx = imt/bx
   nby = jmt/by
   num = nbx*nby

   allocate(proc_rr(nbx,nby))
   allocate(num_rr(nbx,nby))
   proc_rr=-1
   num_rr=-1
   
   jb=0
   do j = 1, jmt, by
      jb = jb +1
      ib=0
      do i = 1, imt, bx
         ib =ib+1

         if (any(ht(i:i+bx-1,j:j+by-1) /= ht_mv )) then
            num_good = num_good + 1
            proc_rr(ib,jb) = mod(num_good-1,num_pes)
            num_rr(ib,jb)=count(ht(i:i+bx-1,j:j+by-1) /= ht_mv)
         endif

      enddo
   enddo

   write(*,*) 'Block layout',nbx,nby
   write(*,*) 'Used  ', num_good, ' of ', num
   write(*,*) 'This requires MAXBLOCKS >=', (num_good-1)/num_pes+1
   write(*,*) 'Worst tile has ', minval(num_rr,mask=num_rr>0_int16), ' points'

   call write_blockfile(proc_rr,num_rr,bx,by,num_pes)
   deallocate(proc_rr, num_rr)

enddo

contains

   subroutine write_blockfile(proc,numpts,nx,ny,npes)
   integer(int16), dimension(:,:), intent(in) :: proc, numpts
   integer,                        intent(in) :: nx,ny,npes

   character(23) :: fname
   integer, dimension(2) :: dids
   integer :: ncid, vid,vid1

   fname= 'layout_XBXxXBY_NPES5.nc'

   write(fname(8:10),'(I3.3)') nx
   write(fname(12:14),'(I3.3)') ny
   write(fname(16:20),'(I5.5)') npes

   call handle_err(nf90_create(fname,ior(nf90_netcdf4,nf90_clobber),ncid))
   call handle_err(nf90_def_dim(ncid,'bx',size(proc,dim=1),dids(1)))
   call handle_err(nf90_def_dim(ncid,'by',size(proc,dim=2),dids(2)))
   call handle_err(nf90_def_var(ncid,'process',nf90_short,dids,vid))
   call handle_err(nf90_put_att(ncid,vid,'long_name','allocated processor for round robin distribution'))
   call handle_err(nf90_put_att(ncid,vid,'missing_value',-1_int16))
   call handle_err(nf90_def_var(ncid,'numpts',nf90_short,dids,vid1))
   call handle_err(nf90_put_att(ncid,vid1,'long_name','number of ocean points on pocessor'))
   call handle_err(nf90_put_att(ncid,vid1,'missing_value',-1_int16))
   call handle_err(nf90_enddef(ncid))
   
   call handle_err(nf90_put_var(ncid,vid,proc))
   call handle_err(nf90_put_var(ncid,vid1,numpts))
   
   call handle_err(nf90_close(ncid))

   end subroutine write_blockfile


subroutine handle_err(ierr,isfatal,icode)
integer, intent(in) :: ierr
logical,intent(in),optional :: isfatal
integer,intent(in),optional :: icode

logical :: fatal

fatal = .true.
if ( present(isfatal) ) fatal = isfatal
   
if ( ierr /= NF90_NOERR ) then
   write(*,*) 'ERROR: ', NF90_STRERROR(ierr)
   if ( present(icode) ) then
      write(*,*) 'Code is ',icode
   endif
   if(fatal) stop
endif
end subroutine handle_err

end program get_num_pes_masked_cice

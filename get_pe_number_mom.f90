program get_pe_number_mom
! Get processor masking informatiomn for MOM5
! Enter the size of each tile. Note that this program wants uniform tiles over the domain. We also assume a halo size of 1.


! Figure out how many blocks with ocean need to be used
! just make sure that you run this with ocean_grid.nc here

! MOM is more conservative than I thought. uses halo for some reason by the looks.

use iso_fortran_env
use netcdf
implicit none

integer, parameter :: imt=3600, jmt=2700
real, dimension(0:imt+1,0:jmt+1) :: ht
real, dimension(imt,jmt) :: ht_in
real                     :: ht_mv
integer(int16), allocatable, dimension(:,:)     :: proc_rr  ! processor allocated according to round robin 

integer :: ncid, id_v, ierr

integer :: bx,by, i,j, ib, jb
integer :: nbx, nby

integer :: num, num_good, num_pes

! Get grid
ierr = nf90_open('ocean_grid.nc', nf90_nowrite, ncid)
if ( ierr /= nf90_noerr ) then
   write(*,*) 'Oi, dickhead, where is ocean_grid.nc?'
   stop
endif
ierr = nf90_inq_varid(ncid, 'ht', id_v)
ierr = nf90_get_var(ncid, id_v, ht_in)
ierr = nf90_get_att(ncid, id_v, 'missing_value', ht_mv)
ierr = nf90_close(ncid)

! Fill halos

ht(1:imt,1:jmt)=ht_in
ht(:,0)=ht_mv
ht(0,jmt+1)=ht_mv
ht(imt+1,jmt+1)=ht_mv
ht(0,1:jmt)=ht_in(imt,:)
ht(imt+1,1:jmt)=ht_in(1,:)
ht(:,jmt+1)=ht(imt+1:0:-1,jmt)

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
   proc_rr=-1
   
   jb=0
   do j = 1, jmt, by
      jb = jb +1
      ib=0
      do i = 1, imt, bx
         ib =ib+1

         !if (any(ht(i:i+bx-1,j:j+by-1) /= ht_mv )) then
         if (any(ht(i-1:i+bx,j-1:j+by) /= ht_mv )) then
            proc_rr(ib,jb) = num_good
            num_good = num_good + 1
            if(all(ht(i:i+bx-1,j:j+by-1) == ht_mv)) then
              write(*,*) 'Block',ib,jb, num_good-1, ' is all land but not masked not masked due to halos',i,j
              if(any(ht(i:i+bx-1,j-1) /= ht_mv)) write(*,*) 'S'
              if(any(ht(i+bx,j:j+by-1) /= ht_mv)) write(*,*) 'E'
              if(any(ht(i:i+bx-1,j+by) /= ht_mv)) write(*,*) 'N'
              if(any(ht(i-1,j:j+by-1) /= ht_mv)) write(*,*) 'W'

              if((ht(i-1,j-1) /= ht_mv)) write(*,*) 'SW'
              if((ht(i+bx,j-1) /= ht_mv)) write(*,*) 'SE'
              if((ht(i+bx,j+by) /= ht_mv)) write(*,*) 'NE'
              if((ht(i-1,j+by) /= ht_mv)) write(*,*) 'NW'
            endif
         endif

      enddo
   enddo

   write(*,*) 'Block layout',nbx,nby
   write(*,*) 'Used  ', num_good, ' of ', num
   write(*,*) 'This requires MAXBLOCKS >=', (num_good-1)/num_pes+1

   call write_blockfile(proc_rr,bx,by,num_pes)
   deallocate(proc_rr)

enddo

contains

   subroutine write_blockfile(proc,nx,ny,npes)
   integer(int16), dimension(:,:), intent(in) :: proc
   integer,                        intent(in) :: nx,ny,npes

   character(23) :: fname
   integer, dimension(2) :: dids
   integer :: ncid, vid

   fname= 'procno_XBXxXBY_NPES5.nc'

   write(fname(8:10),'(I3.3)') nx
   write(fname(12:14),'(I3.3)') ny
   write(fname(16:20),'(I5.5)') npes

   call handle_err(nf90_create(fname,ior(nf90_netcdf4,nf90_clobber),ncid))
   call handle_err(nf90_def_dim(ncid,'bx',size(proc,dim=1),dids(1)))
   call handle_err(nf90_def_dim(ncid,'by',size(proc,dim=2),dids(2)))
   call handle_err(nf90_def_var(ncid,'process',nf90_short,dids,vid))
   call handle_err(nf90_put_att(ncid,vid,'long_name','allocated processor for MOM landmasking'))
   call handle_err(nf90_put_att(ncid,vid,'missing_value',-1_int16))
   call handle_err(nf90_enddef(ncid))
   
   call handle_err(nf90_put_var(ncid,vid,proc))
   
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

end program get_pe_number_mom

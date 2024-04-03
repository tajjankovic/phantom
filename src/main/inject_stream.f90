!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!

!
! Stores simulation data to avoid reading it at every timestep
!
! :References: None
module store_simulation_data
    implicit none
    ! Parameters
    integer :: np
    real :: gamma,mass,vinj1,vinj2,rstream1,zstream1,rstream2,zstream2,&
            inclination,yshift,offset,dtinject_set
    ! Particle arrays
    real, allocatable :: xyzh_cyl1(:,:), vxyzu_cyl1(:,:),xyzh_cyl2(:,:), vxyzu_cyl2(:,:)

    ! Status flag to indicate if data has been initialized
    logical :: is_initialized = .false.
    logical :: identical_streams = .false.

end module store_simulation_data



module inject
!
! Handles stream injection for simulations of the stream self-crossing in TDEs
!
! :References: None
!
! :Owner: Taj Jankovič & Daniel Price
!
! :Runtime parameters:
!   -
!
! :Dependencies: boundary, eos, infile_utils, io, part, partinject,
!   physcon, units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'selfcrossing'
 public :: init_inject,inject_particles,write_options_inject,read_options_inject


contains




!-----------------------------------------------------------------------
!+
!  Subroutine to append an element to a list
!+
!-----------------------------------------------------------------------
subroutine AddToList(list, element)

          IMPLICIT NONE

          integer :: i, isize
          integer, intent(in) :: element
          integer, dimension(:), allocatable, intent(inout) :: list
          integer, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize
              clist(i) = list(i)
              end do
              clist(isize+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(1))
              list(1) = element
          end if


end subroutine AddToList

!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0
end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling particle injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use part,      only:igas,massoftype,kill_particle,shuffle_part
 use partinject,only:add_or_update_particle
 use io,        only:master,iprint
 use sortutils, only:indexx
 use vectorutils,   only:rotatevec
 use physcon,   only:pi,solarm,years,c
 use random,    only:ran2
 use store_simulation_data, only: np,gamma,mass,vinj1,rstream1,zstream1,&
                                  inclination,yshift,offset,dtinject_set,&
                                  is_initialized,identical_streams
 use store_simulation_data, only: vinj2,rstream2,zstream2
 use store_simulation_data, only: xyzh_cyl1,vxyzu_cyl1,xyzh_cyl2,vxyzu_cyl2
 integer, dimension(:), allocatable :: list1,list2
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real,    allocatable   :: xyzh_old(:,:),vxyzu_old(:,:)
 integer,    allocatable   :: list_temp(:)
 real :: y0up,inc
 real :: y0,time2,dz
 integer :: i,nskip,ninj1,ninj2,ierr
 logical :: setexists,store_previous_run
 character(len=512)::setupfile,cylinder_file1,cylinder_file2

 store_previous_run = .false.
 ninj1=0
 ninj2=0

 ! Initialize variables from the setup file
 if (.not. is_initialized) then
     ! store previous values if this is a re-run
     if (npart>0) then
       allocate(xyzh_old(4,npart),stat=ierr)  ! positions + smoothing length
       ! read_dump will overwrite the current particles, so store them in a temporary array
       allocate(vxyzu_old(4,npart),stat=ierr)  ! positions + smoothing length
       if (ierr /= 0) stop ' error allocating memory to store positions'
       !allocate(vxyzu2(maxvxyzu,npart),stat=ierr)  ! velocity + thermal energy
       !if (ierr /= 0) stop ' error allocating memory to store velocity'
       !allocate(eos_vars2(maxeosvars,maxp),stat=ierr) ! temperature
       !if (ierr /= 0) stop ' error allocating memory to store temperature'
       xyzh_old  = xyzh
       vxyzu_old = vxyzu
       store_previous_run = .true.
     endif

     ! Read setup parameters only at the beginning
     setupfile = 'stream.setup'

     inquire(file=setupfile,exist=setexists)
     if (setexists) then
        call read_setup_options_inject(setupfile,np,gamma,vinj1,vinj2,rstream1,zstream1,&
                                        rstream2,zstream2,inclination,yshift,offset,dtinject_set,&
                                        cylinder_file1,cylinder_file2,ierr)
        if (ierr /= 0) then
           stop 'please rerun phantomsetup with revised .setup file'
        endif
        !--Prompt to get inputs and write to file
     else
        print "(a,/)",trim(setupfile)//' not found'
        stop 'please check and edit .setup file and rerun phantom'
     endif

     ! Read and store particle data from the ASCII file
     nskip = 12  ! Number of header lines to skip
     !cylinder_file1 = "cylinder1.ascii"
     call read_stream_ascii(cylinder_file1,nskip,mass,xyzh_cyl1,vxyzu_cyl1)

     !cylinder_file2 = "cylinder1.ascii"
     call read_stream_ascii(cylinder_file1,nskip,mass,xyzh_cyl2,vxyzu_cyl2)

     is_initialized = .true.  ! Set flag to true after initial setup
     if ((cylinder_file1==cylinder_file2) .and. (vinj1==vinj2) .and. (rstream1==rstream2))then
       identical_streams = .true.
     endif
     print *,'Time step and velocity:',dtlast,vinj1
     print *,'Stream offset:',dz
     print *,'Particle mass:',massoftype(igas)
 endif

 !y0 = yshift!stream is from y=[y0-H,y0+H]
 !y0up = y0 +zstream1
 dz = offset*rstream1
 dtinject = dtinject_set
 inc=pi/180.*(inclination/2.) !angle of streams for x-axis

 !advance the first cylinder
 allocate(list_temp(np), stat=ierr)
 call advance_cylinder(list_temp,time,dtlast,np,yshift,rstream1,zstream1,vinj1,xyzh_cyl1(2,:),ninj1)
 do i =1,ninj1
   call AddToList(list1,list_temp(i))
 end do
 deallocate(list_temp)

 !inject "upper" stream
 if (ninj1>1e-2) then
   call inject_streams(list1,npart,time,dtlast,vinj1,dz,inc,ninj1,xyzh_cyl1,&
                           vxyzu_cyl1,npartoftype, xyzh, vxyzu,identical_streams,.true.)
 endif

 !advance the second cylinder if streams are not identical
 if (.not. identical_streams) then
   allocate(list_temp(np), stat=ierr)
   call advance_cylinder(list_temp,time,dtlast,np,yshift,rstream2,zstream2,vinj2,xyzh_cyl2(2,:),ninj2)
   do i =1,ninj2
     call AddToList(list2,list_temp(i))
   end do
   deallocate(list_temp)

   !inject "lower" stream
   if (ninj2>1e-2) then
     call inject_streams(list2,npart,time,dtlast,vinj2,dz,inc,ninj2,xyzh_cyl2,&
                             vxyzu_cyl2,npartoftype, xyzh, vxyzu,identical_streams,.false.)
   endif
 endif



 !
 !-- timestep constraint
 !
 !if (time>0.) deallocate(list)!always 5 columns - xyzmh

end subroutine inject_particles




!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit


end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 character(len=30), parameter :: label = 'read_options_inject'

end subroutine read_options_inject


!-----------------------------------------------------------------------
!+
!  Reads input options from the setup file
!+
!-----------------------------------------------------------------------
subroutine read_setup_options_inject(filename,np,gamma,vinj1,vinj2,rstream1,zstream1,&
                                      rstream2,zstream2,inclination,yshift,offset,dtinject,&
                                      inputfile1,inputfile2,ierr)
  use infile_utils,  only:open_db_from_file,close_db,inopts,read_inopt
  !use setunits,      only:read_options_and_set_units
  !use units,         only:umass,udist,unit_luminosity
  !use physcon,       only:solarm,solarr,solarl
  use part,           only:gr
  character(len=*), intent(in)  :: filename
  character(len=*), intent(inout)  :: inputfile1,inputfile2
  !integer,                   intent(inout) :: ieos
  real,                      intent(inout) :: gamma,inclination,yshift,vinj1,vinj2
  real,                      intent(inout) :: rstream1,rstream2,zstream1,zstream2,offset,dtinject
  integer,                      intent(inout) :: np

  integer,          intent(out) :: ierr
  integer                       :: nerr
  integer,          parameter   :: lu = 21
  type(inopts), allocatable     :: db(:)
  !real :: mstream_msun,rstream_rsun

  call open_db_from_file(db,filename,lu,ierr)
  if (ierr /= 0) return

  nerr = 0

  call read_inopt(np,'np',db,errcount=nerr)
  call read_inopt(gamma,'gamma',db,errcount=nerr)
  call read_inopt(rstream1,'rstream1',db,errcount=nerr)
  call read_inopt(zstream1,'zstream1',db,errcount=nerr)
  call read_inopt(rstream2,'rstream2',db,errcount=nerr)
  call read_inopt(zstream2,'zstream2',db,errcount=nerr)
  call read_inopt(inclination,'inclination',db,errcount=nerr)
  call read_inopt(yshift,'yshift',db,errcount=nerr)
  call read_inopt(vinj1,'vinj1',db,errcount=nerr)
  call read_inopt(vinj2,'vinj2',db,errcount=nerr)
  call read_inopt(offset,'offset',db,errcount=nerr)
  call read_inopt(dtinject,'dtinject',db,errcount=nerr)
  call read_inopt(inputfile1,'inputfile1',db,errcount=nerr)
  call read_inopt(inputfile2,'inputfile2',db,errcount=nerr)

  !call read_inopt(gamma,'gamma',db,errcount=nerr)


  if (nerr /= 0) ierr = ierr + 1


  if (nerr > 0) then
     print "(1x,a,i2,a)",'setup_stream: ',nerr,' error(s) during read of setup file'
     ierr = 1
  endif

  call close_db(db)

end subroutine read_setup_options_inject




!-----------------------------------------------------------------------
!+
!  Reads cylinder properties from an .ascii file
!+
!-----------------------------------------------------------------------
subroutine read_stream_ascii(filename,nskip,mass,xyzh_cyl,vxyzu_cyl)
  character(len=512),    intent(in)::filename
  integer, intent(in) :: nskip
  real, allocatable, intent(out) :: xyzh_cyl(:,:), vxyzu_cyl(:,:)
  real,    intent(inout) :: mass
  real, allocatable :: dataold(:,:)
  integer :: i,n,io,unit

  io = 0
  unit = 10

  open(unit, file=filename,status="old",action="read")
  n = 0
  do
     read(unit,*,iostat=io)
     IF (io/=0) EXIT
     n = n+1
  end do
  rewind(unit) ! Rewind the file to read the data

  ! Skip header lines
  do i = 1, nskip
    read(unit, *)  ! Read and discard header lines
  end do

  ! Allocate based on total lines minus header
  n=n-nskip
  allocate(dataold(5, n))

  ! Allocate arrays based on determined size
  allocate(xyzh_cyl(4, n))
  allocate(vxyzu_cyl(4, n))

  ! Now, read the actual data starting after the header
  do i = 1, n
     read(unit, *) dataold(:, i)
  end do

  ! Assign data to output arrays
  xyzh_cyl = dataold([1, 2, 3, 5], :)
  vxyzu_cyl(4, :) = dataold(10, :)  ! Assuming thermal energy or another property is in the 10th column
  mass = dataold(4, 1)  ! Assuming all particles have the same mass

  close(unit)
  deallocate(dataold)!always 4 columns - xyzh

end subroutine read_stream_ascii


!
! rearrange particles in the cylinder and determine how many are injected
!

subroutine advance_cylinder(list_temp,time,dtlast,np,y0,rstream,zstream,vinj,ycyl,ninj)
  integer, intent(in) :: np
  integer, intent(inout) :: list_temp(:)
  real,    intent(inout) :: ycyl(:)
  real,    intent(in) :: rstream,zstream,vinj,time,dtlast,y0
  integer,    intent(inout) :: ninj
  real  :: tperiod,tcross,shift,y0up
  integer :: i,n,io

  y0up = y0 +zstream ! cylinder "top"
  ninj=0
  tperiod = zstream/vinj
  do i = 1,np
    shift = vinj*time
    tcross = (ycyl(i) -y0)/vinj

    if ((ycyl(i) - shift)<y0) then
      ycyl(i) = y0up - vinj*(time - tcross - floor((time-tcross)/tperiod)*tperiod)
    else
      ycyl(i) = ycyl(i) - shift
    endif

    if ((ycyl(i) - vinj*dtlast)<y0) then
      !call AddToList(list,i)
      ninj = ninj+1
      list_temp(ninj) = i !maybe here?
    endif
  end do
end subroutine advance_cylinder




!
! Inject stream
!
subroutine inject_streams(list,npart,time,dtlast,vinj,dz,inc,ninj,xyzh_cyl,&
                        vxyzu_cyl, npartoftype, xyzh, vxyzu,identical_streams,upper_stream)
  use part,      only:igas
  use physcon,   only:pi
  use vectorutils,   only:rotatevec
  use partinject,only:add_or_update_particle
  logical, intent(in) :: upper_stream,identical_streams
  integer, intent(in) :: ninj
  integer, intent(inout) :: list(:),npartoftype(:)
  real,    intent(inout) :: xyzh_cyl(:,:),vxyzu_cyl(:,:),xyzh(:,:),vxyzu(:,:)
  real,    intent(in) :: vinj,dz,time,dtlast,inc
  integer,    intent(inout) :: npart
  real  :: xyzi1(3),xyzi2(3),vxyz1(3),vxyz2(3),u
  real,allocatable  :: xcyl(:),ycyl(:),zcyl(:),hcyl(:),ucyl(:)
  integer :: i,i_part

  xcyl = xyzh_cyl(1,:)
  ycyl = xyzh_cyl(2,:)
  zcyl = xyzh_cyl(3,:)
  hcyl = xyzh_cyl(4,:)
  !ucyl1 = vxyzu_cyl1(4,:)

  u = 1.e-5*vinj**2. ! change here! - u later
  do i=1,ninj
    if (identical_streams) then
       !upper stream
       xyzi1 = (/xcyl(list(i)),ycyl(list(i))-vinj*dtlast,zcyl(list(i))+ dz/2./)
       vxyz1 = (/0.,-vinj,0./)
       !u = ucyl(list(i))
       call rotatevec(xyzi1,(/0.,0.,1./),pi/2.-inc)
       call rotatevec(vxyz1,(/0.,0.,1./),pi/2.-inc)
       i_part = npart +1
       call add_or_update_particle(igas, xyzi1, vxyz1, hcyl(list(i)),u, i_part, npart, npartoftype, xyzh, vxyzu)

       !lower stream
       xyzi2 = (/xcyl(list(i)),-ycyl(list(i))+vinj*dtlast,zcyl(list(i)) - dz/2./)
       vxyz2 = (/0.,vinj,0./)
       call rotatevec(xyzi2,(/0.,0.,1./),-(pi/2.-inc))
       call rotatevec(vxyz2,(/0.,0.,1./),-(pi/2.-inc))
       i_part = npart +1
       call add_or_update_particle(igas, xyzi2, vxyz2, hcyl(list(i)),u, i_part, npart, npartoftype, xyzh, vxyzu)
     elseif (upper_stream) then
       !upper stream
       xyzi1 = (/xcyl(list(i)),ycyl(list(i))-vinj*dtlast,zcyl(list(i))+ dz/2./)
       vxyz1 = (/0.,-vinj,0./)
       !u = ucyl(list(i))
       call rotatevec(xyzi1,(/0.,0.,1./),pi/2.-inc)
       call rotatevec(vxyz1,(/0.,0.,1./),pi/2.-inc)
       i_part = npart +1
       call add_or_update_particle(igas, xyzi1, vxyz1, hcyl(list(i)),u, i_part, npart, npartoftype, xyzh, vxyzu)
     elseif (.not. upper_stream) then
       !lower stream
       xyzi1 = (/xcyl(list(i)),-ycyl(list(i))+vinj*dtlast,zcyl(list(i)) - dz/2./)
       vxyz1 = (/0.,vinj,0./)
       !u = ucyl(list(i))
       call rotatevec(xyzi1,(/0.,0.,1./),-(pi/2.-inc))
       call rotatevec(vxyz1,(/0.,0.,1./),-(pi/2.-inc))
       i_part = npart +1
       call add_or_update_particle(igas, xyzi1, vxyz1, hcyl(list(i)),u, i_part, npart, npartoftype, xyzh, vxyzu)
     endif
  enddo
end subroutine inject_streams


end module inject

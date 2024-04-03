!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! this module does setup
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon, units
!
 use setstream,            only:stream_t
 use part,           only:gravity,gr
 use eos,                only:ieos

 implicit none
 public :: setpart
 type(stream_t)       :: stream

 private

contains

!----------------------------------------------------------------
!+
!  empty setup for driven simulation
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use units,   only:set_units
 use physcon, only:au,solarm
 use timestep,       only:tmax,dtmax
 use io,             only:master,fatal
 use eos,             only:init_eos,finish_eos,gmw
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer                          :: ierr
 logical            :: iexist,setexists
 character(len=120)               :: setupfile,inname

 !call set_units(dist=100.*au,mass=solarm,G=1.)
 time = 0.
 polyk = 0.
 gamma = 4./3.

 npart = 0
 npartoftype(:) = 0
 massoftype = 1.

 !
 ! determine if the .in file exists
 !
 inname=trim(fileprefix)//'.in'
 inquire(file=inname,exist=iexist)
 if (.not. iexist) then
    tmax  = 100.
    dtmax = 1.0
    ieos  = 2
 endif
 !

 !
 ! determine if the .setup file exists
 !
 setupfile = trim(fileprefix)//'.setup'
 inquire(file=setupfile,exist=setexists)
 if (setexists) then
    call read_setupfile(setupfile,gamma,polyk,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(setupfile,gamma,polyk)
       stop 'please rerun phantomsetup with revised .setup file'
    endif
    !--Prompt to get inputs and write to file
 elseif (id==master) then
    print "(a,/)",trim(setupfile)//' not found: using interactive setup'
    call setup_interactive(polyk,gamma,iexist,id,master,ierr)
    call write_setupfile(setupfile,gamma,polyk)
    stop 'please check and edit .setup file and rerun phantomsetup'
 endif

 !
 ! initialise the equation of state
 !
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('setup','could not initialise equation of state')
 !

 call finish_eos(ieos,ierr)


end subroutine setpart








!-----------------------------------------------------------------------
!+
!  Ask questions of the user to determine which setup to use
!+
!-----------------------------------------------------------------------
subroutine setup_interactive(polyk,gamma,iexist,id,master,ierr)
 use prompting,     only:prompt
 use units,         only:select_unit
 use setstream,       only:set_stream_interactive
 use setunits,      only:set_units_interactive
 real, intent(out)    :: polyk,gamma
 logical, intent(in)  :: iexist
 integer, intent(in)  :: id,master
 integer, intent(out) :: ierr

 ierr = 0

 ! units
 !call set_units_interactive(gr)

 ! stream
 call set_stream_interactive(id,master,stream,ieos,polyk)
 !call prompt('Enter stream radius',stream%rstream,0.)
 !call prompt('Enter stream length',stream%zstream,0.)
 !call prompt('Enter stream mass',stream%mstream,0.)

 ! equation of state
 call prompt('Enter the desired EoS (1=isothermal,2=adiabatic,12=idealplusrad)',ieos)
 select case(ieos)
 !case(15) ! Helmholtz
!    call prompt('Enter temperature',stream%initialtemp,1.0e3,1.0e11)
 case(2)
    call prompt('Enter gamma (adiabatic index)',gamma,1.,7.)
 end select

 !if (need_polyk(stream%iprofile)) then
!    call prompt('Enter polytropic constant (cs^2 if isothermal)',polyk,0.)
 !endif

end subroutine setup_interactive

!-----------------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename,gamma,polyk)
 use infile_utils,  only:write_inopt
 use dim,           only:tagline
 use relaxstar,     only:write_options_relax
 use eos,           only:gmw
 !use setstream,       only:write_options_stream
 use setunits,      only:write_options_units
 real,             intent(in) :: gamma,polyk
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(a)") '# input file for Phantom stream setup'

 call write_options_units(iunit,gr)
 !call write_options_stream(stream,iunit)
 call write_inopt(stream%np,'np','Number of particles',iunit)
 call write_inopt(stream%gamma,'gamma','Adiabatic index',iunit)
 call write_inopt(stream%rstream1,'rstream1','Radius of the 1st stream',iunit)
 call write_inopt(stream%zstream1,'zstream1','Length of the 1st stream along z-axis',iunit)
 call write_inopt(stream%rstream2,'rstream2','Radius of the 2nd stream',iunit)
 call write_inopt(stream%zstream2,'zstream2','Length of the 2nd stream along z-axis',iunit)
 call write_inopt(stream%inclination,'inclination','inclination between streams (degrees)',iunit)
 call write_inopt(stream%yshift,'yshift','"y" coordinate of the shifted stream center (code units)',iunit)
 call write_inopt(stream%vinj1,'vinj1','injection speed (code units)',iunit)
 call write_inopt(stream%vinj2,'vinj2','injection speed (code units)',iunit)
 call write_inopt(stream%offset,'offset','offset between streams (units of rstream)',iunit)
 call write_inopt(stream%dtinject,'dtinject','delta t for particle injection (code units)',iunit)
 call write_inopt(stream%inputfile1,'inputfile1','input file name of the 1st cylinder',iunit)
 call write_inopt(stream%inputfile2,'inputfile2','input file name of the 2nd cylinder',iunit)

 call write_inopt(ieos,'ieos','1=isothermal,2=adiabatic,12=idealplusrad',iunit)

 !write(iunit,"(/,a)") '# equation of state'




 !if (need_polyk(stream%iprofile)) call write_inopt(polyk,'polyk','polytropic constant (cs^2 if isothermal)',iunit)

 close(iunit)

end subroutine write_setupfile
!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,gamma,polyk,ierr)
 use infile_utils,  only:open_db_from_file,inopts,read_inopt,close_db
 use io,            only:error
 use units,         only:select_unit
 use relaxstar,     only:read_options_relax
 use eos,           only:X_in,Z_in,gmw
 use eos_gasradrec, only:irecomb
 !use setstream,       only:read_options_stream
 use setunits,      only:read_options_and_set_units
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma,polyk
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return

 nerr = 0

 ! units
 call read_options_and_set_units(db,nerr,gr)

 ! stream options
 !call read_options_stream(stream,need_iso,ieos,polyk,db,nerr)
 call read_inopt(ieos,'ieos',db,errcount=nerr)
 call read_inopt(stream%np,'np',db,errcount=nerr)
 call read_inopt(stream%gamma,'gamma',db,errcount=nerr)
 call read_inopt(stream%rstream1,'rstream1',db,errcount=nerr)
 call read_inopt(stream%zstream1,'zstream1',db,errcount=nerr)
 call read_inopt(stream%rstream2,'rstream2',db,errcount=nerr)
 call read_inopt(stream%zstream2,'zstream2',db,errcount=nerr)
 call read_inopt(stream%inclination,'inclination',db,errcount=nerr)
 call read_inopt(stream%yshift,'yshift',db,errcount=nerr)
 call read_inopt(stream%vinj1,'vinj1',db,errcount=nerr)
 call read_inopt(stream%vinj2,'vinj2',db,errcount=nerr)
 call read_inopt(stream%offset,'offset',db,errcount=nerr)
 call read_inopt(stream%inputfile1,'inputfile1',db,errcount=nerr)
 call read_inopt(stream%inputfile2,'inputfile2',db,errcount=nerr)

 ! equation of state
 !select case(ieos)
 !case(15) ! Helmholtz
  !  call read_inopt(stream%initialtemp,'initialtemp',db,errcount=nerr)
 !end select

 ! option to write density profile to file
 !call read_inopt(write_rho_to_file,'write_rho_to_file',db)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_stream: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile


end module setup

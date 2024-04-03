!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup procedure for a cylinder with a uniform or Gaussian transverse density profile.
! The cylinder can be shifted from the origin and/or rotated.
! The setup can be run with RADIATION=yes. In this case, temperature has to be provided
! and physical units have to be used.
! At the moment, adiabatic (2) and idealplusrad (12) EOS can be used
!
! :References: None
!
! :Owner: Taj Jankovič
!
! :Runtime parameters:
!   - gamma             : *Adiabatic index*
!   - ieos              : *2=adiabatic,12=idealplusrad*
!   - shift             : *shift and/or rotate cylinder automatically during setup*
!   - relax_cyl         : *relax cylinder automatically during setup*
!   - write_rho_to_file : *write density profile(s) to file*
!
!
 use io,             only:fatal,error,warning,master
 use part,           only:gravity,gr
 use physcon,        only:solarm,solarr,km,pi,c,kb_on_mh,radconst
 use options,        only:nfulldump,iexternalforce
 use timestep,       only:tmax,dtmax
 use eos,                only:ieos
 use externalforces,     only:iext_densprofile
 use setcylinder,            only:cylinder_t
 use setunits,           only:dist_unit,mass_unit
 implicit none
 !
 ! Input parameters
 !
 logical            :: iexist
 logical            :: relax_cylinder_in_setup,write_rho_to_file
 type(cylinder_t)       :: cylinder

 public             :: setpart
 private

contains

!-----------------------------------------------------------------------
!+
!  Setup routine for cylinder calculations
!+
!-----------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use units,           only:set_units,select_unit
 use kernel,          only:hfact_default
 use eos,             only:init_eos,finish_eos,gmw
 use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,eos_vars,rad
 use mpiutils,        only:reduceall_mpi
 use mpidomain,       only:i_belong
 use setup_params,    only:rhozero,npart_total
 use setcylinder,         only:set_cylinder
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
 logical                          :: setexists
 character(len=120)               :: setupfile,inname
 !
 ! Initialise parameters, including those that will not be included in *.setup
 !
 time         = 0.
 polyk        = 1.0
 gamma        = 4./3.
 gmw         = 0.5988
 hfact        = hfact_default
 relax_cylinder_in_setup = .false.
 write_rho_to_file = .false.

 !
 ! set default options
 !
 dist_unit   = 'solarr' !maybe change here?
 mass_unit   = 'solarm'

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
 ! determine if the .setup file exists
 !
 setupfile = trim(fileprefix)//'.setup'
 inquire(file=setupfile,exist=setexists)
 if (setexists) then
    call read_setupfile(setupfile,gamma,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(setupfile,gamma)
       stop 'please rerun phantomsetup with revised .setup file'
    endif
    !--Prompt to get inputs and write to file
 elseif (id==master) then
    print "(a,/)",trim(setupfile)//' not found: using interactive setup'
    call setup_interactive(gamma,id,master)
    call write_setupfile(setupfile,gamma)
    stop 'please check and edit .setup file and rerun phantomsetup'
 endif

 !
 ! Verify correct pre-processor commands
 !
 if (.not.gravity) then
    iexternalforce = iext_densprofile
    write_rho_to_file = .true.
 endif

 !
 ! initialise the equation of state
 !
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('setup','could not initialise equation of state')

 !
 ! set up particles
 !
 npartoftype(:) = 0
 npart          = 0
 nptmass        = 0
 vxyzu          = 0.0
 call set_cylinder(id,master,cylinder,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
               massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
               relax_cylinder_in_setup,write_rho_to_file,&
               rhozero,npart_total,i_belong,ierr)
 !
 ! finish/deallocate equation of state tables
 !
 call finish_eos(ieos,ierr)

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Ask questions of the user to determine which setup to use
!+
!-----------------------------------------------------------------------
subroutine setup_interactive(gamma,id,master)
 use prompting,     only:prompt
 use units,         only:select_unit
 use setcylinder,       only:set_cylinder_interactive
 use setunits,      only:set_units_interactive,dist_unit,mass_unit
 real, intent(out)    :: gamma
 integer, intent(in)  :: id,master

 ! units
 call prompt('Enter the desired units (0="code" or 1="solar")',cylinder%iunits)
 if (cylinder%iunits ==1 ) then
   dist_unit   = 'solarr'
   mass_unit   = 'solarm'
   call set_units_interactive()
   write(*,'(a)') 'Using solar units'
 else
   write(*,'(a)') 'Using code units'
 endif

 ! equation of state
 call prompt('Enter the desired EoS (2=adiabatic,12=idealplusrad)',ieos)
 select case(ieos)
 case(2)
    call prompt('Enter gamma (adiabatic index)',gamma,1.,7.)
 end select

 ! cylinder
 call set_cylinder_interactive(id,master,cylinder,ieos)

 ! relaxation
 call prompt('Relax cylinder automatically during setup?',relax_cylinder_in_setup)

end subroutine setup_interactive

!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,gamma,ierr)
 use infile_utils,  only:open_db_from_file,inopts,read_inopt,close_db
 use relaxstar,     only:read_options_relax
 use units,         only:select_unit
 use setcylinder,       only:read_options_cylinder
 use setunits,      only:read_options_and_set_units
 use physcon,       only:solarr,solarm
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma

 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return

 nerr = 0

 ! set units or keep code units
 call read_inopt(cylinder%iunits,'iunits',db,errcount=nerr)
 if (cylinder%iunits==1)   call read_options_and_set_units(db,nerr)

 ! equation of state
 call read_inopt(ieos,'ieos',db,errcount=nerr)
 call read_inopt(gamma,'gamma',db,errcount=nerr)

 ! cylinder options
 call read_inopt(cylinder%np,'np',db,errcount=nerr)

 ! cylinder options
 call read_options_cylinder(cylinder,ieos,db,nerr)

 ! relax cylinder options
 call read_inopt(relax_cylinder_in_setup,'relax_cylinder',db,errcount=nerr)
 if (relax_cylinder_in_setup) call read_options_relax(db,nerr)
 if (nerr /= 0) ierr = ierr + 1

 ! option to write density profile to file
 call read_inopt(write_rho_to_file,'write_rho_to_file',db)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_cylinder: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile




!-----------------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename,gamma)
 use infile_utils,  only:write_inopt
 use dim,           only:tagline
 use relaxstar,     only:write_options_relax
 use setunits,      only:write_options_units
 use setcylinder,   only:write_options_cylinder
 !use physcon,       only:solarr,solarm
 real,             intent(in) :: gamma
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(a)") '# input file for Phantom cylinder setup'

 ! units
 write(iunit,"(/,a)") '# options for units '
 call write_inopt(cylinder%iunits,'iunits','Code units (0) or solar units (1)',iunit)
 if (cylinder%iunits==1) call write_options_units(iunit)

 ! EOS
 write(iunit,"(/,a)") '# equation of state'
 call write_inopt(ieos,'ieos','2=adiabatic,12=idealplusrad',iunit)
 call write_inopt(gamma,'gamma','Adiabatic index',iunit)

 ! options for cylinder
 call write_options_cylinder(cylinder,iunit,ieos)

 ! relaxation options
 write(iunit,"(/,a)") '# relaxation options'
 call write_inopt(relax_cylinder_in_setup,'relax_cylinder','Relax cylinder automatically during setup',iunit)
 if (relax_cylinder_in_setup) call write_options_relax(iunit)

 call write_inopt(write_rho_to_file,'write_rho_to_file','Write density profile(s) to file',iunit)

 close(iunit)

end subroutine write_setupfile



end module setup

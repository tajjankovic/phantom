!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setcylinder
!
! General routine for setting up a cylinder from a 1D profile
! This is the main functionality from setup_cylinder but in a single routine
! that can also be called from other setups.
!
! The routine creates a cylinder with a desired density profile (obtained from a function)
! and optionally relaxes the particle distribution
! :References: None
!
! :Owner: Taj Jankovič
!
! :Runtime parameters: None
!
!
 use unifdis,    only:mask_prototype,mask_true
 use physcon,    only:pi
 use stretchmap, only:rho_func
 use dim,                only:do_radiation

 implicit none

 !
 ! define a data type with all the options needed
 ! to setup cylinder (these are per-cylinder, not per-simulation options)
 !
 type cylinder_t
    integer :: iunits
    integer :: iprofile
    integer :: np
    logical :: ishift
    real :: mcylinder
    real :: rcylinder
    real :: zcylinder
    real :: gauss_sigma
    real :: rotate_theta
    real :: rotate_phi
    real :: xshift
    real :: yshift
    real :: zshift
    real :: initialtemp
    character(len=120) :: dens_profile
    character(len=2) :: label ! used to rename relax_star snapshots to relax1, relax2 etc.
 end type cylinder_t

 public :: cylinder_t
 public :: set_cylinder,set_defaults_cylinder,shift_cylinder
 public :: read_options_cylinder,write_options_cylinder,set_cylinder_interactive

 real,    private :: tol_ekin = 1.e-7 ! criteria for being converged
 integer, private :: maxits = 1000
 real,    private :: gammaprev,hfactprev
 integer, private :: ieos_prev
 integer, public :: ierr_setup_errors = 1, &
                    ierr_notconverged = 2
 private

contains

!--------------------------------------------------------------------------
!+
!  default options for a particular cylinder
!  (see also set_defaults_given_profile which selects defaults
!   based on the value of iprofile)
!+
!--------------------------------------------------------------------------
subroutine set_defaults_cylinder(cylinder)
 type(cylinder_t), intent(out) :: cylinder
 cylinder%iunits        = 0
 cylinder%np            = 1000
 cylinder%iprofile    = 1
 cylinder%gauss_sigma   = 0.
 cylinder%mcylinder       = 1.0
 cylinder%rcylinder       = 1.0
 cylinder%zcylinder       = 2.0
 cylinder%ishift    = .true.
 cylinder%rotate_theta= 90.
 cylinder%rotate_phi = 0.
 cylinder%xshift        = 0.
 cylinder%yshift        = 4.
 cylinder%zshift        = 0.
 cylinder%initialtemp   = 0. ! temperature is assigned only if ieos=12 or do_radiation
 cylinder%dens_profile         = 'density_out.tab'
 cylinder%label         = ''
end subroutine set_defaults_cylinder



!--------------------------------------------------------------------------
!+
!  default options for a particular cylinder in solar units
! default values are obtained assuming R=10Rsun, v=0.01c, dotM=0.7Msun/yr (see Jankovič et al. 2023)
!+
!--------------------------------------------------------------------------
subroutine set_defaults_cylinder_solar(cylinder)
 use units,   only:udist,umass
 use physcon, only:solarm,solarr
 type(cylinder_t), intent(out) :: cylinder

 cylinder%iunits        = 1
 cylinder%np            = 1000
 cylinder%iprofile    = 1
 cylinder%gauss_sigma   = 0.
 cylinder%mcylinder       = 1.035e-4*real(solarm/umass)
 cylinder%rcylinder       = 10.0*real(solarr/udist)
 cylinder%zcylinder       = 20.0*real(solarr/udist)
 cylinder%ishift    = .true.
 cylinder%rotate_theta= 90.
 cylinder%rotate_phi = 0.
 cylinder%xshift        = 0.
 cylinder%yshift        = 40.0*real(solarr/udist)
 cylinder%zshift        = 0.
 cylinder%initialtemp   = 1.0e4 ! setup is appropriate for a cold cylinder
 cylinder%dens_profile         = 'density_out.tab'
 cylinder%label         = ''
end subroutine set_defaults_cylinder_solar

!--------------------------------------------------------------------------
!+
!  Master routine to setup a cylinder from a specified file or density profile
!+
!--------------------------------------------------------------------------
subroutine set_cylinder(id,master,cylinder,xyzh,vxyzu,eos_vars,rad,&
                    npart,npartoftype,massoftype,hfact,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
                    relax,write_rho_to_file,&
                    rhomean,npart_total,mask,ierr,itype)
 use centreofmass,       only:reset_centreofmass
 use dim,                only:do_radiation,gravity,maxvxyzu
 use io,                 only:fatal,error,warning
 use setstar_utils,      only:set_star_thermalenergy
 use radiation_utils,    only:set_radiation_and_gas_temperature_equal
 use part,               only:igas,set_particle_type,ilum
 use extern_densprofile, only:write_rhotab
 use unifdis,            only:mask_prototype
 use physcon,            only:pi
 use units,              only:utime,udist,umass,unit_density
 use mpiutils,           only:reduceall_mpi
 type(cylinder_t), intent(inout)  :: cylinder
 integer,      intent(in)     :: id,master
 integer,      intent(inout)  :: npart,npartoftype(:),nptmass
 real,         intent(inout)  :: xyzh(:,:),vxyzu(:,:),eos_vars(:,:),rad(:,:)
 real,         intent(inout)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,         intent(inout)  :: massoftype(:)
 real,         intent(in)     :: hfact
 logical,      intent(in)     :: relax,write_rho_to_file
 integer,      intent(in)     :: ieos
 real,         intent(inout)  :: polyk,gamma
 real,         intent(out)    :: rhomean
 integer(kind=8), intent(out) :: npart_total
 integer,      intent(out)    :: ierr
 integer,      intent(in), optional :: itype
 procedure(mask_prototype)    :: mask
 integer                        :: npts,ierr_relax
 integer                        :: npart_old,i
 real, allocatable              :: r(:),den(:),pres(:),temp(:)
 real                           :: rhocentre
 character(len=30)              :: lattice  ! The lattice type='random' is used

 ierr_relax = 0
 rhomean = 0.
 npart_old = npart

 !
 ! do nothing if iprofile is invalid or zero (sink particle)
 !
 if (cylinder%iprofile <= 0) then
    ierr = 1
    return
 endif

 !
 ! get the desired tables of density as a function of transverse radius
 !
 call set_cylinder_profile(cylinder%iprofile,ieos,r,den,npts,&
                        cylinder%rcylinder,cylinder%zcylinder,cylinder%mcylinder,cylinder%gauss_sigma,rhocentre)

 !
 ! set up particles to represent the desired stellar profile
 !
 lattice = 'random'
 if (cylinder%np < 1 .and. npart_old==0) then
    call fatal('set_cylinder','cannot set up a cylinder with zero particles')
    ierr = 2
    return
 endif
 if (cylinder%mcylinder < 0.) then
    call fatal('set_cylinder','cannot set up a cylinder with negative mass!')
    ierr = 2
    return
 endif
 call set_cylinder_density(id,master,cylinder%rcylinder,cylinder%zcylinder,cylinder%mcylinder,&
                       hfact,npts,den,npart,npartoftype,massoftype,xyzh,&
                       cylinder%np,rhomean,npart_total,mask)

 !
 ! Write the desired profile to file (do this before relaxation)
 !
 if (write_rho_to_file) call write_rhotab(cylinder%dens_profile,&
                                          r,den,npts,polyk,gamma,rhocentre,ierr)

 !
 ! relax the density profile to achieve nice hydrostatic equilibrium
 !
 if (relax) then
    if (reduceall_mpi('+',npart)==npart) then
       call relax_cylinder(npts,den,pres,r,npart,xyzh,cylinder%rcylinder,cylinder%zcylinder,&
                       ierr_relax,npin=npart_old,label=cylinder%label)
    else
       call error('setup_cylinder','cannot run relaxation with MPI setup, please run setup on ONE MPI thread')
    endif
 endif

 !
 ! Reset centre of mass (again)
 !
 call reset_centreofmass(npart,xyzh,vxyzu)

 !
 ! set the internal energy and temperature
 !
 !if (maxvxyzu==4) then
  !  call set_star_thermalenergy(ieos,den,pres,r,npts,npart,&
  !                     xyzh,vxyzu,rad,eos_vars,relax,.false.,cylinder%initialtemp,&
  !                     npin=npart_old)
 !endif

 if (do_radiation) then
     call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad,&
                                                    npin=npart_old)
 endif

 !
 ! shift cylinder to requested position and velocity
 !
 if (cylinder%ishift) call shift_cylinder(npart,xyzh,cylinder%rotate_theta,cylinder%rotate_phi,cylinder%xshift,cylinder%yshift,&
                                      cylinder%zshift)

 !
 ! give the particles requested particle type:
 ! this is initially used to tag particles in different cylinder
 ! so the cylinder can later be shifted into position
 !
 if (present(itype)) then
    do i=npart_old+1,npart
       call set_particle_type(i,itype)
    enddo
 endif
 !
 ! Print summary to screen
 !
 if (id==master) then
    write(*,"(70('='))")
    write(*,'(1x,a,f12.5)')       'gamma                                  = ', gamma
    write(*,'(1x,a,i12)')         'Number of particles                    = ', npart_total-npart_old
    if (cylinder%iunits==0) then
      write(*,'(1x,a,2(es12.5,a))') 'Particle mass            = ',massoftype(igas),   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Cylinder mass            = ',cylinder%mcylinder,   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Cylinder radius          = ',cylinder%rcylinder,   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Cylinder height          = ',cylinder%zcylinder,   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Central density          = ',rhocentre,   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Average density          = ',rhomean,   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Cylinder center ("x")                        = ',cylinder%xshift,   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Cylinder center ("y")                        = ',cylinder%yshift,   ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Cylinder center ("z")                        = ',cylinder%zshift,   ' code units'
    else
      call write_mass('Particle mass                          = ',massoftype(igas),umass)
      call write_mass('Cylinder mass                          = ',cylinder%mcylinder,udist)
      call write_dist('Cylinder radius                        = ',cylinder%rcylinder,udist)
      call write_dist('Cylinder height                        = ',cylinder%zcylinder,udist)
      write(*,'(1x,a,2(es12.5,a))') 'Central density                        = ', rhocentre*unit_density,' g/cm^3 ='&
                                                                 , rhocentre,               ' code units'
      write(*,'(1x,a,2(es12.5,a))') 'Average density                        = ', rhomean*unit_density,  ' g/cm^3 = '&
                                                                 , rhomean,               ' code units'
      call write_dist('Cylinder center ("x")                        = ',cylinder%xshift,udist)
      call write_dist('Cylinder center ("y")                        = ',cylinder%yshift,udist)
      call write_dist('Cylinder center ("z")                        = ',cylinder%zshift,udist)
    endif
endif

write(*,'(1x,a,2(es12.5,a))') 'Cylinder rotation angle around "x" axis      = ',cylinder%rotate_theta,   ' degrees'
write(*,'(1x,a,2(es12.5,a))') 'Cylinder rotation angle around "z" axis      = ',cylinder%rotate_phi,   ' degrees'
write(*,"(70('='))")
if (ierr_relax /= 0) call warning('setup_cylinder','ERRORS DURING RELAXATION, SEE ABOVE!!')

end subroutine set_cylinder

!-------------------------------------------------------------------------------
!+
!  set up particles according to the desired density profile
!+
!-------------------------------------------------------------------------------
subroutine set_cylinder_density(id,master,Rcylinder,Zcylinder,Mcylinder,hfact,&
                            npts,den,npart,npartoftype,massoftype,xyzh,&
                            np,rhomean,npart_total,mask)
 use part,         only:set_particle_type,igas
 use unifdis,      only:mask_prototype
 use physcon,      only:pi
 integer,          intent(in)    :: id,master,npts,np
 integer(kind=8),  intent(out)   :: npart_total
 real,             intent(in)    :: Rcylinder,Zcylinder,Mcylinder,hfact
 real,             intent(in)    :: den(npts)
 real,             intent(out)   :: rhomean
 integer,          intent(inout) :: npart,npartoftype(:)
 real,             intent(inout) :: xyzh(:,:),massoftype(:)
 procedure(mask_prototype) :: mask
 integer :: i,ntot,npart_old,n
 real :: vol_cylinder
 logical :: mass_is_set

 ! if this is the first call to set_cylinder_density (npart_old=0),
 ! set the number of particles to the requested np and set
 ! the particle mass accordingly, otherwise assume particle mass
 ! is already set and use Mstar to determine np
 !
 npart_old = npart
 n = np
 mass_is_set = .false.
 if (massoftype(igas) > tiny(0.)) then
    n = nint(Mcylinder/massoftype(igas))
    mass_is_set = .true.
    print "(a,i0)",' WARNING: particle mass is already set, using np = ',n
 endif
 !
 ! place particles in a cylinder
 !
 vol_cylinder  = pi*Rcylinder**2.*Zcylinder
 call set_cylinder_particles(id,master,npts,den,Rcylinder,Zcylinder,hfact,npart,xyzh,npart_total,n,mask)
 !
 ! get total particle number across all MPI threads
 ! and use this to set the particle mass
 !
 ntot = int(npart_total,kind=(kind(ntot)))
 if (.not.mass_is_set) then
    massoftype(igas) = Mcylinder/ntot
 endif
 !
 ! set particle type as gas particles
 !
 npartoftype(igas) = npart   ! npart is number on this thread only
 do i=npart_old+1,npart_old+npart
    call set_particle_type(i,igas)
 enddo
 !
 ! mean density
 !
 rhomean = Mcylinder/vol_cylinder

end subroutine set_cylinder_density



!-----------------------------------------------------------------------
!+
!  This subroutine positions particles on a cylinder - uniform
!  density by default, or if a tabulated array rhotab(:) is passed, with
!  an arbitrary radial density profile rho(r)
!
!  The table rhotab(:) is assumed to contain density values
!  on equally spaced radial bins between rmin and rmax
!+
!-----------------------------------------------------------------------
subroutine set_cylinder_particles(id,master,npts,rhotab,rcylinder,zcylinder,&
                      hfact,np,xyzh,nptot,np_requested,mask)
 use stretchmap, only:set_density_profile,rho_func
 integer,          intent(in)    :: id,master,npts
 integer,          intent(inout) :: np
 real,             intent(in)    :: rhotab(npts),rcylinder,zcylinder,hfact
 real,             intent(out)   :: xyzh(:,:)
 integer,          intent(in),    optional :: np_requested
 integer(kind=8),  intent(inout), optional :: nptot
 procedure(mask_prototype), optional :: mask
 procedure(mask_prototype), pointer  :: my_mask
 integer                         :: ierr

 if (present(mask)) then
    my_mask => mask
 else
    my_mask => mask_true
 endif
 !
 !--Create a cylinder of uniform density
 !
 call set_cylinder_mc(id,master,rcylinder,zcylinder,hfact,np_requested,np,xyzh,ierr,nptot,my_mask)

 !
 !--Stretching the spatial distribution to have the desired density distribution
 !
 call set_density_profile(np,xyzh,rhotab=rhotab,min=0.,max=rcylinder,geom=2)
contains
end subroutine set_cylinder_particles



!-----------------------------------------------------------------------
!+
! set up uniform spherical particle distribution using Monte Carlo particle
! placement. Particles are placed in pairs so the distribution
! is symmetric about the origin
!+
!-----------------------------------------------------------------------
subroutine set_cylinder_mc(id,master,rmax,zmax,hfact,np_requested,np,xyzh, &
                         ierr,nptot,mask)
 use random,     only:ran2
 integer,          intent(in)    :: id,master,np_requested
 integer,          intent(inout) :: np   ! number of actual particles
 real,             intent(in)    :: rmax,zmax,hfact
 real,             intent(out)   :: xyzh(:,:)
 integer,          intent(out)   :: ierr
 integer(kind=8),  intent(inout), optional :: nptot
 procedure(mask_prototype) :: mask
 integer :: i,iseed,maxp
 real    :: vol_cylinder,rr,phi,z,psep,dir(2)
 integer(kind=8) :: iparttot

 iparttot = 0
 vol_cylinder = pi*rmax**2.*zmax
 ! use mean particle spacing to set initial smoothing lengths
 psep = (vol_cylinder/real(np_requested))**(1./3.)
 iseed = -1978
 maxp  = size(xyzh(1,:))
 ierr  = 1

 do i=1,np_requested,2
    !
    ! get a random position on a cylinder
    !
    rr = ran2(iseed)*rmax ! uniform density
    phi = 2. * pi * ran2(iseed)
    z = zmax * (ran2(iseed) - 0.5) ! Generate a random z position within [-zmax/2, zmax/2]
    dir  = (/cos(phi),sin(phi)/)
    !
    ! add TWO particles, symmetric around the origin
    ! CHECK HERE - MAYBE NOT NECESSARY
    !
    iparttot = iparttot + 1
    if (mask(iparttot)) then
       np = np + 1
       if (np > maxp) then
          print*,' ERROR: np > array size: use ./phantomsetup --maxp=',np_requested
          return
       endif
       xyzh(1:2,np) = rr*dir
       xyzh(3,np) = z
       xyzh(4,np)   = hfact*psep
    endif
    iparttot = iparttot + 1
    if (mask(iparttot)) then
       np = np + 1
       if (np > maxp) then
          print*,' ERROR: np > array size: use ./phantomsetup --maxp=',np_requested
          return
       endif
       xyzh(1:2,np) = -rr*dir
       xyzh(3,np) = z
       xyzh(4,np)   = hfact*psep
    endif
 enddo
 ierr = 0
 if (present(nptot)) nptot = iparttot
 if (id==master) write(*,"(1x,a,i10,a)") 'placed ',np,' particles in random-but-symmetric sphere'

end subroutine set_cylinder_mc



!-----------------------------------------------------------------------
!+
!  Shift cylinder origin from (0,0,0) to the desired position.
!  Rotate cylinder for about z-axis for an angle inc (inc=0 -> cylinder axis aligned
!  with z-axis).
!+
!-----------------------------------------------------------------------
subroutine shift_cylinder(npart,xyzh,rot_theta,rot_phi,x0,y0,z0)
 use part,        only:get_particle_type,set_particle_type,igas
 use vectorutils, only:rotatevec
 integer, intent(in) :: npart
 real, intent(inout) :: xyzh(:,:)
 real, intent(in)    :: x0,y0,z0,rot_theta,rot_phi
 integer             :: i

 do i=1,npart
    !rotate about x-axis
    call rotatevec(xyzh(1:3,i),(/1.,0.,0./),rot_theta/180.*pi)

    !rotate about z-axis
    call rotatevec(xyzh(1:3,i),(/0.,0.,1./),rot_phi/180.*pi)

    !shift cylinder
    xyzh(1,i) = xyzh(1,i) + x0
    xyzh(2,i) = xyzh(2,i) + y0
    xyzh(3,i) = xyzh(3,i) + z0
 enddo
end subroutine shift_cylinder



!-------------------------------------------------------------------------------
!+
!  set cylinder profile from a function
!+
!-------------------------------------------------------------------------------
subroutine set_cylinder_profile(iprofile,ieos,r,den,npts,Rcylinder,Zcylinder,Mcylinder,&
                             gauss_sigma,rhocentre)
 use extern_densprofile, only      : nrhotab
 integer,           intent(in)    :: iprofile,ieos
 real, allocatable, intent(out)   :: r(:),den(:)
 integer,           intent(out)   :: npts
 real,              intent(inout) :: Rcylinder,Zcylinder,Mcylinder,rhocentre,gauss_sigma
 real                             :: dr
 integer, parameter               :: ng_max = nrhotab
 integer                          :: i

 !
 ! set up tabulated density profile
 !
 allocate(r(ng_max),den(ng_max))

 dr      = Rcylinder/real(ng_max)
 do i=1,ng_max
    r(i)   = i*dr
    den(i)   = densfunc(r(i))
 enddo
 npts = ng_max
 contains
 !-------------------------------------------------
 !+
 !  callback function for desired density profile
 !+
 !-------------------------------------------------
  real function densfunc(r)
   real, intent(in) :: r

   if (iprofile==1) then
     rhocentre = 1/3.14 !from dotM=2*pi*v*rho_0*int_0^H r*exp(-r^2/H^2)dr
     densfunc = rhocentre
   elseif (iprofile==2) then
     rhocentre = mcylinder/(zcylinder*0.5*gauss_sigma**2.*(1-EXP(-rcylinder**2./gauss_sigma**2.)))
     !from mcylinder=2*pi*rho_0*int_0^Rcylinder r*exp(-r^2/H^2)dr
     !https://www.wolframalpha.com/input?i=integrate+%5Bx*exp%28-x%5E2%2F%28S%5E2%29%29%2C+%7Bx%2C0%2CH%7D%5D
     densfunc = rhocentre*EXP(-r**2./gauss_sigma**2.)
   endif
 end function densfunc

end subroutine set_cylinder_profile




!----------------------------------------------------------------
!+
!  relax a cylinder to hydrostatic equilibrium. We run the main
!  code but with a fake equation of state, low neighbour number
!  and fixing the entropy as a function of r
!
!  IN:
!    rho(nt)   - tabulated density as function of r (in code units)
!    r(nt)     - radius for each point in the table
!
!  IN/OUT:
!    xyzh(:,:) - positions and smoothing lengths of all particles
!+
!----------------------------------------------------------------
subroutine relax_cylinder(nt,rho,pr,r,npart,xyzh,Rcylinder,Zcylinder,ierr,npin,label)
 use table_utils,     only:yinterp
 use deriv,           only:get_derivs_global
 use dim,             only:maxp,maxvxyzu,gravity
 use part,            only:vxyzu,rad,massoftype,igas
 use step_lf_global,  only:init_step,step
 use initial,         only:initialise
 use memory,          only:allocate_memory
 use energies,        only:compute_energies,ekin,epot,etherm
 use checksetup,      only:check_setup
 use io,              only:error,warning,fatal,id,master
 use fileutils,       only:getnextfilename
 use readwrite_dumps, only:write_fulldump,init_readwrite_dumps
 use physcon,         only:pi
 use io_summary,      only:summary_initialise
 use setstar_utils,   only:set_star_thermalenergy
 integer, intent(in)    :: nt
 integer, intent(inout) :: npart
 real,    intent(in)    :: rho(nt),pr(nt),r(nt),Rcylinder,Zcylinder
 real,    intent(inout) :: xyzh(:,:)
 integer, intent(out)   :: ierr
 integer, intent(in), optional :: npin
 character(len=*), intent(in), optional :: label
 integer :: nits,nerr,nwarn,iunit,i1
 real    :: t,dt,rmserr
 real    :: mr(nt),rmax,dtext
 logical :: converged,restart
 logical, parameter :: write_files = .true.
 character(len=20) :: filename,mylabel

 i1 = 0
 if (present(npin)) i1 = npin  ! starting position in particle array
 !
 ! label for relax_star snapshots
 !
 mylabel = ''
 if (present(label)) mylabel = label
 !
 ! save settings and set a bunch of options
 !
 ierr = 0
 mr = get_mr(rho,r)
 !
 ! see if we can restart or skip the relaxation process
 ! based on a previous run
 !
 filename = 'relax'//trim(mylabel)//'_00000'
 if (write_files) call init_readwrite_dumps()
 call check_for_existing_file(filename,npart,massoftype(igas),&
                              xyzh,vxyzu,restart,ierr)
 !
 ! quit with fatal error if non-matching file found, otherwise
 ! this will be overwritten
 !
 if (write_files .and. ierr /= 0) then
    if (id==master) print "(a)",' ERROR: pre-existing relaxed star dump(s) do not match'
    call fatal('relax_star','please delete relax'//trim(mylabel)//'_* and restart...')
 endif

 call set_options_for_relaxation()
 call summary_initialise()
 !
 ! check particle setup is sensible
 !
 call check_setup(nerr,nwarn,restart=.true.) ! restart=T allows accreted/masked particles
 if (nerr > 0) then
    call error('relax_star','cannot relax star because particle setup contains errors')
    call restore_original_options(i1,npart)
    ierr = ierr_setup_errors
    return
 endif

 call set_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,r,rmax,rmserr)
 !
 ! compute derivatives the first time around (needed if using actual step routine)
 !
 t = 0.
 call allocate_memory(int(min(2*npart,maxp),kind=8))
 call get_derivs_global()
 call set_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,r,rmax,rmserr)
 call compute_energies(t)
 !
 ! perform sanity checks
 !

 if (id==master) print "(/,3(a,1pg11.3),/,a,1pg11.3,a,i4)",&
   ' RELAX-A-STAR-O-MATIC: Etherm:',etherm,' Ekin:',ekin, ' R*:',maxval(r), &
   '       WILL stop when Iter=',maxits

 if (write_files) then
    if (.not.restart) call write_fulldump(t,filename)
    open(newunit=iunit,file='relax'//trim(mylabel)//'.ev',status='replace')
    write(iunit,"(a)") '# nits,rmax,etherm,epot,ekin,L2_{err}'
 endif
 converged = .false.
 dt = epsilon(0.) ! To avoid error in sink-gas substepping
 dtext = huge(dtext)

 nits = 0
 do while (.not. converged .and. nits < maxits)
    nits = nits + 1
    !
    ! shift particles by one "timestep"
    !
    call shift_particles(i1,npart,xyzh,vxyzu,Rcylinder,Zcylinder,dt)

    !
    ! reset thermal energy and calculate information
    !
    call set_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,r,rmax,rmserr)
    !
    ! compute energies and check for convergence
    !
    call compute_energies(t)
    converged = (ekin > 0. .and. ekin/abs(epot) < tol_ekin) !.and. rmserr < 0.01*tol_dens)
    !
    ! print information to screen
    !
    if (id==master .and. mod(nits,10)==0 .or. nits==1) then
      print "(a,i4,a,i4,a,2pf6.2,2(a,1pg11.3))",&
            ' Relaxing star: Iter',nits,'/',maxits, &
            ', dens error:',rmserr,'%, R*:',rmax,' Ekin/Epot:',ekin/abs(epot)
    endif
    !
    ! additional diagnostic output, mainly for debugging/checking
    !
    if (write_files) then
       !
       ! write information to the relax.ev file
       !
       write(iunit,*) nits,rmax,etherm,epot,ekin,rmserr
       !
       ! write dump files
       !
       if (mod(nits,100)==0 .or. ((nits==maxits .or. converged).and.nits > 1)) then
          filename = getnextfilename(filename)
          !
          ! before writing a file, set the real thermal energy profile
          ! so the file is useable as a starting file for the main calculation
          !

          !if (maxvxyzu==4) call set_star_thermalenergy(ieos_prev,rho,pr,&
          !                      r,nt,npart,xyzh,vxyzu,rad,eos_vars,.true.,&
          !                      use_var_comp=.false.,initialtemp=1.e3,npin=i1)

          ! write relaxation snapshots
          if (write_files) call write_fulldump(t,filename)

          ! flush the relax.ev file
          call flush(iunit)

          ! restore the fake thermal energy profile
          call set_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,r,rmax,rmserr)
       endif
    endif
 enddo
 if (write_files) close(iunit)
 !
 ! warn if relaxation finished due to hitting nits=nitsmax
 !
 if (.not.converged) then
    call warning('relax_star','relaxation did not converge, just reached max iterations')
    ierr = ierr_notconverged
 else
    if (id==master) print "(5(a,/))",&
    "                             _                    _ ",&
    "                    _ __ ___| | __ ___  _____  __| |",&
    "                   | '__/ _ \ |/ _` \ \/ / _ \/ _` |",&
    "                   | | |  __/ | (_| |>  <  __/ (_| |",&
    "          o  o  o  |_|  \___|_|\__,_/_/\_\___|\__,_|  o  o  o"
 endif
 !
 ! unfake some things
 !
 call restore_original_options(i1,npart)

end subroutine relax_cylinder



!----------------------------------------------------------------
!+
!  shift particles: this is like timestepping but done
!  asynchronously. Each particle shifts by dx = 0.5*dt^2*a
!  where dt is the local courant timestep, i.e. h/c_s
!+
!----------------------------------------------------------------
subroutine shift_particles(i1, npart, xyzh, vxyzu, r_cylinder, h_cylinder,dtmin)
    use deriv, only:get_derivs_global
    use part,  only:fxyzu,fext,rhoh,massoftype,igas
    use ptmass,only:get_accel_sink_gas
    use eos,   only:get_spsound
    use options, only:ieos
    integer, intent(in) :: i1, npart
    real, intent(inout) :: xyzh(:,:), vxyzu(:,:)
    real, intent(in) :: r_cylinder, h_cylinder
    real, intent(out)   :: dtmin
    real :: dti, rhoi, cs, hi,force_const
    real :: dx(3), distance_from_axis, distance_from_plane
    integer :: i,nlargeshift

    ! shift particles asynchronously
    dtmin = huge(dtmin)
    nlargeshift = 0
    force_const=1e-3
    fext(1:3,:) = 0.0

    do i = i1 + 1, npart
        hi = xyzh(4,i)
        rhoi = rhoh(hi, massoftype(igas))
        cs = get_spsound(ieos, xyzh(:,i), rhoi, vxyzu(:,i))
        dti = 0.3 * hi / cs   ! local Courant timestep, i.e., h/cs

        ! Initialize external force to zero

        ! Calculate distance from the axis of the cylinder
        distance_from_axis = sqrt(xyzh(1, i)**2 + xyzh(2, i)**2)
        ! Calculate distance from the mid-plane of the cylinder
        distance_from_plane = abs(xyzh(3, i)) - h_cylinder / 2.0

        ! Wrap around the z-coordinate if it goes beyond the cylinder half-height
        if (distance_from_plane > 0.0) then
            xyzh(3, i) = xyzh(3, i) - sign(h_cylinder, xyzh(3, i))
        endif

        ! Apply restoring force towards the cylinder's boundary
        if (distance_from_axis>(r_cylinder-hi)) then
        !if (distance_from_axis>(0.9*r_cylinder)) then
          fext(1, i) = -force_const * (xyzh(1, i) / distance_from_axis)
          fext(2, i) = -force_const * (xyzh(2, i) / distance_from_axis)
        endif

        ! Calculate the asynchronous shift
        dx = 0.5 * dti**2 * (fxyzu(1:3, i) + fext(1:3, i))
        if (dot_product(dx,dx) > hi**2) then
           dx = dx / sqrt(dot_product(dx,dx)) * hi  ! Avoid large shift in particle position
           nlargeshift = nlargeshift + 1
        endif
        ! Apply the shift
        xyzh(1:3, i) = xyzh(1:3, i) + dx
        ! Update the velocities based on the shift
        vxyzu(1:3, i) = dx / dti
        dtmin = min(dtmin,dti)   ! used to print a "time" in the output (but it is fake)
    enddo
    if (nlargeshift > 0) print*,'Warning: Restricted dx for ', nlargeshift, 'particles'
    !
    ! get forces on particles
    !
    call get_derivs_global()
end subroutine shift_particles



!----------------------------------------------------------------
!+
!  set code options specific to relaxation calculations
!+
!----------------------------------------------------------------
subroutine set_options_for_relaxation()
 use eos,  only:ieos,gamma
 use part, only:hfact,maxvxyzu
 use damping, only:damp
 use options, only:idamp

 gammaprev = gamma
 hfactprev = hfact
 ieos_prev = ieos
 !
 ! turn on settings appropriate to relaxation
 !
 if (maxvxyzu >= 4) ieos = 2
 idamp = 1
 damp = 0.05
end subroutine set_options_for_relaxation


!----------------------------------------------------------------
!+
!  reset the thermal energy to be exactly p(r)/((gam-1)*rho(r))
!  according to the desired p(r) and rho(r)
!  also compute error between true rho(r) and desired rho(r)
!+
!----------------------------------------------------------------
subroutine set_u_and_get_errors(i1, npart, xyzh, vxyzu, rad, nt, mr, rhotab,rtab,rmax, rmserr)
    use table_utils, only:yinterp
    use sortutils,   only:find_rank,r2func
    use part,        only:rhoh, massoftype, igas, maxvxyzu, iorder=>ll
    use dim,         only:do_radiation
    use eos,         only:gamma
    integer, intent(in) :: i1, npart, nt
    real, intent(in)    :: xyzh(:,:), mr(nt), rhotab(nt), rtab(nt)
    real, intent(inout) :: vxyzu(:,:), rad(:,:)
    real, intent(out)   :: rmax, rmserr
    !logical, intent(in) :: fix_entrop
    real :: ri, rhotarget, rhoi,rho1, Ptarget, Pfactor, mstar
    integer :: i

    ! Define or get the target density and pressure functions
    ! These could be external functions or could be defined elsewhere
    ! real function return_user_desired_target_density(r)
    ! real function return_user_desired_target_pressure(r)
    rho1 = yinterp(rhotab,rtab,0.)


    Pfactor = 5e-6 ! Target pressure at position ri

    rmax = 0.
    rmserr = 0.
    call find_rank(npart-i1,r2func,xyzh(1:3,i1+1:npart),iorder)

    mstar = mr(nt)

    do i = i1+1, npart
        ri = sqrt(dot_product(xyzh(1:2,i), xyzh(1:2,i)))  ! Calculate radial distance from the cylinder axis
        !rhotarget = yinterp(rho,mr,ri)  ! Target density at position ri
        rhotarget = yinterp(rhotab,rtab,ri)  ! Target density at position ri

        rhoi = rhoh(xyzh(4,i), massoftype(igas))  ! Actual density

        ! Calculate target pressure and adjust internal energy accordingly
        Ptarget = Pfactor * (rhoi / rhotarget)
        !if (fix_entrop) then
        !    vxyzu(4,i) = (Pratio * rhoi**(gamma-1.0)) / (gamma - 1.0)
        vxyzu(4,i) = Ptarget / (rhoi * (gamma - 1.0))
        !endif

        ! Update rms error calculation
        rmserr = rmserr + (rhotarget - rhoi)**2
        rmax   = max(rmax, ri)
    end do

    if (do_radiation) rad = 0.
    rmserr = sqrt(rmserr / real(npart)) / rho1

end subroutine set_u_and_get_errors



!--------------------------------------------------
!+
!  get mass coordinate m(r) = \int 2.*pi*rho*r dr
!+
!--------------------------------------------------
function get_mr(rho,r) result(mr)
 use physcon, only:pi
 real, intent(in)  :: rho(:),r(:)
 real :: mr(size(r))
 integer :: i

 mr(1) = 0.
 do i=2,size(rho)
    mr(i) = mr(i-1) + 2.*pi*rho(i) * (r(i) - r(i-1)) * ((r(i) + r(i-1)) / 2.)
 enddo

end function get_mr



!----------------------------------------------------------------
!+
!  restore previous settings
!+
!----------------------------------------------------------------
subroutine restore_original_options(i1,npart)
 use eos,     only:ieos,gamma
 use damping, only:damp
 use options, only:idamp
 use part,    only:hfact,vxyzu
 integer, intent(in) :: i1,npart

 gamma = gammaprev
 hfact = hfactprev
 ieos  = ieos_prev
 idamp = 0
 damp = 0.
 vxyzu(1:3,i1+1:npart) = 0.

end subroutine restore_original_options


!----------------------------------------------------------------
!+
!  check if a previous snapshot exists of the relaxed star
!+
!----------------------------------------------------------------
subroutine check_for_existing_file(filename,npart,mgas,xyzh,vxyzu,restart,ierr)
 use dump_utils, only:open_dumpfile_r,read_header,dump_h,lenid,extract
 use fileutils,  only:getnextfilename
 use io,         only:idump,idisk1,id,nprocs,iprint
 use readwrite_dumps, only:read_dump
 character(len=*), intent(inout) :: filename
 integer,          intent(in)    :: npart
 real,             intent(in)    :: mgas
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical,          intent(out)   :: restart
 integer,          intent(out)   :: ierr
 logical :: iexist
 character(len=len(filename)) :: restart_file,filetmp
 character(len=lenid) :: fileid
 type(dump_h) :: hdr
 integer :: npartfile
 real :: hfactfile,tfile,mfile
 !
 ! check for the last file in the list relax_00000, relax_00001 etc
 !
 ierr = 0
 iexist = .true.
 filetmp = filename
 restart_file = ''
 restart = .false.
 do while (iexist)
    inquire(file=filetmp,exist=iexist)
    if (iexist) then
       restart_file = filetmp
       filetmp = getnextfilename(filetmp)
    endif
 enddo
 if (len_trim(restart_file) <= 0) return

 print "(/,1x,a)",'>> RESTARTING relaxation from '//trim(restart_file)
 call open_dumpfile_r(idump,restart_file,fileid,ierr)
 call read_header(idump,hdr,ierr)
 close(idump)
 if (ierr /= 0) then
    print "(a)",' ERROR: could not read file header'
    return
 else
    call extract('nparttot',npartfile,hdr,ierr,default=0)
    if (npartfile /= npart) then
       print "(a,i0,a,i0)",' ERROR: np=',npartfile,' in '//&
         trim(restart_file)//' differs from current np=',npart
       ierr = 2
       return
    else
       call extract('massoftype',mfile,hdr,ierr,default=0.)
       if (abs(mfile-mgas) > epsilon(0.)) then
          print "(a,es10.3,a,es10.3)",' ERROR: M=',npart*mfile,' in '//&
            trim(restart_file)//' differs from current M=',npart*mgas
          ierr = 3
          return
       else
          restart = .true.
          call read_dump(restart_file,tfile,hfactfile,&
                         idisk1,iprint,id,nprocs,ierr)
          filename = restart_file
       endif
    endif
 endif

end subroutine check_for_existing_file


!-----------------------------------------------------------------------
!+
!  print a distance in both code units and physical units
!+
!-----------------------------------------------------------------------
subroutine write_dist(item_in,dist_in,udist)
 use physcon, only:solarr,km
 real,             intent(in) :: dist_in
 real(kind=8),     intent(in) :: udist
 character(len=*), intent(in) :: item_in

 write(*,'(1x,2(a,es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in*udist/solarr,' R_sun'

end subroutine write_dist
!-----------------------------------------------------------------------
!+
!  print a mass in both code units and physical units
!+
!-----------------------------------------------------------------------
subroutine write_mass(item_in,mass_in,umass)
 use physcon, only:solarm
 real,             intent(in) :: mass_in
 real(kind=8),     intent(in) :: umass
 character(len=*), intent(in) :: item_in
 write(*,'(1x,2(a,es12.5),a)') item_in, mass_in*umass,' g      = ',mass_in*umass/solarm,' M_sun'

end subroutine write_mass

!-----------------------------------------------------------------------
!+
!  interactive prompting for setting up a cylinder
!+
!-----------------------------------------------------------------------
subroutine set_cylinder_interactive(id,master,cylinder,ieos)
 use prompting,     only:prompt
 use setstar_utils, only:profile_opt,need_inputprofile
 use units,         only:in_solarm,in_solarr,udist,umass
 !use units,         only:in_solarl,unit_luminosity
 use physcon,       only:solarr,solarm,solarl
 type(cylinder_t), intent(out)   :: cylinder
 integer,      intent(in)    :: id,master
 integer,      intent(inout) :: ieos
 real :: mcylinder_msun,rcylinder_rsun,zcylinder_rsun

 ! set defaults
 if (cylinder%iunits==0) then
    call set_defaults_cylinder(cylinder)
 else
   call set_defaults_cylinder_solar(cylinder)
   mcylinder_msun = real(in_solarm(cylinder%mcylinder))
   rcylinder_rsun = real(in_solarr(cylinder%rcylinder))
   zcylinder_rsun = real(in_solarr(cylinder%zcylinder))
 endif

 ! Select cylinder & set default values
 write(*,"(i2,')',1x,a)") 1, 'Uniform density profile     '
 write(*,"(i2,')',1x,a)") 2, 'Gaussian density profile     '

 call prompt('Enter which density profile to use',cylinder%iprofile,1,2)
 call prompt('Enter sigma (=0 if uniform density profile)',cylinder%gauss_sigma)
 if (id==master) write(*,"('Setting up ',a)") trim(profile_opt(cylinder%iprofile))

 ! resolution
 if (cylinder%iprofile > 0) then
    call prompt('Enter the approximate number of particles in the sphere ',cylinder%np,0)
 endif

 ! cylinder properties
 if (cylinder%iunits==0) then
   call prompt('Enter the mass of the cylinder',cylinder%mcylinder)
   call prompt('Enter the radius of the cylinder',cylinder%rcylinder)
   call prompt('Enter the height of the cylinder',cylinder%zcylinder)
 else
   call prompt('Enter the mass of the cylinder',mcylinder_msun)
   call prompt('Enter the radius of the cylinder',rcylinder_rsun)
   call prompt('Enter the height of the cylinder',zcylinder_rsun)
   cylinder%mcylinder = mcylinder_msun*real(solarm/umass)
   cylinder%rcylinder = rcylinder_rsun*real(solarr/udist)
   cylinder%zcylinder = zcylinder_rsun*real(solarr/udist)
 endif

 ! shift and/or rotate cylinder
 call prompt('Shift the center of the cylinder?',cylinder%ishift)
 if (cylinder%ishift) then
   call prompt('Enter the rotation angle of the cylinder along "x" axis (degrees)',cylinder%rotate_theta)
   call prompt('Enter the rotation angle of the cylinder along "z" axis (degrees)',cylinder%rotate_phi)
   call prompt('Enter the "x" coordinate of the shifted cylinder center',cylinder%xshift)
   call prompt('Enter the "y" coordinate of the shifted cylinder center',cylinder%yshift)
   call prompt('Enter the "z" coordinate of the shifted cylinder center',cylinder%zshift)
 endif

end subroutine set_cylinder_interactive



!-----------------------------------------------------------------------
!+
!  read setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine read_options_cylinder(cylinder,ieos,db,nerr,label)
 use infile_utils,  only:inopts,read_inopt
 use setstar_utils, only:need_inputprofile,nprofile_opts
 use units,         only:umass,udist
 !use units,         only:unit_luminosity
 use physcon,       only:solarm,solarr
 !use physcon,       only:solarl
 type(cylinder_t),              intent(out)   :: cylinder
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: ieos
 integer,                   intent(inout) :: nerr
 character(len=*),          intent(in), optional :: label
 character(len=10) :: c
 real :: mcylinder_msun,rcylinder_rsun,zcylinder_rsun
 real :: xshift_rsun,yshift_rsun,zshift_rsun


 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))
 cylinder%label = trim(c)
 call read_inopt(cylinder%iunits,'iunits'//trim(c),db,errcount=nerr)

 ! set defaults
 if (cylinder%iunits==0) then
   call set_defaults_cylinder(cylinder)
 else
   call set_defaults_cylinder_solar(cylinder)
 endif

 call read_inopt(cylinder%iprofile,'iprofile'//trim(c),db,errcount=nerr,min=0,max=nprofile_opts)

 ! cylinder properties
 if (cylinder%iunits==0) then
   call read_inopt(cylinder%mcylinder,'mcylinder'//trim(c),db,errcount=nerr,min=0.)
   call read_inopt(cylinder%rcylinder,'rcylinder'//trim(c),db,errcount=nerr,min=0.)
   call read_inopt(cylinder%zcylinder,'zcylinder'//trim(c),db,errcount=nerr,min=0.)
 else
   call read_inopt(mcylinder_msun,'mcylinder'//trim(c),db,errcount=nerr,min=0.)
   call read_inopt(rcylinder_rsun,'rcylinder'//trim(c),db,errcount=nerr,min=0.)
   call read_inopt(zcylinder_rsun,'zcylinder'//trim(c),db,errcount=nerr,min=0.)
   cylinder%mcylinder = mcylinder_msun*real(solarm/umass)
   cylinder%rcylinder = rcylinder_rsun*real(solarr/udist)
   cylinder%zcylinder = zcylinder_rsun*real(solarr/udist)
 endif
 ! set temperature
 if ((ieos==12) .or. do_radiation) call read_inopt(cylinder%initialtemp,'initialtemp',db,errcount=nerr)

 ! cylinder density profile
 call read_inopt(cylinder%iprofile,'iprofile',db,errcount=nerr)
 call read_inopt(cylinder%gauss_sigma,'gauss_sigma',db,errcount=nerr)

 ! shift and/or rotate cylinder
 call read_inopt(cylinder%ishift,'shift cylinder',db,errcount=nerr)
 if (cylinder%ishift) then
   call read_inopt(cylinder%rotate_theta,'rotate_theta',db,errcount=nerr)
   call read_inopt(cylinder%rotate_phi,'rotate_phi',db,errcount=nerr)
   if (cylinder%iunits==0) then
     call read_inopt(cylinder%xshift,'xshift',db,errcount=nerr)
     call read_inopt(cylinder%yshift,'yshift',db,errcount=nerr)
     call read_inopt(cylinder%zshift,'zshift',db,errcount=nerr)
   else
     call read_inopt(xshift_rsun,'xshift',db,errcount=nerr)
     call read_inopt(yshift_rsun,'yshift',db,errcount=nerr)
     call read_inopt(zshift_rsun,'zshift',db,errcount=nerr)
     cylinder%xshift = xshift_rsun*real(solarm/umass)
     cylinder%yshift = yshift_rsun*real(solarm/umass)
     cylinder%zshift = zshift_rsun*real(solarm/umass)
   endif
 endif

 ! resolution
 if (cylinder%iprofile > 0 .and. (len_trim(c)==0 .or. c(1:1)=='1')) then
    call read_inopt(cylinder%np,'np'//trim(c),db,errcount=nerr,min=10)
 endif

end subroutine read_options_cylinder




!-----------------------------------------------------------------------
!+
!  write setupfile options needed for a cylinder
!+
!-----------------------------------------------------------------------
subroutine write_options_cylinder(cylinder,iunit,ieos,label)
 use infile_utils,  only:write_inopt,get_optstring
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile
 use units,         only:in_solarm,in_solarr
 use physcon,       only:solarr,solarm
 type(cylinder_t),     intent(in) :: cylinder
 integer,          intent(in) :: iunit,ieos
 character(len=*), intent(in), optional :: label
 !real :: mcylinder_msun,rcylinder_rsun,zcylinder_rsun
 !real :: xshift_rsun,yshift_rsun,zshift_rsun
 character(len=120) :: string
 character(len=10) :: c

 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))

 write(iunit,"(/,a)") '# options for cylinder '//trim(c)
 call get_optstring(nprofile_opts,profile_opt,string,4)

 ! options for cylinder
 call write_inopt(cylinder%np,'np','Number of particles',iunit)
 if (cylinder%iunits==0) then
   call write_inopt(cylinder%mcylinder,'mcylinder','Mass of the cylinder (code units)',iunit)
   call write_inopt(cylinder%rcylinder,'rcylinder','Radius of the cylinder (code units)',iunit)
   call write_inopt(cylinder%zcylinder,'zcylinder','Height of the cylinder along z-axis (code units)',iunit)
 else
   call write_inopt(cylinder%mcylinder/solarm,'mcylinder','Mass of the cylinder (M_sun)',iunit)
   call write_inopt(cylinder%rcylinder/solarr,'rcylinder','Radius of the cylinder (R_sun)',iunit)
   call write_inopt(cylinder%zcylinder/solarr,'zcylinder','Height of the cylinder along z-axis (R_sun)',iunit)
   ! set temperature
   if ((ieos==12) .or. do_radiation) then
     call write_inopt(cylinder%initialtemp,'initialtemp','Initial temperature of the cylinder (K)',iunit)
   endif
 endif


 ! options for cylinder transverse density profile
 write(iunit,"(/,a)") '# options for cylinder transverse density profile'
 call write_inopt(cylinder%iprofile,'iprofile','Choose density profile: 1=uniform, 2=Gaussian',iunit)
 call write_inopt(cylinder%gauss_sigma,'gauss_sigma','Standard deviation of the density profile',iunit)

 !options for shifting and/or rotating cylinder
 write(iunit,"(/,a)") '# options for shifting and/or rotating cylinder'
 call write_inopt(cylinder%ishift,'shift cylinder','Shift cylinder automatically during setup',iunit)
 if (cylinder%ishift) then
   call write_inopt(cylinder%rotate_theta,'rotate_theta','Rotation angle of the cylinder along "x" axis (degrees)',iunit)
   call write_inopt(cylinder%rotate_phi,'rotate_phi','Rotation angle of the cylinder along "z" axis (degrees)',iunit)
   if (cylinder%iunits==0) then
     call write_inopt(cylinder%xshift,'xshift','"x" coordinate of the shifted cylinder center',iunit)
     call write_inopt(cylinder%yshift,'yshift','"y" coordinate of the shifted cylinder center',iunit)
     call write_inopt(cylinder%zshift,'zshift','"z" coordinate of the shifted cylinder center',iunit)
   else
     call write_inopt(cylinder%xshift/solarm,'xshift','"x" coordinate of the shifted cylinder center (M_sun)',iunit)
     call write_inopt(cylinder%yshift/solarr,'yshift','"y" coordinate of the shifted cylinder center (R_sun)',iunit)
     call write_inopt(cylinder%zshift/solarr,'zshift','"z" coordinate of the shifted cylinder center (R_sun)',iunit)
   endif
 endif



end subroutine write_options_cylinder


end module setcylinder

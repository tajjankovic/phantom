!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setstream
!
! General routine for setting up a stream (cylinder) from a 1D profile
! This is the main functionality from setup_stream but in a single routine
! that can also be called from other setups.
!
! :References: None
!
! :Owner: Taj Jankovič
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dim, eos, extern_densprofile, infile_utils,
!   io, mpiutils, part, physcon, prompting, radiation_utils, relaxstar,
!   setstar_utils, unifdis, units, vectorutils
!
 use unifdis,    only:set_unifdis,mask_prototype,mask_true
 use physcon,    only:pi
 use stretchmap, only:rho_func

 implicit none

 !
 ! define a data type with all the options needed
 ! to setup stream (these are per-stream, not per-simulation options)
 !
 type stream_t
    integer :: np
    real :: gamma
    real :: rstream1
    real :: zstream1
    real :: rstream2
    real :: zstream2
    real :: inclination
    real :: yshift
    real :: vinj1
    real :: vinj2
    real :: offset
    real :: dtinject
    character(len=120) :: inputfile1 ! outputfilename is the path to the cored profile
    character(len=120) :: inputfile2 ! outputfilename is the path to the cored profile
 end type stream_t

 public :: stream_t
 public :: set_defaults_stream
 public :: write_options_stream,read_options_stream,set_stream_interactive


 real,    private :: tol_ekin = 1.e-7 ! criteria for being converged
 integer, private :: maxits = 5000
 real,    private :: gammaprev,hfactprev,mass1prev
 integer, private :: ieos_prev
 integer, public :: ierr_setup_errors = 1, &
                    ierr_notconverged = 2
 private

contains

!--------------------------------------------------------------------------
!+
!  default options for a particular stream
!  (see also set_defaults_given_profile which selects defaults
!   based on the value of iprofile)
!+
!--------------------------------------------------------------------------
subroutine set_defaults_stream(stream)
 !use units,   only:udist,umass
 !use physcon, only:solarm,solarr
 type(stream_t), intent(out) :: stream
 stream%np            = 1000
 stream%gamma    = 4./3.
 !stream%mstream       = 1.0*real(solarm/umass)
 !stream%rstream       = 1.0*real(solarr/udist)
 !stream%zstream       = 1.0*real(solarr/udist)
 stream%rstream1       = 1.0
 stream%zstream1       = 2.0
 stream%rstream2       = 1.0
 stream%zstream2       = 2.0
 stream%inclination   = 180.0
 stream%yshift        = 4.
 stream%vinj1         = 1.0
 stream%vinj2         = 1.0
 stream%offset        = 0.0
 stream%dtinject        = 0.01
 stream%inputfile1  = 'cylinder1.ascii'
 stream%inputfile2  = 'cylinder1.ascii'
end subroutine set_defaults_stream


!----------------------------------------------------------------
!+
!  restore previous settings
!+
!----------------------------------------------------------------
subroutine restore_original_options(i1,npart)
 use eos,     only:ieos,gamma
 use damping, only:damp
 use options, only:idamp
 use part,    only:hfact,vxyzu,gr
 use externalforces, only:mass1
 integer, intent(in) :: i1,npart

 gamma = gammaprev
 hfact = hfactprev
 ieos  = ieos_prev
 idamp = 0
 damp = 0.
 vxyzu(1:3,i1+1:npart) = 0.
 if (gr) mass1 = mass1prev

end subroutine restore_original_options


!-----------------------------------------------------------------------
!+
!  interactive prompting for setting up a stream
!+
!-----------------------------------------------------------------------
subroutine set_stream_interactive(id,master,stream,ieos,polyk)
 use prompting,     only:prompt
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar
 !use units,         only:in_solarm,in_solarr,in_solarl,udist,umass,unit_luminosity
 use physcon,       only:solarr,solarm,solarl
 type(stream_t), intent(out)   :: stream
 integer,      intent(in)    :: id,master
 integer,      intent(inout) :: ieos
 real,         intent(inout) :: polyk
 !real :: mstream_msun,rstream_rsun,zstream_rsun
 !real :: mstream,rstream,zstream,rotation_angle,xshift,yshift,zshift,dotM,vinj

 ! set defaults
 call set_defaults_stream(stream)
 !mstream_msun = real(in_solarm(stream%Mstream))
 !rstream_rsun = real(in_solarr(stream%rstream))
 !zstream_rsun = real(in_solarr(stream%zstream))
 !mstream = real(stream%mstream)
 !rstream = real(stream%rstream)
 !zstream = real(stream%zstream)

 ! Select stream & set default values
 !write(*,"(i2,')',1x,a)") 1, 'Uniform density profile     '
 !write(*,"(i2,')',1x,a)") 2, 'Gaussian density profile     '
 !if (id==master) write(*,"('Setting up ',a)") trim(profile_opt(stream%iprofile))
 ! resolution

 call prompt('Enter the number of particles in a cylinder',stream%np,1)
 !call prompt('Enter the number of particles in a cylinder',stream%np,1)

 call prompt('Enter the adiabatic index gamma',stream%gamma)

 !
 ! set default file output parameters
 !
 !call set_defaults_given_profile(stream%iprofile,stream%input_profile,&
  !                               need_iso,ieos,stream%mstream,polyk)


 ! stream properties
 call prompt('Enter the radius of the 1st stream (code units)',stream%rstream1)
 call prompt('Enter the length of the 1st stream (code units)',stream%zstream1)
 call prompt('Enter the radius of the 2nd stream (code units)',stream%rstream2)
 call prompt('Enter the length of the 2nd stream (code units)',stream%zstream2)
 call prompt('Enter inclination between streams (degrees)',stream%inclination,0.)
 call prompt('Enter the "y" coordinate of the shifted cylinder center (code units)',stream%yshift)
 call prompt('Enter the injection speed of the first stream (code units)',stream%vinj1)
 call prompt('Enter the injection speed of the second stream(code units)',stream%vinj2)
 call prompt('Enter offset between streams (units of rstream)',stream%offset)
 call prompt('Enter delta t injection (units of rstream)',stream%dtinject)
 call prompt('Enter input file name of the 1st cylinder:',stream%inputfile1)
 call prompt('Enter input file name of the 2nd cylinder:',stream%inputfile2)


end subroutine set_stream_interactive

!-----------------------------------------------------------------------
!+
!  write setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine write_options_stream(stream,iunit)
 use infile_utils,  only:write_inopt,get_optstring
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar
 use units,         only:in_solarm,in_solarr,in_solarl
 type(stream_t),     intent(in) :: stream
 integer,          intent(in) :: iunit
 character(len=120) :: string
 character(len=10) :: c

 ! append optional label e.g. '1', '2'
 c = ''

 write(iunit,"(/,a)") '# options for stream '//trim(c)
 call get_optstring(nprofile_opts,profile_opt,string,4)
 !call write_inopt(stream%iprofile,'iprofile'//trim(c),'0=Sink,'//trim(string(1:40)),iunit)

 call write_inopt(stream%np,'np'//trim(c),'number of particles',iunit)

 !if (stream%iprofile > 0 .and. (len_trim(c)==0 .or. c(1:1)=='1')) then
 !endif

end subroutine write_options_stream

!-----------------------------------------------------------------------
!+
!  read setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine read_options_stream(stream,ieos,polyk,db,nerr)
 use infile_utils,  only:inopts,read_inopt
 use setstar_utils, only:need_inputprofile,need_rstar,nprofile_opts
 !use units,         only:umass,udist,unit_luminosity
 use physcon,       only:solarm,solarr,solarl
 type(stream_t),              intent(out)   :: stream
 type(inopts), allocatable, intent(inout) :: db(:)
 !integer,                   intent(out)   :: need_iso
 integer,                   intent(inout) :: ieos
 real,                      intent(inout) :: polyk
 integer,                   intent(inout) :: nerr
 character(len=10) :: c
 !real :: mstream_msun,rstream_rsun

 ! set defaults
 call set_defaults_stream(stream)


 !call read_inopt(stream%iprofile,'iprofile'//trim(c),db,errcount=nerr,min=0,max=nprofile_opts)

! if (need_inputprofile(stream%iprofile)) then
!    call read_inopt(stream%input_profile,'input_profile'//trim(c),db,errcount=nerr)
 !endif

 ! resolution
 !if (stream%iprofile > 0 .and. (len_trim(c)==0 .or. c(1:1)=='1')) then
  !  call read_inopt(stream%np,'np'//trim(c),db,errcount=nerr,min=10)
 !endif
end subroutine read_options_stream


end module setstream

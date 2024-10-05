!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! moddump to add a sink particle
!
! :References: None
!
! :Owner: Taj Jankovič & Aleksej Jurca
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part, physcon, prompting, setflyby, units,
!   vectorutils
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:nptmass,xyzmh_ptmass,vxyz_ptmass!,igas,ihacc,ihsoft
 use prompting,         only:prompt
 use units,             only:umass,utime,udist,print_units
 use physcon,           only:au,solarm,pi,years
 use centreofmass,      only:reset_centreofmass,get_centreofmass
 use vectorutils,       only:rotatevec
 use setflyby,          only:get_T_flyby
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: msink, rsink, z0, incl, vz0
 real                   :: xp(3), vp(3), rot_axis(3)

 ! sink particle properties
 msink = 6.7d-5 ! mass of the perturber
 z0 = 10. ! initial distance of the perturber
 incl = 0. ! inclination
 vz0   = 1.0
 rsink = 0.5
 print *,'-----Sink particle mass is:', msink
 print *,'-----Sink particle radius is:', rsink
 print *,'-----Sink particle distance alog z-direction is:', z0
 print *,'-----Sink particle velocity is:', vz0

 nptmass = nptmass + 1

 ! perturber initial position
 xp = (/0.0,0.0,z0/)

 ! perturber initial velocity
 vp  = (/0.0,0.0,-vz0/)

 ! assigning position, velocity and accretion radius to the sink particle
 xyzmh_ptmass(1:3,nptmass)   = xp
 xyzmh_ptmass(4,nptmass)     = msink
 xyzmh_ptmass(5,nptmass)     = rsink
 xyzmh_ptmass(6,nptmass)     = rsink ! softening length
 vxyz_ptmass(:,nptmass)      = vp

! call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 return

end subroutine modify_dump

end module moddump


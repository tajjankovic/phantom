!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Convert non-radiation dump (assuming LTE, ieos=12) to radiation dump
!
! :References: None
!
! :Owner: Taj Jankovič
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, eos_idealplusrad, eos_mesa, io,
!   mesa_microphysics, part, radiation_utils, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use partinject, only:add_or_update_particle
 use random,  only:ran2
 use units,            only:unit_density,unit_opacity,unit_ergg,get_radconst_code,get_kbmh_code,get_c_code
 use units,            only:utime
 use dim,              only:do_radiation
 use io,               only:fatal
 use eos,              only:gmw,gamma,get_cv
 use eos_idealplusrad, only:get_idealplusrad_temp
 use eos_mesa,         only:init_eos_mesa
 use part,             only:igas,idust,rad,iradxi,ikappa,rhoh,radprop,ithick
 use radiation_utils,  only:radiation_and_gas_temperature_equal,ugas_from_Tgas
 use timestep,          only:tmax,dtmax
 use mesa_microphysics,only:get_kappa_mesa
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 !real, dimension(:,:), allocatable :: xyzh_add,vxyzu_add(:,:)
 integer                :: i,ipart,iseed,npart0
 real                   :: xi0,ugas0,Mstar,Rstar, dist, tdyn, rhoi, Tbath
 real                   :: cv, a_code, rho_mean, Teq,kbmh_code, erad0,hpart
 real                   :: V0, Eej, c_code,erad0_QPE
 logical                :: thermal_equil,uniform

 if (.not. do_radiation) call fatal("moddump_LTE_to_rad","Not compiled with radiation")

 a_code = get_radconst_code()
 c_code = get_c_code()
 kbmh_code = get_kbmh_code()
 print*,'gamma=',gamma,', mu=',gmw,', a_code=',a_code,', kb/mh_code=',kbmh_code

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!! Initialize values !!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 thermal_equil = .true. ! true, false
 uniform = .true. ! true -> uniform radial rho profile, false -> polytropic rho profile
 Mstar = 0.0
 Rstar = 0.0
 ! Sum total mass Mstar and find maximum radial distance Rstar
 do i = 1, npart
    Mstar = Mstar + massoftype(igas)
    
    ! Calculate radial distance from origin for each particle
    dist = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
    
    ! Update Rstar with the maximum distance found
    if (dist > Rstar) then
        Rstar = dist
    end if
 end do
 tdyn = sqrt(4*3.14*Rstar**3/3/Mstar) ! sets timescale
 rho_mean = 3*Mstar/(4*3.14*Rstar**3) ! for xi
 if (uniform) then
   rho_mean = rhoh(xyzh(4, 1), massoftype(igas)) ! if calculate as rho_mean, then a bit wrong
 endif

 ! set timesteps
 tmax = 3600/utime ! 1h in code units
 dtmax = tmax/1000.
 print*,'Mstar is:',Mstar
 print*,'Rstar is:',Rstar
 print*,'tdyn is:',tdyn
 print*,'rho_mean is:',rho_mean

 ! set values expected for QPEs (Linial 2023)
 V0 = 4.*3.14/3.*Rstar**3. ! volume of the shocked gas
 Eej = 0.5*Mstar*(2**0.5*0.1*c_code)**2. ! Energy of ejecta assuming stellar velocity 0.1c
 erad0_QPE = 1e-2*Eej/V0 ! rad. energy density
 print*,'Injected energy is:',Eej
 print*,'Injected erad is:',erad0_QPE
 print*,'Shocked volume is',V0

 do i=1,npart
    ! assing u and xi
    if (uniform) then
      if (thermal_equil) then
         !!!!! get T from Prad=1/3*erad=a/3*T^4 !!!!!!
         Teq = (erad0_QPE/a_code)**(1./4.)
         xi0 = a_code*Teq**4/rho_mean
         ugas0 = kbmh_code/(gmw*(gamma-1.))*Teq
         rad(iradxi,i) = xi0
         vxyzu(4,i) = ugas0
         if (i==1) then
            print*,'Sphere is isothermal, uniform rho profile, Teq:',Teq
            print*,'xi0=',xi0,', ugas0=',ugas0
         endif
      else
         !!!!! set initial xi (E/mass) from erad (E/V)=const;
         !!!!! assume erad = 0.01*Ebind/V= 3*G*M^2/(4*pi*R^4)
         erad0 = 0.01*3*Mstar**2/(4*3.14*Rstar**4) ! G=1
         xi0 = erad0/rho_mean
         ugas0 = xi0/(300*(gamma-1.)) ! from Pgas=100*Prad
         rad(iradxi,i) = xi0
         vxyzu(4,i) = ugas0
         if (i==1) then
            print*,'Sphere is not isothermal, uniform rho profile, Prad=100Pgas:',Teq
            print*,'xi0=',xi0,', ugas0=',ugas0
         endif
         print*,'Not working!'
      endif
    else
      if (thermal_equil) then
         !!!!! T(300*k/(mu*mp*a)*rho_mean)**1/3 - from Prad=100*Pgas !!!!!!
         rhoi = rhoh(xyzh(4, i), massoftype(igas)) ! no longer uniform -> calculate for every rho
         Teq = (300*kbmh_code/gmw/a_code*rhoi)**(1./3.)
         xi0 = a_code*Teq**4/rhoi
         ugas0 = kbmh_code/(gmw*(gamma-1.))*Teq
         rad(iradxi,i) = xi0
         vxyzu(4,i) = ugas0
         if (i==1) then
            print*,'Sphere is in thermal equil., poly rho profile:',Teq
         endif
         print*,'Not working!'
      else
         !!!!! set initial xi (E/mass) from erad (E/V)=const;
         !!!!! assume erad = 0.01*Ebind/V= 3*G*M^2/(4*pi*R^4)
         erad0 = 0.01*3*Mstar**2/(4*3.14*Rstar**4) ! G=1
         rhoi = rhoh(xyzh(4, i), massoftype(igas)) ! no longer uniform -> calculate for every rho
         xi0 = erad0/rhoi
         ugas0 = xi0/(300*(gamma-1.)) ! from Pgas=100*Prad
         rad(iradxi,i) = xi0
         vxyzu(4,i) = ugas0
         if (i==1) then
            print*,'Sphere is in not in thermal equil., poly rho profile:',Teq
         endif
         print*,'Not working!'
      endif
    endif


    ! method 4: Trad=Tgas
    !rad(iradxi, i) = a_code * T0**4. * (rhoi / rho0)**(4.0 / 3.0) / rhoi  ! smooth transition
    !cv = get_cv(rhoi,vxyzu(4,i),0) ! cv_type=0
    !vxyzu(4,i) = cv * T0

    radprop(ithick,i) = 1. ! assume optically thick gas
 
 enddo
 
end subroutine modify_dump

end module moddump



!=================================================================

function FouriercomponentCoulombinteractionFZ ( Atmdc, Qx, Qy )
!~ function FouriercomponentCoulombinteractionFZ ( Atmdc, kappavalue, Qx, Qy )
use moduleone
implicit none
real(kind=DP):: FouriercomponentCoulombinteractionFZ
real(kind=DP):: Qx, Qy, QQ, kappavalue, Atmdc, carboncarbondistance


QQ=sqrt(Qx**2+Qy**2)

carboncarbondistance=Atmdc/(sq3)  !! Q^(-1) in unit of ( m )


FouriercomponentCoulombinteractionFZ= 4.0d0*pi*Fine_structure_constant*Planck_constant*speedoflight/( QQ*carboncarbondistance )   !! System International 
!~ FouriercomponentCoulombinteractionFZ= 4.0d0*pi*Fine_structure_constant*Planck_constant*speedoflight/( QQ*kappavalue*carboncarbondistance )   !! System International 
end 


!===========================================================================

subroutine Chargedistribution( RNprime, Doping_way, Nlayer, RNprimelayer )
use moduleone
Implicit real (kind=DP) (a-h,o-z) 
integer(kind=4), intent(in):: Doping_way, Nlayer

real(kind=DP), intent(in) :: RNprime
real(kind=DP), dimension(1:3*Nlayer), intent(out) :: RNprimelayer
RNprimelayer(:)=0.d0

if      ( Doping_way==1 ) then
    RNprimelayer(1)=RNprime*0.30d0
    RNprimelayer(2)=RNprime*0.35d0
    RNprimelayer(3)=RNprime-RNprimelayer(2)-RNprimelayer(1)
else if ( Doping_way==2 ) then
    RNprimelayer(1)=RNprime*0.20d0
    RNprimelayer(2)=RNprime*0.60d0
    RNprimelayer(3)=RNprime*0.20d0
else if ( Doping_way==3 ) then
    print*, "This part is not complete. "

else
    print*, "what's the value of Doping_way?"
    stop
end if

end	subroutine Chargedistribution


!===========================================================================

subroutine Kappavalue ( Doping_way, kappa )
use moduleone
Implicit real (kind=DP) (a-h,o-z) 
integer(kind=4), intent(in):: Doping_way
real(kind=DP), intent(inout)::   kappa
real(kind=DP) :: A

A=kappa

if      ( Doping_way==1 ) then
    kappa=(A+1.0)/2.0_DP

else if ( Doping_way==2 ) then
    kappa= A
else if ( Doping_way==3 ) then
    print*, "This part is not complete. "

else
    print*, "what's the value of Doping_way?"
    stop
end if

end	subroutine Kappavalue



!===========================================================================
subroutine Counter(A, probability, deltax, Ndimension,  Densityofstate)
use moduleone
Implicit real (kind=DP) (a-h,o-z) 

real(kind=DP) :: A,  probability
real(kind=DP), dimension(0:Ndimension-1), intent(INout) :: Densityofstate

!~ real(kind=DP), dimension(0:Ndeltax) :: Densityofstate

Densityofstate(:)=0.d0

If ( (A-Al) .Le. 0.d0 ) then     
    print*     
    print*, "A = ", A     
    print*     
    print*, "please have a check the first parameter in subrountine_Counter( .A, ... )"     
    stop 
end if

kc=Nint((A- Al)/deltaX) 

Densityofstate(kc)=Densityofstate(kc)+probability

end	subroutine Counter

!===========================================================================

subroutine   area( Ndimension, Densityofstate, N0, Nend, S )
use moduleone
Implicit real (kind=DP) (a-h,o-z) 
real(kind=DP), dimension(0:Ndimension-1), intent(IN) :: Densityofstate

integer(kind=4):: N0, Nend

real(kind=DP) :: S

   
         S=0.0_DP 
              do i=N0+1, Nend-1
                  S=S+Densityofstate(i)
              end do
	 S=S+ (Densityofstate(N0)+Densityofstate(Nend))/2.0_DP

end  subroutine area

!===========================================================================

subroutine  paoguang( Ndimension, Densityofstate )
use moduleone
Implicit real (kind=DP) (a-h,o-z) 
real(kind=DP), dimension(0:Ndimension-1), intent(INout) :: Densityofstate
real(kind=DP) :: delta

    delta=1e-6
    do i= 0, Ndimension-1
        y=Densityofstate(i)
        if ( abs(y) .le. 1e-6 ) then
            Densityofstate(i)=0.d0
        end if        
    end do

end  subroutine paoguang

!===========================================================================
  subroutine Rnorm(nstates,psi_R,r0)
  use moduleone
        implicit real(kind=DP) (a-h,o-z)
        real(kind=DP) psi_R(0:nstates-1)

        r0=0.0
        do j=0,nstates-1
        r0=r0+psi_R(j)
        enddo
end  subroutine Rnorm


!===========================================================================

subroutine Vfun( kappa, nprime, carriernumber, distanceofinterlayer, Sunitcell,  Vdeta)
use moduleone
implicit none
 
 integer(KIND=4), intent(in):: carriernumber
 
 real(kind=DP), intent(in):: nprime, Sunitcell, kappa
 
 real(kind=DP):: Vdeta,  distanceofinterlayer
 
  !~ electron=1.60217653*10.0**(-19) !!   unit ( C )
  !~ kappa = (3.90d0+1.0d0)/2.0d0 ;     
  !~ epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V

Vdeta=0.d0
Vdeta=electron*distanceofinterlayer*carriernumber*nprime/(epsilong*kappa*Sunitcell)    ! UNIT ~ eV

!~ AA    Vdeta=electron*distanceofinterlayer*nprime*carriernumber*10.d0**(-17)/(epsilong*kappa*Sunitcell) ! UNIT ~ eV  
end subroutine Vfun


!===========================================================================

subroutine Vfunx( kappa, panduanway ,nprime, carriernumber, distanceofinterlayer, Sunitcell,  Vdeta)
use moduleone
implicit none
 
 integer(KIND=4), intent(in):: carriernumber, panduanway
 
 real(kind=DP), intent(in):: nprime, Sunitcell, kappa
 
 real(kind=DP):: Vdeta,  distanceofinterlayer
 
 
 if ( panduanway==1 ) then
        Vdeta=electron*distanceofinterlayer*carriernumber*nprime/(epsilong*kappa*Sunitcell)    ! UNIT ~ eV
    else if ( panduanway==2 ) then
        Vdeta=0.5d0*electron*distanceofinterlayer*carriernumber*nprime/(epsilong*kappa*Sunitcell)    ! UNIT ~ eV
 end if
 
  !~ electron=1.60217653*10.0**(-19) !!   unit ( C )
  !~ kappa = (3.90d0+1.0d0)/2.0d0 ;     
  !~ epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V
    
Vdeta=electron*distanceofinterlayer*carriernumber*nprime/(epsilong*kappa*Sunitcell)    ! UNIT ~ eV

!~ AA    Vdeta=electron*distanceofinterlayer*nprime*carriernumber*10.d0**(-17)/(epsilong*kappa*Sunitcell) ! UNIT ~ eV  
end subroutine Vfunx


!===========================================================================

subroutine Efun(kappa, nprime, Efiled)    !! nprime  UNIT   (     10^{13} cm^{-2}   )
use moduleone
 implicit none
 real(kind=DP), intent(in):: kappa, nprime     !! UNIT   (     10^{13} cm^{-2}   )
 
 real(kind=DP), intent(out):: Efiled
  
  !~ electron=1.60217653*10.0**(-19) !!   unit ( C )

  !~ kappa = (3.90d0+1.0d0)/2.0d0 ;     

  !~ epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V
  
   !! nprime  UNIT   (     10^{13} cm^{-2}   )
  
  
Efiled=electron*nprime*10.d0**(8)/(2.0d0*epsilong*kappa)    !! UNIT  (    V/nm   )


end subroutine Efun


!===========================================================================

subroutine EfunIIT( kappa, nprime, Efiled )    !! nprime  UNIT   (     10^{13} cm^{-2}   )
use moduleone
 implicit none
 real(kind=DP), intent(in):: nprime, kappa     !! UNIT   (     10^{13} cm^{-2}   )
 
 real(kind=DP), intent(out):: Efiled
  
  !~ electron=1.60217653*10.0**(-19) !!   unit ( C )

  !~ kappa = (3.90d0+1.0d0)/2.0d0 ;     

  !~ epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V
  
  
Efiled=electron*nprime*10.d0**(8)/(2.0d0*epsilong*kappa)    !! UNIT  (    V/nm   )


end subroutine EfunIIT

!===========================================================================

function FermiDistribution (E,Ef,Beta)
use moduleone
    implicit none
    real(kind=DP):: FermiDistribution,E,Ef,Beta
    integer imax
    imax=100
    FermiDistribution=Beta*(E-Ef)
    if (FermiDistribution>imax) then
    FermiDistribution=0.
    else if (FermiDistribution<-imax) then
    FermiDistribution=1.
    else
    FermiDistribution=1.+dexp(FermiDistribution)
    FermiDistribution=1./FermiDistribution
    end if
    end


!=================================================================

function FouriercomponentCoulombinteractionFME ( Atmdc, kappavalue, Qx, Qy )
use moduleone
implicit none
    real(kind=DP):: FouriercomponentCoulombinteractionFME, Qx, Qy, QQ, kappavalue
    real(kind=DP):: electron_hat_2, Atmdc, AtmdcII
    real(kind=DP):: carboncarbondistance

    QQ=sqrt(Qx**2+Qy**2)
    !! shengjun note 78page due to the CGS units
    electron_hat_2=1.4399644d0
    AtmdcII=Atmdc*10.0**(10)       
    carboncarbondistance=AtmdcII/(sq3)     !!Q^(-1) in unit of ( nm )

    FouriercomponentCoulombinteractionFME= 2.0d0*pi*electron_hat_2/(QQ*kappavalue*carboncarbondistance)   !! System International 

end 

!!~ FouriercomponentCoulombinteraction= 2.0d0*pi*electron_hat_2/(Qx*kappa*carboncarbondistance)   !! System International 
!!~ real(kind=DP) :: electron=1.60217653*10.0**(-19) !!   unit ( C )
!!~ real(kind=DP) :: kappa = (3.90d0+1.0d0)/2.0d0 ;   
!!~ real(kind=DP) :: epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V 
! 1.4399644d0/( t(eV) a(nm) )
!=================================================================


!~ function CGSzhuanhuan( Atmdc)
!~ use moduleone
!~ implicit none
!~ real(kind=DP):: CGSzhuanhuan, Atmdc
!~ real(kind=DP):: e_cgs
    !~ Atmdc=Atmdc*10.0**(10)
    !~ carboncarbondistance=Atmdc/(sq3)
    !~ e_cgs  = 4.80320427_dp*10.0_dp**(-10)  ! statcoulomb
    !~ ev_cgs = 1.0_dp/( 299.792458_dp*4.80320427_dp*10.0_dp**(-10) )
    !~ FouriercomponentCoulombinteraction= 2.0d0*pi*electron_hat_2/(Qx*kappa)   !! System International 
!~ end 
!=================================================================





!~ ============================================================
subroutine Bond_integrals ( Tswitch, Delta0, Delta1, Delta2, DeltaP, DeltaZ, &  
        & V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, &
		& V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp)
!~ ============================================================
	use moduleone
	use moduletwo
	implicit real (kind=DP) (a-h, o-z)
	
	integer (kind=4), intent(in) :: Tswitch
	
	real (kind=DP) :: Delta0, Delta1, Delta2, DeltaP, DeltaZ
	real (kind=DP) :: V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi
	real (kind=DP) :: V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp

	
	!~ ================================
	!	MoS2
	if (Tswitch==1) then
		Delta0=Delta0_MoS2
		Delta1=Delta1_MoS2
		Delta2=Delta2_MoS2
		DeltaP=DeltaP_MoS2
		DeltaZ=DeltaZ_MoS2
		V_pd_sigma=V_pd_sigma_MoS2
		V_pd_pi=V_pd_pi_MoS2
		V_dd_sigma=V_dd_sigma_MoS2
		V_dd_pi=V_dd_pi_MoS2
		V_dd_delta=V_dd_delta_MoS2
		V_pp_sigma=V_pp_sigma_MoS2
		V_pp_pi=V_pp_pi_MoS2
		V_pp_sigma_interlayer=V_pp_sigma_interlayer_MoS2
		V_pp_pi_interlayer=V_pp_pi_interlayer_MoS2
		SOC_dd=SOC_dd_MoS2
		SOC_pp=SOC_pp_MoS2	
	!	WS2
	else if (Tswitch==2) then
		Delta0=Delta0_WS2
		Delta1=Delta1_WS2
		Delta2=Delta2_WS2
		DeltaP=DeltaP_WS2
		DeltaZ=DeltaZ_WS2
		V_pd_sigma=V_pd_sigma_WS2
		V_pd_pi=V_pd_pi_WS2
		V_dd_sigma=V_dd_sigma_WS2
		V_dd_pi=V_dd_pi_WS2
		V_dd_delta=V_dd_delta_WS2
		V_pp_sigma=V_pp_sigma_WS2
		V_pp_pi=V_pp_pi_WS2
		V_pp_sigma_interlayer=V_pp_sigma_interlayer_WS2
		V_pp_pi_interlayer=V_pp_pi_interlayer_WS2
		SOC_dd=SOC_dd_WS2
		SOC_pp=SOC_pp_WS2	
	!	MoSe2
	else if (Tswitch==3) then
		Delta0=Delta0_MoSe2
		Delta1=Delta1_MoSe2
		Delta2=Delta2_MoSe2
		DeltaP=DeltaP_MoSe2
		DeltaZ=DeltaZ_MoSe2
		V_pd_sigma=V_pd_sigma_MoSe2
		V_pd_pi=V_pd_pi_MoSe2
		V_dd_sigma=V_dd_sigma_MoSe2
		V_dd_pi=V_dd_pi_MoSe2
		V_dd_delta=V_dd_delta_MoSe2
		V_pp_sigma=V_pp_sigma_MoSe2
		V_pp_pi=V_pp_pi_MoSe2
		V_pp_sigma_interlayer=V_pp_sigma_interlayer_MoSe2
		V_pp_pi_interlayer=V_pp_pi_interlayer_MoSe2
		SOC_dd=SOC_dd_MoSe2
		SOC_pp=SOC_pp_MoSe2	
	!	WSe2
	else if (Tswitch==4) then
		Delta0=Delta0_WSe2
		Delta1=Delta1_WSe2
		Delta2=Delta2_WSe2
		DeltaP=DeltaP_WSe2
		DeltaZ=DeltaZ_WSe2
		V_pd_sigma=V_pd_sigma_WSe2
		V_pd_pi=V_pd_pi_WSe2
		V_dd_sigma=V_dd_sigma_WSe2
		V_dd_pi=V_dd_pi_WSe2
		V_dd_delta=V_dd_delta_WSe2
		V_pp_sigma=V_pp_sigma_WSe2
		V_pp_pi=V_pp_pi_WSe2
		V_pp_sigma_interlayer=V_pp_sigma_interlayer_WSe2
		V_pp_pi_interlayer=V_pp_pi_interlayer_WSe2
		SOC_dd=SOC_dd_WSe2
		SOC_pp=SOC_pp_WSe2
	end if	
	
!~ ============================================
end subroutine Bond_integrals
!~ ============================================	




!~ ============================================================
subroutine PP_Bonds_Family( Tswitch, Sx, Sy, Sz, &  
        & E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
		& E_z_x,  E_z_y , E_z_z)
!~ ============================================================
	use moduleone
	use moduletwo
	implicit real (kind=DP) (a-h, o-z)
	
	integer (kind=4), intent(in) :: Tswitch
	
	real (kind=DP) :: Delta0, Delta1, Delta2, DeltaP, DeltaZ
	real (kind=DP) :: V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi
	real (kind=DP) :: V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp
	real (kind=DP) :: Vn, Vx2, Vy2, Vz2	
	
	real (kind=DP), intent(in) :: Sx, Sy, Sz	

	real (kind=DP), intent(out) :: E_x_x, E_x_y, E_x_z, &
		& E_y_x,  E_y_y , E_y_z ,  E_z_x,  E_z_y , E_z_z

	Vn=sqrt( Sx**2+Sy**2+Sz**2  )
	Vx=  Sx/Vn  ;       Vy= Sy/Vn  ;      Vz= Sz/Vn     
    Vx2=Vx**2
	Vy2=Vy**2
	Vz2=Vz**2
	
	call Bond_integrals ( Tswitch, Delta0, Delta1, Delta2, DeltaP, DeltaZ, &  
        & V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, &
		& V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp)	
	
	E_x_x = Vx2*V_pp_sigma+(1.d0-Vx2)*V_pp_pi
	
	E_x_y = Vx*Vy*( V_pp_sigma - V_pp_pi )

	E_x_z = Vx*Vz*( V_pp_sigma - V_pp_pi )
	
	E_y_x = E_x_y 
	
	E_y_y = Vy2*V_pp_sigma+(1.d0-Vy2)*V_pp_pi
	
	E_y_z = Vy*Vz*( V_pp_sigma - V_pp_pi )	

	E_z_x = E_x_z
	
	E_z_y = E_y_z 
	
	E_z_z = Vz2*V_pp_sigma+(1.d0-Vz2)*V_pp_pi
	
!~ ============================================
end subroutine PP_Bonds_Family
!~ ============================================	


!~ ============================================================
subroutine PPINTERLAYER_Bonds_Family( Tswitch, Sx, Sy, Sz, &  
        & E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
		& E_z_x,  E_z_y , E_z_z)
!~ ============================================================
	use moduleone
	use moduletwo
	implicit real (kind=DP) (a-h, o-z)
	
	integer (kind=4), intent(in) :: Tswitch
	
	real (kind=DP) :: Delta0, Delta1, Delta2, DeltaP, DeltaZ
	real (kind=DP) :: V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi
	real (kind=DP) :: V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp
	real (kind=DP) :: Vn, Vx2, Vy2, Vz2	
	
	real (kind=DP), intent(in) :: Sx, Sy, Sz	

	real (kind=DP), intent(out) :: E_x_x, E_x_y, E_x_z, &
		& E_y_x,  E_y_y , E_y_z ,  E_z_x,  E_z_y , E_z_z

	Vn=sqrt( Sx**2+Sy**2+Sz**2  )
	Vx=  Sx/Vn  ;       Vy= Sy/Vn  ;      Vz= Sz/Vn     
    Vx2=Vx**2
	Vy2=Vy**2
	Vz2=Vz**2
	
	call Bond_integrals ( Tswitch, Delta0, Delta1, Delta2, DeltaP, DeltaZ, &  
        & V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, &
		& V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp)	
	
	E_x_x = Vx2*V_pp_sigma_interlayer+(1.d0-Vx2)*V_pp_pi_interlayer
	
	E_x_y = Vx*Vy*( V_pp_sigma_interlayer - V_pp_pi_interlayer )

	E_x_z = Vx*Vz*( V_pp_sigma_interlayer - V_pp_pi_interlayer )
	
	E_y_x = E_x_y 
	
	E_y_y = Vy2*V_pp_sigma_interlayer+(1.d0-Vy2)*V_pp_pi_interlayer
	
	E_y_z = Vy*Vz*( V_pp_sigma_interlayer - V_pp_pi_interlayer )	

	E_z_x = E_x_z
	
	E_z_y = E_y_z 
	
	E_z_z = Vz2*V_pp_sigma_interlayer+(1.d0-Vz2)*V_pp_pi_interlayer
	
!~ ============================================
end subroutine PPINTERLAYER_Bonds_Family
!~ ============================================	




!~ ============================================================
subroutine PD_Bonds_Family( Tswitch, Sx, Sy, Sz, &  
        & E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
		& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
		& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)
!~ ============================================================
	use moduleone
	use moduletwo
	implicit real (kind=DP) (a-h, o-z)
	
	integer (kind=4), intent(in) :: Tswitch
	
	real (kind=DP) :: Delta0, Delta1, Delta2, DeltaP, DeltaZ
	real (kind=DP) :: V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi
	real (kind=DP) :: V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp
	real (kind=DP) :: Vn, Vx2, Vy2, Vz2	
	
	real (kind=DP), intent(in) :: Sx, Sy, Sz	

	real (kind=DP), intent(out) :: E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2, &
		& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
		& E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2

	Vn=sqrt( Sx**2+Sy**2+Sz**2  )
	Vx=  Sx/Vn  ;       Vy= Sy/Vn  ;      Vz= Sz/Vn     
    Vx2=Vx**2
	Vy2=Vy**2
	Vz2=Vz**2
	
	call Bond_integrals ( Tswitch, Delta0, Delta1, Delta2, DeltaP, DeltaZ, &  
        & V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, &
		& V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp)	
	
	E_x_xy =sq3*Vx2*Vy*V_pd_sigma+Vy*(1-2.0d0*Vx2)*V_pd_pi	
	
	E_x_yz =Vx*Vy*Vz*(sq3*V_pd_sigma-2.0d0*V_pd_pi)	

	E_x_xz =sq3*Vx2*Vz*V_pd_sigma+Vz*(1-2.0d0*Vx2)*V_pd_pi

	E_x_x2y2 =0.5d0*sq3*Vx*(Vx2-Vy2)*V_pd_sigma +Vx*(1.0d0-Vx2+Vy2)*V_pd_pi

	E_x_z2 =Vx*( Vz2-0.5d0*(Vx2+Vy2) )*V_pd_sigma - sq3*Vx*Vz2*V_pd_pi


	E_y_xy =sq3*Vy2*Vx*V_pd_sigma+Vx*(1-2.0d0*Vy2)*V_pd_pi
	
	E_y_yz =sq3*Vy2*Vz*V_pd_sigma+Vz*(1-2.0d0*Vy2)*V_pd_pi	
	
	E_y_xz =Vy*Vx*Vz*(sq3*V_pd_sigma-2.0d0*V_pd_pi)
	
	E_y_x2y2 =0.5d0*sq3*Vy*(Vx2-Vy2)*V_pd_sigma - Vy*(1.0d0+Vx2-Vy2)*V_pd_pi

	E_y_z2 =Vy*( Vz2-0.5d0*(Vx2+Vy2) )*V_pd_sigma - sq3*Vy*Vz2*V_pd_pi
	


	E_z_xy =Vz*Vx*Vy*(sq3*V_pd_sigma-2.0d0*V_pd_pi)
	
	E_z_yz =sq3*Vz2*Vy*V_pd_sigma+Vy*(1-2.0d0*Vz2)*V_pd_pi	

	E_z_xz =sq3*Vz2*Vx*V_pd_sigma+Vx*(1-2.0d0*Vz2)*V_pd_pi
	
	E_z_x2y2 =Vz*(Vx2-Vy2)*(0.5d0*sq3*V_pd_sigma-V_pd_pi)	

	E_z_z2 =Vz*( Vz2-0.5d0*(Vx2+Vy2) )*V_pd_sigma + sq3*Vz*(Vx2+Vy2)*V_pd_pi	

!~ ============================================
end subroutine PD_Bonds_Family
!~ ============================================	




!~ ============================================================
subroutine DD_Bonds_Family( Tswitch, Sx, Sy, Sz, &  
	& E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2, &
	& E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2, &
	& E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2, &
	& E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2, &
	& E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2 )
!~ ============================================================
	use moduleone
	use moduletwo
	implicit real (kind=DP) (a-h, o-z)
	
	integer (kind=4), intent(in) :: Tswitch
	
	real (kind=DP) :: Delta0, Delta1, Delta2, DeltaP, DeltaZ
	real (kind=DP) :: V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi
	real (kind=DP) :: V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp
	real (kind=DP) :: Vn, Vx2, Vy2, Vz2	
	
	real (kind=DP), intent(in) :: Sx, Sy, Sz	

	real (kind=DP), intent(out) :: E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2, &
		& E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2, &
		& E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2, &
		& E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2, &
		& E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2

	Vn=sqrt( Sx**2+Sy**2+Sz**2  )
	Vx=  Sx/Vn  ;       Vy= Sy/Vn  ;      Vz= Sz/Vn     
    Vx2=Vx**2
	Vy2=Vy**2
	Vz2=Vz**2
	
	call Bond_integrals ( Tswitch, Delta0, Delta1, Delta2, DeltaP, DeltaZ, &  
        & V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, &
		& V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp)	
	
	E_xy_xy =3.0d0*Vx2*Vy2*V_dd_sigma + (Vx2+Vy2-4.0*Vx2*Vy2)*V_dd_pi + (Vz2+Vx2*Vy2)*V_dd_delta
	E_xy_yz =3.0d0*Vx*Vy2*Vz*V_dd_sigma + Vx*Vz*(1.0d0-4.0*Vy2)*V_dd_pi + Vx*Vz*(Vy2-1.0d0)*V_dd_delta
	E_xy_xz =3.0d0*Vx2*Vy*Vz*V_dd_sigma + Vy*Vz*(1.0d0-4.0*Vx2)*V_dd_pi + Vy*Vz*(Vx2-1.0d0)*V_dd_delta
	E_xy_x2y2 =1.5d0*Vx*Vy*(Vx2-Vy2)*V_dd_sigma + 2.0d0*Vx*Vy*(Vy2-Vx2)*V_dd_pi + 0.5d0*Vx*Vy*(Vx2-Vy2)*V_dd_delta
	E_xy_z2 =sq3*Vx*Vy*(Vz2-0.5d0*Vx2-0.5d0*Vy2)*V_dd_sigma + 2.0d0*sq3*Vx*Vy*Vz2*V_dd_pi + 0.5d0*sq3*Vx*Vy*(1.0d0+Vz2)*V_dd_delta
	
	
	E_yz_xy =3.0d0*Vx*Vy2*Vz*V_dd_sigma + Vx*Vz*(1.0d0-4.0*Vy2)*V_dd_pi + Vx*Vz*(Vy2-1.0d0)*V_dd_delta
	E_yz_yz =3.0d0*Vy2*Vz2*V_dd_sigma + (Vy2+Vz2-4.0*Vy2*Vz2)*V_dd_pi + (Vx2+Vy2*Vz2)*V_dd_delta
	E_yz_xz =3.0d0*Vy*Vz2*Vx*V_dd_sigma + Vy*Vx*(1.0d0-4.0*Vz2)*V_dd_pi + Vy*Vx*(Vz2-1.0d0)*V_dd_delta
	E_yz_x2y2 =1.5d0*Vy*Vz*(Vx2-Vy2)*V_dd_sigma - Vy*Vz*(1.0d0+2.0d0*(Vx2-Vy2) )*V_dd_pi + Vy*Vz*(1.d0+0.5d0*(Vx2-Vy2) )*V_dd_delta
	E_yz_z2 =sq3*Vy*Vz*(Vz2-0.5d0*Vx2-0.5d0*Vy2)*V_dd_sigma + sq3*Vy*Vz*(Vx2+Vy2-Vz2)*V_dd_pi - 0.5d0*sq3*Vy*Vz*(Vx2+Vy2)*V_dd_delta	
	
	
	E_xz_xy =E_xy_xz
	E_xz_yz =E_yz_xz 
	E_xz_xz =3.0d0*Vz2*Vx2*V_dd_sigma + (Vz2+Vx2-4.0*Vx2*Vz2)*V_dd_pi + (Vy2+Vx2*Vz2)*V_dd_delta
	E_xz_x2y2 =1.5d0*Vx*Vz*(Vx2-Vy2)*V_dd_sigma + Vz*Vx*(1.0d0-2.0d0*(Vx2-Vy2))*V_dd_pi - Vx*Vz*(1.0d0-0.50d0*(Vx2-Vy2))*V_dd_delta
	E_xz_z2 =sq3*Vx*Vz*(Vz2-0.5d0*Vx2-0.5d0*Vy2)*V_dd_sigma + sq3*Vx*Vz*(Vx2+Vy2-Vz2)*V_dd_pi - 0.5d0*sq3*Vx*Vz*(Vx2+Vy2)*V_dd_delta	


	E_x2y2_xy =E_xy_x2y2
	E_x2y2_yz =E_yz_x2y2 
	E_x2y2_xz =E_xz_x2y2
	E_x2y2_x2y2 =3.0d0/4.0d0*(Vx2-Vy2)**2*V_dd_sigma + (Vx2+Vy2-(Vx2-Vy2)**2)*V_dd_pi + (Vz2+0.25d0*(Vx2-Vy2)**2 )*V_dd_delta
	E_x2y2_z2 =0.5d0*sq3*(Vx2-Vy2)*(Vz2-0.50d0*(Vx2+Vy2)) *V_dd_sigma + sq3*Vz2*(Vy2-Vx2)*V_dd_pi + 0.25d0*sq3*(1.0d0+Vz2)*(Vx2-Vy2)*V_dd_delta	

	
	E_z2_xy =E_xy_z2
	E_z2_yz =E_yz_z2 
	E_z2_xz =E_xz_z2
	E_z2_x2y2 =E_x2y2_z2
	E_z2_z2 =  (Vz2-0.50d0*(Vx2+Vy2) )**2*V_dd_sigma + 3.0d0*Vz2*(Vx2+Vy2)*V_dd_pi +0.75d0*(Vx2+Vy2)**2*V_dd_delta	

!~ ============================================
end subroutine DD_Bonds_Family
!~ ============================================	






	
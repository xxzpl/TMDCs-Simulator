!!!!!!!!-----------------------------Subrountine: Initialize the Exchanging Interaction-------one layer---------------------------------------------------------------------------

subroutine Hoppingvalues(Tswitch, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ)
use moduleone
use moduletwo
	
Implicit real (kind=DP) (a-h,o-z) 

integer(KIND=4)::Tswitch

real(kind=DP):: tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)
real(kind=DP)::Delta0,Delta1,Delta2,DeltaP,DeltaZ
real(kind=DP)::V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, V_pp_sigma_interlayer
real(kind=DP)::V_pp_pi_interlayer, SOC_dd,SOC_pp

	tdp(1:3, 1:5, 1:6) = 0.d0
	tdd(1:5, 1:5, 1:6) = 0.d0 
	tpp(1:3, 1:3, 1:10) = 0.d0
	
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
	
!-------------------------------------------------------------------------------------------------------------------------------
	
!	Next Nearest Hopping
!	Mo-Mo  / W-W
!	Z2-1, X2-2, Y2-3
!		6(alpha)	5(-gamma)
!	1(-beta)			4(beta)
!		2(gamma)	3(-alpha)
!	


	tdd(1,1,6)=1./16*(9*V_dd_sigma+4*V_dd_pi+3*V_dd_delta)
	tdd(1,2,6)=0.
	tdd(1,3,6)=0.
	tdd(1,4,6)=sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(1,5,6)=3./8*(V_dd_sigma-V_dd_delta)
	
	tdd(2,1,6)=0.
	tdd(2,2,6)=1./4*(3*V_dd_pi+V_dd_delta)
	tdd(2,3,6)=-sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(2,4,6)=0.
	tdd(2,5,6)=0.
	
	tdd(3,1,6)=0.
	tdd(3,2,6)=-sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(3,3,6)=1./4*(V_dd_pi+ 3*V_dd_delta)
	tdd(3,4,6)=0.
	tdd(3,5,6)=0.	

	tdd(4,1,6)=sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(4,2,6)=0.
	tdd(4,3,6)=0.
	tdd(4,4,6)=1./16*(3*V_dd_sigma+12*V_dd_pi+V_dd_delta)
	tdd(4,5,6)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	
	tdd(5,1,6)=3./8*(V_dd_sigma-V_dd_delta)
	tdd(5,2,6)=0.
	tdd(5,3,6)=0.
	tdd(5,4,6)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	tdd(5,5,6)=1./4*(V_dd_sigma+3*V_dd_delta)

	
!	beta
	tdd(1,1,4)=V_dd_pi
	tdd(1,2,4)=0.
	tdd(1,3,4)=0.
	tdd(1,4,4)=0.
	tdd(1,5,4)=0.
	
	tdd(2,1,4)=0.
	tdd(2,2,4)=V_dd_delta
	tdd(2,3,4)=0.
	tdd(2,4,4)=0.
	tdd(2,5,4)=0.
	
	tdd(3,1,4)=0.
	tdd(3,2,4)=0.
	tdd(3,3,4)=V_dd_pi
	tdd(3,4,4)=0.
	tdd(3,5,4)=0.	

	tdd(4,1,4)=0.
	tdd(4,2,4)=0.
	tdd(4,3,4)=0.
	tdd(4,4,4)=1./4*(3*V_dd_sigma+V_dd_delta)
	tdd(4,5,4)=-sqrt(3.)/4*(V_dd_sigma-V_dd_delta)
	
	tdd(5,1,4)=0.
	tdd(5,2,4)=0.
	tdd(5,3,4)=0.
	tdd(5,4,4)=-sqrt(3.)/4*(V_dd_sigma-V_dd_delta)
	tdd(5,5,4)=1./4*(V_dd_sigma+3*V_dd_delta)
	
	
!	gamma
	tdd(1,1,2)=1./16*(9*V_dd_sigma+4*V_dd_pi+3*V_dd_delta)
	tdd(1,2,2)=0.
	tdd(1,3,2)=0.
	tdd(1,4,2)=-sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(1,5,2)=-3./8*(V_dd_sigma-V_dd_delta)
	
	tdd(2,1,2)=0.
	tdd(2,2,2)=1./4*(3*V_dd_pi+V_dd_delta)
	tdd(2,3,2)=sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(2,4,2)=0.
	tdd(2,5,2)=0.
	
	tdd(3,1,2)=0.
	tdd(3,2,2)=sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(3,3,2)=1./4*(V_dd_pi+ 3*V_dd_delta)
	tdd(3,4,2)=0.
	tdd(3,5,2)=0.	

	tdd(4,1,2)=-sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(4,2,2)=0.
	tdd(4,3,2)=0.
	tdd(4,4,2)=1./16*(3*V_dd_sigma+12*V_dd_pi+V_dd_delta)
	tdd(4,5,2)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	
	tdd(5,1,2)=-3./8*(V_dd_sigma-V_dd_delta)
	tdd(5,2,2)=0.
	tdd(5,3,2)=0.
	tdd(5,4,2)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	tdd(5,5,2)=1./4*(V_dd_sigma+3*V_dd_delta)

	tdd(:,:,1)=tdd(:,:,4)
	tdd(:,:,3)=tdd(:,:,6)
	tdd(:,:,5)=tdd(:,:,2)
	
	
!=============================================================	
	
	
!	Nearest Hoping
!~ W-S
!~ W-Se
!~ M0-S

		    !~ Mo
		 !~ 1(Beta)
              !~ S  Top
!~ 3 (Alpha)      2(Gamma)



	a=sqrt(1./7)/7.d0
!	beta	
	tdp(1,1,1)=14*a*V_pd_pi
	tdp(1,2,1)=0.
	tdp(1,3,1)=-7*a*sqrt(3.)*V_pd_pi
	tdp(1,4,1)=0.
	tdp(1,5,1)=0.
	
	tdp(2,1,1)=0.
	tdp(2,2,1)=-sqrt(3.)*a*(4*sqrt(3.)*V_pd_sigma-V_pd_pi )
	tdp(2,3,1)=0.
	tdp(2,4,1)=-1*a*(4*sqrt(3.)*V_pd_sigma+6*V_pd_pi )
	tdp(2,5,1)=2*a*(V_pd_sigma- 3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,1)=0.
	tdp(3,2,1)=2*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,1)=0.
	tdp(3,4,1)=2*a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,1)=-sqrt(3.)*a*(V_pd_sigma+ 4*sqrt(3.)*V_pd_pi )
	
!	gamma	
	tdp(1,1,2)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,2)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,2)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,2)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,2)=sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,2)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,2)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,2)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,2)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,2)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,2)=3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,2)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,2)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,2)=-a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,2)=-sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )
	
!	alpha
	tdp(1,1,3)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,3)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,3)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,3)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,3)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,3)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,3)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,3)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,3)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,3)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,3)=-3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,3)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,3)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,3)=-a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,3)=-sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )
	

!	Nearest Hoping
!~ W-S
!~ W-Se
!~ M0-S
!	Mo-Se	Bottom layer


		    !~ Mo
		 !~ 4(Beta)
              !~ S  Bottom
!~ 6 (Alpha)      5(Gamma)


	a=sqrt(1./7)/7
!	4beta	
	tdp(1,1,4)=14*a*V_pd_pi
	tdp(1,2,4)=0.
	tdp(1,3,4)=7*a*sqrt(3.)*V_pd_pi
	tdp(1,4,4)=0.
	tdp(1,5,4)=0.
	
	tdp(2,1,4)=0.
	tdp(2,2,4)=sqrt(3.)*a*(4*sqrt(3.)*V_pd_sigma-V_pd_pi )
	tdp(2,3,4)=0.
	tdp(2,4,4)=-1*a*(4*sqrt(3.)*V_pd_sigma+6*V_pd_pi )
	tdp(2,5,4)=2*a*(V_pd_sigma- 3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,4)=0.
	tdp(3,2,4)=2*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,4)=0.
	tdp(3,4,4)=-2*a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,4)=sqrt(3.)*a*(V_pd_sigma+ 4*sqrt(3.)*V_pd_pi )
	
!	gamma	
	tdp(1,1,5)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,5)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,5)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,5)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,5)=sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,5)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,5)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,5)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,5)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,5)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,5)=-3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,5)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,5)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,5)=a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,5)=sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )
	
!	alpha
	tdp(1,1,6)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,6)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,6)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,6)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,6)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,6)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,6)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,6)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,6)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,6)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,6)=3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,6)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,6)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,6)=a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,6)=sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )	



	
!=============================================================
	
	
!	S-S   / Se-Se
!	x-1, y-2, z-3
!		6(alpha)	5(-gamma)
!	1(-beta)			4(beta)
!		2(gamma)	3(-alpha)
!	alpha
	tpp(1,1,6)=1./4*(3*V_pp_pi+V_pp_sigma)
	tpp(1,2,6)=sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(1,3,6)=0.
	tpp(2,1,6)=sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(2,2,6)=1./4*(V_pp_pi+3*V_pp_sigma)
	tpp(2,3,6)=0.
	tpp(3,1,6)=0.
	tpp(3,2,6)=0.
	tpp(3,3,6)=V_pp_pi

!	beta
	tpp(1,1,4)=V_pp_sigma
	tpp(1,2,4)=0.
	tpp(1,3,4)=0.
	tpp(2,1,4)=0.
	tpp(2,2,4)=V_pp_pi
	tpp(2,3,4)=0.
	tpp(3,1,4)=0.
	tpp(3,2,4)=0.
	tpp(3,3,4)=V_pp_pi
	
!	gamma
	tpp(1,1,2)=1./4*(3*V_pp_pi+V_pp_sigma)
	tpp(1,2,2)=-sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(1,3,2)=0.
	tpp(2,1,2)=-sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(2,2,2)=1./4*(V_pp_pi+3*V_pp_sigma)
	tpp(2,3,2)=0.
	tpp(3,1,2)=0.
	tpp(3,2,2)=0.
	tpp(3,3,2)=V_pp_pi 

	tpp(:,:,1)=tpp(:,:,4)
	tpp(:,:,3)=tpp(:,:,6)
	tpp(:,:,5)=tpp(:,:,2)
	

!-------------------------------------- hopping  S(PT)-S(PB) ------------------------------- 	
!	
!	S-S	
!                    S
!     10
!		S
!	
	tpp(1,1,10)=V_pp_pi
	tpp(1,2,10)=0.
	tpp(1,3,10)=0.
	tpp(2,1,10)=0.
	tpp(2,2,10)=V_pp_pi
	tpp(2,3,10)=0.
	tpp(3,1,10)=0.
	tpp(3,2,10)=0.
	tpp(3,3,10)=V_pp_sigma	

!=============================================================


	
!	MoS2
	if (Tswitch==1) then	
!--------------------------------------------------------------------------------------------------------

!                   ----------------------- interlayer hopping  S-S -----------------------------
 
!                        -----------------------	          MoS2         ------------------------------- 	

!-------------------------------------------------------------------------------------------------------- 	
!	Interlayer Hoping
!	S-S	
!    8 (u/ gamma)		9 (v/ alpha)
!		S
!		7 (lambda /beta) 
!	(Mo,S,Direction)    


!	lambda
	tpp(1,1,7)=V_pp_pi_interlayer
	tpp(1,2,7)=0.
	tpp(1,3,7)=0.
	tpp(2,1,7)=0.
	tpp(2,2,7)=0.274234*V_pp_sigma_interlayer+0.725766*V_pp_pi_interlayer
	tpp(2,3,7)=-0.446128*V_pp_sigma_interlayer+0.446128*V_pp_pi_interlayer
	tpp(3,1,7)=0.
	tpp(3,2,7)=-0.446128*V_pp_sigma_interlayer+0.446128*V_pp_pi_interlayer
	tpp(3,3,7)=0.725766*V_pp_sigma_interlayer+0.274234*V_pp_pi_interlayer
	
!	u
	tpp(1,1,8)=0.205676*V_pp_sigma_interlayer+0.794324*V_pp_pi_interlayer
	tpp(1,2,8)=-0.118747*V_pp_sigma_interlayer+ 0.118747*V_pp_pi_interlayer
	tpp(1,3,8)=-0.386358*V_pp_sigma_interlayer+ 0.386358*V_pp_pi_interlayer
	tpp(2,1,8)=-0.118747*V_pp_sigma_interlayer+ 0.118747*V_pp_pi_interlayer
	tpp(2,2,8)=0.0685586*V_pp_sigma_interlayer+0.931441*V_pp_pi_interlayer
	tpp(2,3,8)=0.223064*V_pp_sigma_interlayer- 0.223064*V_pp_pi_interlayer
	tpp(3,1,8)=-0.386358*V_pp_sigma_interlayer+ 0.386358*V_pp_pi_interlayer
	tpp(3,2,8)=0.223064*V_pp_sigma_interlayer- 0.223064*V_pp_pi_interlayer
	tpp(3,3,8)=0.725766*V_pp_sigma_interlayer+0.274234*V_pp_pi_interlayer
	
!	v
	tpp(1,1,9)=0.205676*V_pp_sigma_interlayer+0.794324*V_pp_pi_interlayer
	tpp(1,2,9)=0.118747*V_pp_sigma_interlayer-0.118747*V_pp_pi_interlayer
	tpp(1,3,9)=0.386358*V_pp_sigma_interlayer-0.386358*V_pp_pi_interlayer
	tpp(2,1,9)=0.118747*V_pp_sigma_interlayer-0.118747*V_pp_pi_interlayer
	tpp(2,2,9)=0.0685586*V_pp_sigma_interlayer+0.931441*V_pp_pi_interlayer
	tpp(2,3,9)=0.223064*V_pp_sigma_interlayer- 0.223064*V_pp_pi_interlayer
	tpp(3,1,9)=0.386358*V_pp_sigma_interlayer-0.386358*V_pp_pi_interlayer
	tpp(3,2,9)=0.223064*V_pp_sigma_interlayer- 0.223064*V_pp_pi_interlayer
	tpp(3,3,9)=0.725766*V_pp_sigma_interlayer+0.274234*V_pp_pi_interlayer	
	
!~ !	lambda
	!~ tpp(1,1,7)=V_pp_pi_interlayer
	!~ tpp(1,2,7)=0.
	!~ tpp(1,3,7)=0.
	!~ tpp(2,1,7)=0.
	!~ tpp(2,2,7)=0.214*V_pp_sigma_interlayer+0.786*V_pp_pi_interlayer
	!~ tpp(2,3,7)=0.410*V_pp_sigma_interlayer-0.410*V_pp_pi_interlayer
	!~ tpp(3,1,7)=0.
	!~ tpp(3,2,7)=0.410*V_pp_sigma_interlayer-0.410*V_pp_pi_interlayer
	!~ tpp(3,3,7)=0.785*V_pp_sigma_interlayer+0.215*V_pp_pi_interlayer
	
!~ !	u
	!~ tpp(1,1,8)=0.161*V_pp_sigma_interlayer+0.839*V_pp_pi_interlayer
	!~ tpp(1,2,8)=-0.093*V_pp_sigma_interlayer+ 0.093*V_pp_pi_interlayer
	!~ tpp(1,3,8)=0.355*V_pp_sigma_interlayer- 0.355*V_pp_pi_interlayer
	!~ tpp(2,1,8)=-0.093*V_pp_sigma_interlayer+ 0.093*V_pp_pi_interlayer
	!~ tpp(2,2,8)=0.054*V_pp_sigma_interlayer+0.946*V_pp_pi_interlayer
	!~ tpp(2,3,8)=-0.206*V_pp_sigma_interlayer+ 0.206*V_pp_pi_interlayer
	!~ tpp(3,1,8)=0.355*V_pp_sigma_interlayer- 0.355*V_pp_pi_interlayer
	!~ tpp(3,2,8)=-0.206*V_pp_sigma_interlayer+ 0.206*V_pp_pi_interlayer
	!~ tpp(3,3,8)=0.785*V_pp_sigma_interlayer+0.215*V_pp_pi_interlayer
	
!~ !	v
	!~ tpp(1,1,9)=0.161*V_pp_sigma_interlayer+0.839*V_pp_pi_interlayer
	!~ tpp(1,2,9)=0.093*V_pp_sigma_interlayer-0.093*V_pp_pi_interlayer
	!~ tpp(1,3,9)=-0.355*V_pp_sigma_interlayer+0.355*V_pp_pi_interlayer
	!~ tpp(2,1,9)=0.093*V_pp_sigma_interlayer-0.093*V_pp_pi_interlayer
	!~ tpp(2,2,9)=0.054*V_pp_sigma_interlayer+0.946*V_pp_pi_interlayer
	!~ tpp(2,3,9)=-0.206*V_pp_sigma_interlayer+ 0.206*V_pp_pi_interlayer
	!~ tpp(3,1,9)=-0.355*V_pp_sigma_interlayer+0.355*V_pp_pi_interlayer
	!~ tpp(3,2,9)=-0.206*V_pp_sigma_interlayer+ 0.206*V_pp_pi_interlayer
	!~ tpp(3,3,9)=0.785*V_pp_sigma_interlayer+0.215*V_pp_pi_interlayer	
	
	else if (Tswitch==2) then
!--------------------------------------------------------------------------------------------------------
 
!                   -----------------------	WS2     ------------------------------- 	

!--------------------------------------------------------------------------------------------------------
!	Interlayer Hoping
!	S-S	
!    8 (u/ gamma)		9 (v/ alpha)
!		S
!		7 (lambda /beta) 
!	(Mo,S,Direction)  
	
!	lambda
	tpp(1,1,7)=V_pp_pi_interlayer
	tpp(1,2,7)=0.
	tpp(1,3,7)=0.
	tpp(2,1,7)=0.
	tpp(2,2,7)=0.266766*V_pp_sigma_interlayer+0.733234*V_pp_pi_interlayer
	tpp(2,3,7)=-0.442269*V_pp_sigma_interlayer+0.442269*V_pp_pi_interlayer
	tpp(3,1,7)=0.
	tpp(3,2,7)=-0.442269*V_pp_sigma_interlayer+0.442269*V_pp_pi_interlayer
	tpp(3,3,7)=0.733234*V_pp_sigma_interlayer+0.266766*V_pp_pi_interlayer
	
!	u
	tpp(1,1,8)=0.200075*V_pp_sigma_interlayer+0.799925*V_pp_pi_interlayer
	tpp(1,2,8)=-0.115513*V_pp_sigma_interlayer+0.115513 *V_pp_pi_interlayer
	tpp(1,3,8)=-0.383016*V_pp_sigma_interlayer+ 0.383016*V_pp_pi_interlayer
	tpp(2,1,8)=-0.115513*V_pp_sigma_interlayer+ 0.115513*V_pp_pi_interlayer
	tpp(2,2,8)=0.0666915*V_pp_sigma_interlayer+0.933308*V_pp_pi_interlayer
	tpp(2,3,8)=0.221135*V_pp_sigma_interlayer- 0.221135*V_pp_pi_interlayer
	tpp(3,1,8)=-0.383016*V_pp_sigma_interlayer+0.383016*V_pp_pi_interlayer
	tpp(3,2,8)=0.221135*V_pp_sigma_interlayer- 0.221135*V_pp_pi_interlayer
	tpp(3,3,8)=0.733234*V_pp_sigma_interlayer+0.266766*V_pp_pi_interlayer
	
!	v
	tpp(1,1,9)=0.200075*V_pp_sigma_interlayer+0.799925*V_pp_pi_interlayer
	tpp(1,2,9)=0.115513*V_pp_sigma_interlayer-0.115513 *V_pp_pi_interlayer
	tpp(1,3,9)=0.383016*V_pp_sigma_interlayer-0.383016*V_pp_pi_interlayer
	tpp(2,1,9)=0.115513*V_pp_sigma_interlayer-0.115513*V_pp_pi_interlayer
	tpp(2,2,9)=0.0666915*V_pp_sigma_interlayer+0.933308*V_pp_pi_interlayer
	tpp(2,3,9)=0.221135*V_pp_sigma_interlayer- 0.221135*V_pp_pi_interlayer
	tpp(3,1,9)=0.383016*V_pp_sigma_interlayer-0.383016*V_pp_pi_interlayer
	tpp(3,2,9)=0.221135*V_pp_sigma_interlayer- 0.221135*V_pp_pi_interlayer
	tpp(3,3,9)=0.733234*V_pp_sigma_interlayer+0.266766*V_pp_pi_interlayer

	else if (Tswitch==3) then
!--------------------------------------------------------------------------------------------------------
 
!                   -----------------------	MoSe2     ------------------------------- 	

!--------------------------------------------------------------------------------------------------------
!	Interlayer Hoping
!	S-S	
!    8 (u/ gamma)		9 (v/ alpha)
!		S
!		7 (lambda /beta) 
!	(Mo,S,Direction)  
	
!	lambda
	tpp(1,1,7)=V_pp_pi_interlayer
	tpp(1,2,7)=0.
	tpp(1,3,7)=0.
	tpp(2,1,7)=0.
	tpp(2,2,7)=0.269799*V_pp_sigma_interlayer+0.730201*V_pp_pi_interlayer
	tpp(2,3,7)=-0.443855*V_pp_sigma_interlayer+0.443855*V_pp_pi_interlayer
	tpp(3,1,7)=0.
	tpp(3,2,7)=-0.443855*V_pp_sigma_interlayer+0.443855*V_pp_pi_interlayer
	tpp(3,3,7)=0.730201*V_pp_sigma_interlayer+0.269799*V_pp_pi_interlayer
	
!	u
	tpp(1,1,8)=0.20235*V_pp_sigma_interlayer+0.79765*V_pp_pi_interlayer
	tpp(1,2,8)=-0.116827*V_pp_sigma_interlayer+0.116827*V_pp_pi_interlayer
	tpp(1,3,8)=-0.38439*V_pp_sigma_interlayer+0.38439*V_pp_pi_interlayer
	tpp(2,1,8)=-0.116827*V_pp_sigma_interlayer+0.116827*V_pp_pi_interlayer
	tpp(2,2,8)=0.0674498*V_pp_sigma_interlayer+0.93255*V_pp_pi_interlayer
	tpp(2,3,8)=0.221928*V_pp_sigma_interlayer-0.221928*V_pp_pi_interlayer
	tpp(3,1,8)=-0.38439*V_pp_sigma_interlayer+0.38439*V_pp_pi_interlayer
	tpp(3,2,8)=0.221928*V_pp_sigma_interlayer-0.221928*V_pp_pi_interlayer
	tpp(3,3,8)=0.730201*V_pp_sigma_interlayer+ 0.269799*V_pp_pi_interlayer
	
!	v
	tpp(1,1,9)=0.20235*V_pp_sigma_interlayer+0.79765*V_pp_pi_interlayer
	tpp(1,2,9)=0.116827*V_pp_sigma_interlayer-0.116827*V_pp_pi_interlayer
	tpp(1,3,9)=0.38439*V_pp_sigma_interlayer-0.38439*V_pp_pi_interlayer
	tpp(2,1,9)=0.116827*V_pp_sigma_interlayer-0.116827*V_pp_pi_interlayer
	tpp(2,2,9)=0.0674498*V_pp_sigma_interlayer+0.93255*V_pp_pi_interlayer
	tpp(2,3,9)=0.221928*V_pp_sigma_interlayer-0.221928*V_pp_pi_interlayer
	tpp(3,1,9)=0.38439*V_pp_sigma_interlayer-0.38439*V_pp_pi_interlayer
	tpp(3,2,9)=0.221928*V_pp_sigma_interlayer-0.221928*V_pp_pi_interlayer
	tpp(3,3,9)=0.730201*V_pp_sigma_interlayer+ 0.269799*V_pp_pi_interlayer
	
	else if (Tswitch==4) then
!--------------------------------------------------------------------------------------------------------
 
!                   -----------------------	WSe2     ------------------------------- 	

!--------------------------------------------------------------------------------------------------------
!	Interlayer Hoping
!	S-S	
!    8 (u/ gamma)		9 (v/ alpha)
!		S
!		7 (lambda /beta) 
!	(Mo,S,Direction)  
	
!	lambda
	tpp(1,1,7)=V_pp_pi_interlayer
	tpp(1,2,7)=0.
	tpp(1,3,7)=0.
	tpp(2,1,7)=0.
	tpp(2,2,7)=0.268328*V_pp_sigma_interlayer+ 0.731672*V_pp_pi_interlayer
	tpp(2,3,7)=-0.44309*V_pp_sigma_interlayer+ 0.44309*V_pp_pi_interlayer
	tpp(3,1,7)=0.
	tpp(3,2,7)=-0.44309*V_pp_sigma_interlayer+0.44309*V_pp_pi_interlayer
	tpp(3,3,7)=0.731672*V_pp_sigma_interlayer+ 0.268328*V_pp_pi_interlayer
	
!	u
	tpp(1,1,8)=0.201246*V_pp_sigma_interlayer+0.798754*V_pp_pi_interlayer
	tpp(1,2,8)=-0.11619*V_pp_sigma_interlayer+0.11619*V_pp_pi_interlayer
	tpp(1,3,8)=-0.383727*V_pp_sigma_interlayer+0.383727*V_pp_pi_interlayer
	tpp(2,1,8)=-0.11619*V_pp_sigma_interlayer+0.11619*V_pp_pi_interlayer
	tpp(2,2,8)=0.067082*V_pp_sigma_interlayer+0.932918*V_pp_pi_interlayer
	tpp(2,3,8)=0.221545*V_pp_sigma_interlayer-0.221545*V_pp_pi_interlayer
	tpp(3,1,8)=-0.383727*V_pp_sigma_interlayer+0.383727*V_pp_pi_interlayer
	tpp(3,2,8)=0.221545*V_pp_sigma_interlayer-0.221545*V_pp_pi_interlayer
	tpp(3,3,8)=0.731672*V_pp_sigma_interlayer+0.268328*V_pp_pi_interlayer
	
!	v
	tpp(1,1,9)=0.201246*V_pp_sigma_interlayer+0.798754*V_pp_pi_interlayer
	tpp(1,2,9)=0.11619*V_pp_sigma_interlayer-0.11619*V_pp_pi_interlayer
	tpp(1,3,9)=0.383727*V_pp_sigma_interlayer-0.383727*V_pp_pi_interlayer
	tpp(2,1,9)=0.11619*V_pp_sigma_interlayer-0.11619*V_pp_pi_interlayer
	tpp(2,2,9)=0.067082*V_pp_sigma_interlayer+0.932918*V_pp_pi_interlayer
	tpp(2,3,9)=0.221545*V_pp_sigma_interlayer-0.221545*V_pp_pi_interlayer
	tpp(3,1,9)=0.383727*V_pp_sigma_interlayer-0.383727*V_pp_pi_interlayer
	tpp(3,2,9)=0.221545*V_pp_sigma_interlayer-0.221545*V_pp_pi_interlayer
	tpp(3,3,9)=0.731672*V_pp_sigma_interlayer+0.268328*V_pp_pi_interlayer
		
	
end if


!~ ++++++++++++++++++++++++++++++  check 

	open (unit =UNIT3, file = FILENAME3,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( 8x, a )'  ) 'A_Mo_Mo_nnb_alpha = np.zeros((5, 5), dtype = complex)    ' 	
	write(  UNIT3,  '( 8x, a )'  ) 'A_Mo_Mo_nnb_beta = np.zeros((5, 5), dtype = complex)      '
	write(  UNIT3,  '( 8x, a )'  ) 'A_Mo_Mo_nnb_gamma = np.zeros((5, 5), dtype = complex)	  '
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 5
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdd( i, j, 6 )
		end do
	end do
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 5
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdd( i, j, 4 )
		end do
	end do
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 5
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdd( i, j, 2 )
		end do
	end do

	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( 8x, a )'  ) '# S_top-Mo, nearest neighbours ' 	
	write(  UNIT3,  '( 8x, a )'  ) '# in three directions (alpha, beta, gamma)   '
	write(  UNIT3,  '( 8x, a )'  ) '#	  '	
	write(  UNIT3,  '( 8x, a )'  ) '#      Mo_beta    ' 	
	write(  UNIT3,  '( 8x, a )'  ) '#         S   '
	write(  UNIT3,  '( 8x, a )'  ) '# Mo_alpha  Mo_gamma  '
	write(  UNIT3,  '( 8x, a )'  ) '# HE TMDC XIANG FAN BOTTOM	  '	
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( 8x, a )'  ) 'A_Stop_Mo_nb_alpha = np.zeros((3, 5), dtype = complex)  ' 	
	write(  UNIT3,  '( 8x, a )'  ) 'A_Stop_Mo_nb_beta = np.zeros((3, 5), dtype = complex)      '
	write(  UNIT3,  '( 8x, a )'  ) 'A_Stop_Mo_nb_gamma = np.zeros((3, 5), dtype = complex)	  '	
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 3
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdp( i, j, 1 )
		end do
	end do
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 3
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdp( i, j, 2 )
		end do
	end do
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 3
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdp( i, j, 3 )
		end do
	end do
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( 8x, a )'  ) '# S_bottom-Mo, nearest neighbours ' 	
	write(  UNIT3,  '( 8x, a )'  ) '# in three directions (alpha, beta, gamma)   '
	write(  UNIT3,  '( 8x, a )'  ) '#	  '	
	write(  UNIT3,  '( 8x, a )'  ) '#      Mo_beta    ' 	
	write(  UNIT3,  '( 8x, a )'  ) '#         S   '
	write(  UNIT3,  '( 8x, a )'  ) '# Mo_alpha  Mo_gamma  '
	write(  UNIT3,  '( 8x, a )'  ) '# HE TMDC XIANG FAN BOTTOM	  '	
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 3
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdp( i, j, 4 )
		end do
	end do	
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 3
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdp( i, j, 5 )
		end do
	end do	
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	do  i=1, 3
		do j=1, 5
			write( UNIT3, '( F20.13)' ) tdp( i, j, 6 )
		end do
	end do		
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( 8x, a )'  ) '# S_top-S_top / S_bottom-S_bottom, next-nearest neighbours  ' 	
	write(  UNIT3,  '( 8x, a )'  ) '# in three directions (alpha, beta, gamma)     '
	write(  UNIT3,  '( 8x, a )'  ) '#	  '	
	write(  UNIT3,  '( 8x, a )'  ) '# S_alpha    ' 	
	write(  UNIT3,  '( 8x, a )'  ) '#      S  S_beta    '
	write(  UNIT3,  '( 8x, a )'  ) '# S_gamma	  '
	write(  UNIT3,  '( 8x, a )'  ) '#	  '		
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 	
	do  i=1, 3
		do j=1, 3
			write( UNIT3, '( F20.13)' ) tpp( i, j, 6 )
		end do
	end do		
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 	
	do  i=1, 3
		do j=1, 3
			write( UNIT3, '( F20.13)' ) tpp( i, j, 4 )
		end do
	end do		
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 	
	do  i=1, 3
		do j=1, 3
			write( UNIT3, '( F20.13)' ) tpp( i, j, 2 )
		end do
	end do			

	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( 8x, a )'  ) '# S-S interlayer  ' 	
	write(  UNIT3,  '( 8x, a )'  ) '# in three directions (alpha, beta, gamma)    '
	write(  UNIT3,  '( 8x, a )'  ) '#	  '	
	write(  UNIT3,  '( 8x, a )'  ) '# Stop_gamma   Stop_alpha    ' 	
	write(  UNIT3,  '( 8x, a )'  ) '#             S    '
	write(  UNIT3,  '( 8x, a )'  ) '#         Stop_beta	  '
	write(  UNIT3,  '( 8x, a )'  ) '#	  '	
	write(  UNIT3,  '( a )'  ) 	
	do  i=1, 3
		do j=1, 3
			write( UNIT3, '( F20.13)' ) tpp( i, j, 7 )
		end do
	end do		

	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 	
	do  i=1, 3
		do j=1, 3
			write( UNIT3, '( F20.13)' ) tpp( i, j, 8 )
		end do
	end do		
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 	
	do  i=1, 3
		do j=1, 3
			write( UNIT3, '( F20.13)' ) tpp( i, j, 9 )
		end do
	end do			
	
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 
	write(  UNIT3,  '( a )'  ) 	
	do  i=1, 3
		do j=1, 3
			write( UNIT3, '( F20.13)' ) tpp( i, j, 10 )
		end do
	end do		

!~ ====================================================================

!~ open (unit =UNIT3, file = FILENAME3,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)
!~ write(UNIT3, '( 60("-") )')
!~ write (UNIT3, '("#TDP(1:3,1:5,1:6)", 3(" "), "tdp(1,k,i)", 3(" "), "tdp(2,k,i)", 3(" "), "tdp(3,k,i)" )')
!~ write(UNIT3, '( 60("-") )')
!~ write(UNIT3, '(a)') 
 
   !~ Do i=1,6
          !~ Do k=1,5
           !~ write(UNIT3, '(1x,3(3x2E18.11) )') tdp(1,k,i), tdp(2,k,i), tdp(3,k,i)
	  !~ end do
  !~ end do

!~ write(UNIT3, '(a)')   
!~ write(UNIT3, '( 60("-") )')
!~ write (UNIT3, '("#TPP(1:3,1:3,6|4|2)", 3(" "), "tpp(1,k,i)", 3(" "), "tpp(2,k,i)", 3(" "), "tpp(3,k,i)" )')
!~ write(UNIT3, '( 60("-") )')
!~ write(UNIT3, '(a)') 
 
   !~ Do i=6, 2, -2
          !~ Do k=1,3
           !~ write(UNIT3, '(1x,3(3x2E18.11) )') tpp(1,k,i), tpp(2,k,i), tpp(3,k,i)
	  !~ end do
  !~ end do 
  
!~ write(UNIT3, '(a)')   
!~ write(UNIT3, '( 60("-") )')
!~ write (UNIT3, '("#TDD(1:5,1:5,6|4|2)", 3(" "), "tdd(1,k,i)", 3(" "), "tdd(2,k,i)", 3(" "), "tdd(3,k,i)" &
                         !~ ,3(" "), "tdd(4,k,i)" , 3(" "), "tdd(5,k,i)" )')
!~ write(UNIT3, '( 60("-") )')
!~ write(UNIT3, '(a)') 
!~ ! 
   !~ Do i=6, 2, -2
          !~ Do k=1,5
           !~ write(UNIT3, '(1x,5(3x2E18.11) )') tdd(1,k,i), tdd(2,k,i), tdd(3,k,i), tdd(4,k,i), tdd(5,k,i)
	  !~ end do
  !~ end do 

close(UNIT3)
end 




	
	
	
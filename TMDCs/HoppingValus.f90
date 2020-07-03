!~ ============================================
program InterAtomshopping
!~ ============================================	
	use moduleone
	use moduletwo
	implicit real (kind=DP) (a-h, o-z)
	
	integer (kind=4) :: Tswitch
	
	real (kind=DP) :: Delta0, Delta1, Delta2, DeltaP, DeltaZ
	real (kind=DP) :: V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi
	real (kind=DP) :: V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp
	real (kind=DP) :: Vx, Vy, Vz, Vn, a, Vx2, Vy2, Vz2
	
	real (kind=DP) :: tdp(1:3,1:5,1:6),  tdd(1:5,1:5,1:6),  tpp(1:3,1:3,1:10)

	Tswitch=1
	a1=avalue	
	tdp(1:3, 1:5, 1:6) = 0.d0
	tdd(1:5, 1:5, 1:6) = 0.d0 
	tpp(1:3, 1:3, 1:10) = 0.d0


	
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


	
!=============================================================	
	
	
!	Nearest Hoping
!~ W-S
!~ W-Se
!~ M0-S

		    !~ Mo
		 !~ 1(Beta)
              !~ S  Top
!~ 3 (Alpha)      2(Gamma)



	!~ =============================================	
	!~                                  PD            V-Beta Stop
	!~ =============================================	

	Vx=  0.d0  ;       Vy=   1.d0/sq3*a1   ;      Vz= -0.5d0*a1   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)

	!~ =============================================

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
	
	print*,''
	print*, E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2
	print*,''
	print*, tdp(1, 1:5, 1)

	print*,''	
	print*, E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2
	print*,''
	print*, tdp(2, 1:5, 1)	


	print*,''	
	print*, E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2
	print*,''
	print*, tdp(3, 1:5, 1)	
	

	!~ =============================================	
	!~                                  PD            V-Gamma Stop
	!~ =============================================	

	Vx=  0.5d0*a1  ;       Vy=   -0.5d0/sq3*a1   ;      Vz= -0.5d0*a1   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)

	!~ =============================================

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
	
	print*,''
	print*, E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2
	print*,''
	print*, tdp(1, 1:5, 2)

	print*,''	
	print*, E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2
	print*,''
	print*, tdp(2, 1:5, 2)	


	print*,''	
	print*, E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2
	print*,''
	print*, tdp(3, 1:5, 2)	


	!~ =============================================	
	!~                                  PD            V-alpha Stop
	!~ =============================================	

	Vx=  -0.5d0*a1  ;       Vy=   -0.5d0/sq3*a1   ;      Vz= -0.5d0*a1   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)

	!~ =============================================

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
	
	print*,''
	print*, E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2
	print*,''
	print*, tdp(1, 1:5, 3)

	print*,''	
	print*, E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2
	print*,''
	print*, tdp(2, 1:5, 3)	


	print*,''	
	print*, E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2
	print*,''
	print*, tdp(3, 1:5, 3)	
	

!!~ Bottom


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
	!~ tdp(1,1,6)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	!~ tdp(1,2,6)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	!~ tdp(1,3,6)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	!~ tdp(1,4,6)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	!~ tdp(1,5,6)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	!~ tdp(2,1,6)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	!~ tdp(2,2,6)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	!~ tdp(2,3,6)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	!~ tdp(2,4,6)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	!~ tdp(2,5,6)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	!~ tdp(3,1,6)=3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	!~ tdp(3,2,6)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	!~ tdp(3,3,6)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	!~ tdp(3,4,6)=a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	!~ tdp(3,5,6)=sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )	




	!~ =============================================	
	!~                                  PD            V-Beta S bottom
	!~ =============================================	

	Vx=  0.d0  ;       Vy=   1.d0/sq3*a1   ;      Vz= 0.5d0*a1   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)



	print*,''
	print*, E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2
	print*,''
	print*, tdp(1, 1:5, 4)

	print*,''
	print*,  E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2
	print*,''
	print*, tdp(2, 1:5, 4)
	
	print*,''
	print*, E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2
	print*,''
	print*, tdp(3, 1:5, 4)	


	!~ =============================================	
	!~                                  PD            V-Gamma S bottom
	!~ =============================================	

	Vx=  0.5d0*a1  ;       Vy=   -0.5d0/sq3*a1   ;      Vz= 0.5d0*a1   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)


	print*,''
	print*, E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2
	print*,''
	print*, tdp(1, 1:5, 5)

	print*,''
	print*,  E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2
	print*,''
	print*, tdp(2, 1:5, 5)
	
	print*,''
	print*, E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2
	print*,''
	print*, tdp(3, 1:5, 5)



	!~ =============================================	
	!~                                  PD            V-alpha S bottom
	!~ =============================================	

	Vx=  -0.5d0*a1  ;       Vy=   -0.5d0/sq3*a1   ;      Vz= 0.5d0*a1   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)

	!~ =============================================

    tdp(1, 1:5, 6)= [ E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2  ]
	
    tdp(2, 1:5, 6)= [ E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2  ]

    tdp(3, 1:5, 6)= [ E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2  ]

	print*,''
	print*, E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2
	print*,''
	print*, tdp(1, 1:5, 6)

	print*,''
	print*,  E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2
	print*,''
	print*, tdp(2, 1:5, 6)
	
	print*,''
	print*, E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2
	print*,''
	print*, tdp(3, 1:5, 6)
	
	
	


!~ ============================================
end program InterAtomshopping
!~ ============================================	



	



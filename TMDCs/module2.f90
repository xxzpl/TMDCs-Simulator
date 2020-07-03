	module moduletwo
	use moduleone
	implicit none
	
                                                                                
                                                                                !~ MoS2

	!SOC
	real(kind=DP), parameter :: SOC_dd_MoS2=0.086d0, SOC_pp_MoS2=0.052d0
	!Crystal Fields
	real(kind=DP), parameter :: Delta0_MoS2=-1.094d0, Delta1_MoS2=-0.050d0, Delta2_MoS2=-1.511d0, DeltaP_MoS2=-3.559d0, DeltaZ_MoS2=-6.886d0
	!Mo-S
	real(kind=DP), parameter :: V_pd_sigma_MoS2=3.689d0, V_pd_pi_MoS2=-1.241d0
	!Mo-Mo
	real(kind=DP), parameter :: V_dd_sigma_MoS2=-0.895d0, V_dd_pi_MoS2=0.252d0, V_dd_delta_MoS2=0.228d0
	!S-S
	real(kind=DP), parameter :: V_pp_sigma_MoS2=1.225d0, V_pp_pi_MoS2=-0.467d0
	!S-S Interlayer
	real(kind=DP), parameter :: V_pp_sigma_interlayer_MoS2=-0.774d0, V_pp_pi_interlayer_MoS2=0.123d0
	
	!~ real(kind=DP), parameter :: V_pp_sigma_interlayer_MoS2=-0.55d0, V_pp_pi_interlayer_MoS2=-0.6d0
!------------------------------------------------------------------------------------------------------------------------------------------------------------------	

                                                                                 
																				 !~ WS2

	!SOC
	real(kind=DP), parameter :: SOC_dd_WS2=0.271d0, SOC_pp_WS2=0.057d0
	!Crystal Fields
	real(kind=DP), parameter :: Delta0_WS2=-1.155d0, Delta1_WS2=-0.650d0, Delta2_WS2=-2.279d0, DeltaP_WS2=-3.864d0, DeltaZ_WS2=-7.327d0
	!W-S
	real(kind=DP), parameter :: V_pd_sigma_WS2=4.911d0, V_pd_pi_WS2=-1.220d0
	!W-W
	real(kind=DP), parameter :: V_dd_sigma_WS2=-1.328d0, V_dd_pi_WS2=0.121d0, V_dd_delta_WS2=0.442d0
	!S-S
	real(kind=DP), parameter :: V_pp_sigma_WS2=1.178d0, V_pp_pi_WS2=-0.273d0
	!S-S Interlayer
	real(kind=DP), parameter :: V_pp_sigma_interlayer_WS2=-0.774d0, V_pp_pi_interlayer_WS2=0.123d0
!------------------------------------------------------------------------------------------------------------------------------------------------------------------	



                                                                                 
																				 !~ MoSe2
																				 
	!SOC
	real(kind=DP), parameter :: SOC_dd_MoSe2=0.089d0, SOC_pp_MoSe2=0.256d0
	!Crystal Fields	
	real(kind=DP), parameter :: Delta0_MoSe2=-1.144d0, Delta1_MoSe2=-0.250d0, Delta2_MoSe2=-1.488d0, DeltaP_MoSe2=-4.931d0, DeltaZ_MoSe2=-7.503d0
	!Mo-Se
	real(kind=DP), parameter :: V_pd_sigma_MoSe2=3.728d0, V_pd_pi_MoSe2=-1.222d0
	!Mo-Mo
	real(kind=DP), parameter :: V_dd_sigma_MoSe2=-0.823d0, V_dd_pi_MoSe2=0.215d0, V_dd_delta_MoSe2=0.192d0
	!Se-Se
	real(kind=DP), parameter :: V_pp_sigma_MoSe2=1.256d0, V_pp_pi_MoSe2=-0.205d0

	!Se-Se Interlayer
	real(kind=DP), parameter :: V_pp_sigma_interlayer_MoSe2=-0.774d0, V_pp_pi_interlayer_MoSe2=0.123d0
!------------------------------------------------------------------------------------------------------------------------------------------------------------------	


																				 !~ WSe2

	!SOC
	real(kind=DP), parameter :: SOC_dd_WSe2=0.251d0, SOC_pp_WSe2=0.439d0
	!Crystal Fields
	real(kind=DP), parameter :: Delta0_WSe2=-0.935d0, Delta1_WSe2=-1.250d0, Delta2_WSe2=-2.321d0, DeltaP_WSe2=-5.629d0, DeltaZ_WSe2=-6.759d0
	!W-Se
	real(kind=DP), parameter :: V_pd_sigma_WSe2=5.083d0, V_pd_pi_WSe2=-1.081d0
	!W-W
	real(kind=DP), parameter :: V_dd_sigma_WSe2=-1.129d0, V_dd_pi_WSe2=0.094d0, V_dd_delta_WSe2=0.317d0
	!Se-Se
	real(kind=DP), parameter :: V_pp_sigma_WSe2=1.530d0, V_pp_pi_WSe2=-0.123d0
	!Se-Se Interlayer
	real(kind=DP), parameter :: V_pp_sigma_interlayer_WSe2=-0.55d0, V_pp_pi_interlayer_WSe2=-0.6d0
!------------------------------------------------------------------------------------------------------------------------------------------------------------------		



!~!!  Reference: Electronic band structure of transition metal dichalcogenides from ab initio and slater-kpster tight-binding model
	End module moduletwo
    
    
    

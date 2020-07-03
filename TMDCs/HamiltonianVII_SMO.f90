!~ ==========================================================

        !~ Hamiltonian of XS2    VII version

!~ ==========================================================
subroutine Hamiltonianofxs2VII( neigenvalues,  XK, Yk, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eig, vec) !  ,v1,v2,v3
use moduleone
Implicit real (kind=DP) (a-h,o-z) 
!----------------------------------------------------------------------------------------------------------------------
        !~ integer(kind=4):: Numlayers
        !~ integer(kind=4):: neigenvalues=Numlayers*11
        integer(kind=4):: neigenvalues
        real(kind=DP),dimension(0:neigenvalues-1):: eig
        complex(kind=DP):: vec(0:neigenvalues-1,0:neigenvalues-1)
        complex(kind=DP):: H(0:(neigenvalues*(neigenvalues+1))/2-1)
        complex(kind=DP):: cwork(0:2*neigenvalues-1)
        real(kind=DP):: aux(0:7*neigenvalues-1)
        integer:: iwork(0:5*neigenvalues-1),ifail(0:neigenvalues-1)
        real(kind=DP),parameter:: abstol=2*tiny(abstol)  ! see ZHPEVX

        !!! =========  Define more parameters     =====
        integer(kind=4):: NdimensionVonset,  lneigenvalues
        real(kind=DP):: XK, Yk
        real(kind=DP):: tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)
        real(kind=DP)::Delta0,Delta1,Delta2,DeltaP,DeltaZ
        real(kind=DP):: v12, v23
        real(kind=DP), dimension ( 1:3, 1:(neigenvalues/11) ):: Vonset    !~! Allocate(  Vonset  ( 1:3, Nlayer  )     )
    
        complex(kind=DP), dimension(neigenvalues,neigenvalues)::Hamiltonian, Hamiltonian0, Hamiltonian1   
        
        complex(kind=DP), dimension(9:11, 12:14)::  Hamiltonian1inter,  Hamiltonian0inter
        
        complex(kind=DP)::Expon_KRalpha, Expon_KRbeta,Expon_KRgemma,  phase 		
    
!======================================================================================

!~ CALL MKL_SET_NUM_THREADS(1)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      Hamiltonian0
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xRalpha   =   -Vec2X        ;      yRalpha    =   -Vec2Y
xRbeta    =    Vec1X+Vec2X  ;      yRbeta     =    Vec1Y+Vec2Y
xRgemma   =   -Vec1X        ;      yRgemma    =   -Vec1Y

CS_KRalpha      =   COS(xk*xRalpha+yk*yRalpha)
CS_KRbeta       =   COS(xk*xRbeta+yk*yRbeta)
CS_KRgemma      =   COS(xk*xRgemma+yk*yRgemma)

!!%%%%%%%%%%%%%%        ST-----ST            %%%%%%%%%%%%%%%%

Hamiltonian0(1,1)=2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP  !+v23
Hamiltonian0(1,2)=2.0_DP*( tpp(1,2,6)*CS_KRalpha+tpp(1,2,4)*CS_KRbeta+tpp(1,2,2)*CS_KRgemma )
Hamiltonian0(1,3)=2.0_DP*( tpp(1,3,6)*CS_KRalpha+tpp(1,3,4)*CS_KRbeta+tpp(1,3,2)*CS_KRgemma )

Hamiltonian0(2,2)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP   !+v23
Hamiltonian0(2,3)=2.0_DP*( tpp(2,3,6)*CS_KRalpha+tpp(2,3,4)*CS_KRbeta+tpp(2,3,2)*CS_KRgemma )
Hamiltonian0(3,3)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ   !+v23

!!%%%%%%%%%%%%%%        SB-----SB            %%%%%%%%%%%%%%%%

Hamiltonian0(9,9) =2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP   !-v12
Hamiltonian0(9,10)=Hamiltonian0(1,2)
Hamiltonian0(9,11)=Hamiltonian0(1,3)

Hamiltonian0(10,10)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP  !-v12
Hamiltonian0(10,11)=Hamiltonian0(2,3)
Hamiltonian0(11,11)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ  !-v12

!!%%%%%%%%%%%%%%        Mo-----Mo            %%%%%%%%%%%%%%%%

Hamiltonian0(4,4)=2.0_DP*( tdd(1,1,6)*CS_KRalpha+tdd(1,1,4)*CS_KRbeta+tdd(1,1,2)*CS_KRgemma )+Delta2  !  v2=0 eV
Hamiltonian0(4,5)=2.0_DP*( tdd(1,2,6)*CS_KRalpha+tdd(1,2,4)*CS_KRbeta+tdd(1,2,2)*CS_KRgemma )
Hamiltonian0(4,6)=2.0_DP*( tdd(1,3,6)*CS_KRalpha+tdd(1,3,4)*CS_KRbeta+tdd(1,3,2)*CS_KRgemma )
Hamiltonian0(4,7)=2.0_DP*( tdd(1,4,6)*CS_KRalpha+tdd(1,4,4)*CS_KRbeta+tdd(1,4,2)*CS_KRgemma )
Hamiltonian0(4,8)=2.0_DP*( tdd(1,5,6)*CS_KRalpha+tdd(1,5,4)*CS_KRbeta+tdd(1,5,2)*CS_KRgemma )

Hamiltonian0(5,5)=2.0_DP*( tdd(2,2,6)*CS_KRalpha+tdd(2,2,4)*CS_KRbeta+tdd(2,2,2)*CS_KRgemma )+Delta1  !  v2=0 eV
Hamiltonian0(5,6)=2.0_DP*( tdd(2,3,6)*CS_KRalpha+tdd(2,3,4)*CS_KRbeta+tdd(2,3,2)*CS_KRgemma )
Hamiltonian0(5,7)=2.0_DP*( tdd(2,4,6)*CS_KRalpha+tdd(2,4,4)*CS_KRbeta+tdd(2,4,2)*CS_KRgemma )
Hamiltonian0(5,8)=2.0_DP*( tdd(2,5,6)*CS_KRalpha+tdd(2,5,4)*CS_KRbeta+tdd(2,5,2)*CS_KRgemma )

Hamiltonian0(6,6)=2.0_DP*( tdd(3,3,6)*CS_KRalpha+tdd(3,3,4)*CS_KRbeta+tdd(3,3,2)*CS_KRgemma )+Delta1    !  v2=0 eV
Hamiltonian0(6,7)=2.0_DP*( tdd(3,4,6)*CS_KRalpha+tdd(3,4,4)*CS_KRbeta+tdd(3,4,2)*CS_KRgemma )
Hamiltonian0(6,8)=2.0_DP*( tdd(3,5,6)*CS_KRalpha+tdd(3,5,4)*CS_KRbeta+tdd(3,5,2)*CS_KRgemma )

Hamiltonian0(7,7)=2.0_DP*( tdd(4,4,6)*CS_KRalpha+tdd(4,4,4)*CS_KRbeta+tdd(4,4,2)*CS_KRgemma )+Delta2    !  v2=0 eV
Hamiltonian0(7,8)=2.0_DP*( tdd(4,5,6)*CS_KRalpha+tdd(4,5,4)*CS_KRbeta+tdd(4,5,2)*CS_KRgemma )
Hamiltonian0(8,8)=2.0_DP*( tdd(5,5,6)*CS_KRalpha+tdd(5,5,4)*CS_KRbeta+tdd(5,5,2)*CS_KRgemma )+Delta0    !  v2=0 eV

!!%%%%%%%%%%%%%%        Mo-----> ST    dega+        %%%%%%%%%%%%%%%%

xRalpha   =  -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
xRbeta    =   0                                        ; yRbeta    =    carbondistance    
xRgemma   =   0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma   =   -0.5*carbondistance 

Expon_KRalpha      =exp(  -AI*(xk*xRalpha+yk*yRalpha )  )
Expon_KRbeta       =exp(  -AI*(xk*xRbeta+yk*yRbeta )    )
Expon_KRgemma      =exp(  -AI*(xk*xRgemma+yk*yRgemma)   )

phase =exp(AI*(xk*xRbeta+yk*yRbeta ))     !!!!   Bloch theory in one Unit Cell

Hamiltonian0(1,4)=tdp(1,1,3)*Expon_KRalpha+tdp(1,1,2)*Expon_KRgemma+tdp(1,1,1)*Expon_KRbeta
Hamiltonian0(1,5)=tdp(1,2,3)*Expon_KRalpha+tdp(1,2,2)*Expon_KRgemma+tdp(1,2,1)*Expon_KRbeta
Hamiltonian0(1,6)=tdp(1,3,3)*Expon_KRalpha+tdp(1,3,2)*Expon_KRgemma+tdp(1,3,1)*Expon_KRbeta 
Hamiltonian0(1,7)=tdp(1,4,3)*Expon_KRalpha+tdp(1,4,2)*Expon_KRgemma+tdp(1,4,1)*Expon_KRbeta  
Hamiltonian0(1,8)=tdp(1,5,3)*Expon_KRalpha+tdp(1,5,2)*Expon_KRgemma+tdp(1,5,1)*Expon_KRbeta 

Hamiltonian0(2,4)=tdp(2,1,3)*Expon_KRalpha+tdp(2,1,2)*Expon_KRgemma+tdp(2,1,1)*Expon_KRbeta
Hamiltonian0(2,5)=tdp(2,2,3)*Expon_KRalpha+tdp(2,2,2)*Expon_KRgemma+tdp(2,2,1)*Expon_KRbeta  
Hamiltonian0(2,6)=tdp(2,3,3)*Expon_KRalpha+tdp(2,3,2)*Expon_KRgemma+tdp(2,3,1)*Expon_KRbeta
Hamiltonian0(2,7)=tdp(2,4,3)*Expon_KRalpha+tdp(2,4,2)*Expon_KRgemma+tdp(2,4,1)*Expon_KRbeta 
Hamiltonian0(2,8)=tdp(2,5,3)*Expon_KRalpha+tdp(2,5,2)*Expon_KRgemma+tdp(2,5,1)*Expon_KRbeta 

Hamiltonian0(3,4)=tdp(3,1,3)*Expon_KRalpha+tdp(3,1,2)*Expon_KRgemma+tdp(3,1,1)*Expon_KRbeta  
Hamiltonian0(3,5)=tdp(3,2,3)*Expon_KRalpha+tdp(3,2,2)*Expon_KRgemma+tdp(3,2,1)*Expon_KRbeta 
Hamiltonian0(3,6)=tdp(3,3,3)*Expon_KRalpha+tdp(3,3,2)*Expon_KRgemma+tdp(3,3,1)*Expon_KRbeta 
Hamiltonian0(3,7)=tdp(3,4,3)*Expon_KRalpha+tdp(3,4,2)*Expon_KRgemma+tdp(3,4,1)*Expon_KRbeta  
Hamiltonian0(3,8)=tdp(3,5,3)*Expon_KRalpha+tdp(3,5,2)*Expon_KRgemma+tdp(3,5,1)*Expon_KRbeta
Hamiltonian0(1:3 , 4:8 )=Hamiltonian0(1:3 , 4:8 )*phase

!!%%%%%%%%%%%%%%        SB -----> Mo            %%%%%%%%%%%%%%%%

xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
xRbeta    =  0.d0                                     ; yRbeta    =    carbondistance   
xRgemma   =  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma   =   -0.5*carbondistance 

Expon_KRalpha      =exp(AI*(xk*xRalpha+yk*yRalpha ))
Expon_KRbeta       =exp(AI*(xk*xRbeta+yk*yRbeta))
Expon_KRgemma      =exp(AI*(xk*xRgemma+yk*yRgemma))

phase =exp(-AI*(xk*xRbeta+yk*yRbeta ))    !!!!   Bloch theory in one Unit Cell

Hamiltonian0(4,9)  =tdp(1,1,6)*Expon_KRalpha+tdp(1,1,5)*Expon_KRgemma+tdp(1,1,4)*Expon_KRbeta
Hamiltonian0(5,9)  =tdp(1,2,6)*Expon_KRalpha+tdp(1,2,5)*Expon_KRgemma+tdp(1,2,4)*Expon_KRbeta 
Hamiltonian0(6,9)  =tdp(1,3,6)*Expon_KRalpha+tdp(1,3,5)*Expon_KRgemma+tdp(1,3,4)*Expon_KRbeta 
Hamiltonian0(7,9)  =tdp(1,4,6)*Expon_KRalpha+tdp(1,4,5)*Expon_KRgemma+tdp(1,4,4)*Expon_KRbeta 
Hamiltonian0(8,9)  =tdp(1,5,6)*Expon_KRalpha+tdp(1,5,5)*Expon_KRgemma+tdp(1,5,4)*Expon_KRbeta 

Hamiltonian0(4,10)  =tdp(2,1,6)*Expon_KRalpha+tdp(2,1,5)*Expon_KRgemma+tdp(2,1,4)*Expon_KRbeta 
Hamiltonian0(5,10)  =tdp(2,2,6)*Expon_KRalpha+tdp(2,2,5)*Expon_KRgemma+tdp(2,2,4)*Expon_KRbeta 
Hamiltonian0(6,10)  =tdp(2,3,6)*Expon_KRalpha+tdp(2,3,5)*Expon_KRgemma+tdp(2,3,4)*Expon_KRbeta 
Hamiltonian0(7,10)  =tdp(2,4,6)*Expon_KRalpha+tdp(2,4,5)*Expon_KRgemma+tdp(2,4,4)*Expon_KRbeta 
Hamiltonian0(8,10)  =tdp(2,5,6)*Expon_KRalpha+tdp(2,5,5)*Expon_KRgemma+tdp(2,5,4)*Expon_KRbeta 

Hamiltonian0(4,11)  =tdp(3,1,6)*Expon_KRalpha+tdp(3,1,5)*Expon_KRgemma+tdp(3,1,4)*Expon_KRbeta 
Hamiltonian0(5,11)  =tdp(3,2,6)*Expon_KRalpha+tdp(3,2,5)*Expon_KRgemma+tdp(3,2,4)*Expon_KRbeta 
Hamiltonian0(6,11)  =tdp(3,3,6)*Expon_KRalpha+tdp(3,3,5)*Expon_KRgemma+tdp(3,3,4)*Expon_KRbeta 
Hamiltonian0(7,11)  =tdp(3,4,6)*Expon_KRalpha+tdp(3,4,5)*Expon_KRgemma+tdp(3,4,4)*Expon_KRbeta 
Hamiltonian0(8,11)  =tdp(3,5,6)*Expon_KRalpha+tdp(3,5,5)*Expon_KRgemma+tdp(3,5,4)*Expon_KRbeta 
Hamiltonian0(4:8 , 9:11)  =Hamiltonian0(4:8 , 9:11)*phase 

!!%%%%%%%%%%%%%%        SB -----> ST           %%%%%%%%%%%%%%%%
!!    SB----> ST

    Hamiltonian0(1,9) =tpp(1,1,10)
    Hamiltonian0(1,10)=tpp(1,2,10)
    Hamiltonian0(1,11)=tpp(1,3,10)

    Hamiltonian0(2,9) =tpp(2,1,10)
    Hamiltonian0(2,10)=tpp(2,2,10)
    Hamiltonian0(2,11)=tpp(2,3,10)

    Hamiltonian0(3,9) =tpp(3,1,10)
    Hamiltonian0(3,10)=tpp(3,2,10)
    Hamiltonian0(3,11)=tpp(3,3,10)

!!%%%%%%%%%%%%%%        ST one layer -----> SB up layer    daga+       %%%%%%%%%%%%%%%%
 
    xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
    xRbeta    =  0.d0                                     ; yRbeta    =    carbondistance          
    xRgemma   =  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma   =   -0.5*carbondistance
 
    Expon_KRalpha      =exp(AI*(xk*xRalpha+yk*yRalpha ))
    Expon_KRbeta       =exp(AI*(xk*xRbeta+yk*yRbeta ))
    Expon_KRgemma      =exp(AI*(xk*xRgemma+yk*yRgemma))
    
    phase =exp(-AI*(xk*xRbeta+yk*yRbeta ))    !!!!   Bloch theory in one Unit Cell
    
    Hamiltonian0inter(9,12)=tpp(1,1,9)*Expon_KRalpha+tpp(1,1,8)*Expon_KRgemma+tpp(1,1,7)*Expon_KRbeta
    Hamiltonian0inter(9,13)=tpp(1,2,9)*Expon_KRalpha+tpp(1,2,8)*Expon_KRgemma+tpp(1,2,7)*Expon_KRbeta
    Hamiltonian0inter(9,14)=tpp(1,3,9)*Expon_KRalpha+tpp(1,3,8)*Expon_KRgemma+tpp(1,3,7)*Expon_KRbeta	

    Hamiltonian0inter(10,12)=tpp(2,1,9)*Expon_KRalpha+tpp(2,1,8)*Expon_KRgemma+tpp(2,1,7)*Expon_KRbeta
    Hamiltonian0inter(10,13)=tpp(2,2,9)*Expon_KRalpha+tpp(2,2,8)*Expon_KRgemma+tpp(2,2,7)*Expon_KRbeta
    Hamiltonian0inter(10,14)=tpp(2,3,9)*Expon_KRalpha+tpp(2,3,8)*Expon_KRgemma+tpp(2,3,7)*Expon_KRbeta	

    Hamiltonian0inter(11,12)=tpp(3,1,9)*Expon_KRalpha+tpp(3,1,8)*Expon_KRgemma+tpp(3,1,7)*Expon_KRbeta
    Hamiltonian0inter(11,13)=tpp(3,2,9)*Expon_KRalpha+tpp(3,2,8)*Expon_KRgemma+tpp(3,2,7)*Expon_KRbeta
    Hamiltonian0inter(11,14)=tpp(3,3,9)*Expon_KRalpha+tpp(3,3,8)*Expon_KRgemma+tpp(3,3,7)*Expon_KRbeta
    Hamiltonian0inter( 9:11, 12:14)=phase*Hamiltonian0inter( 9:11, 12:14)

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         Hamiltonian1
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
!======================================================================================
    xRalpha   =   -Vec2X            ;       yRalpha    =   -Vec2Y
    xRbeta    =    Vec1X+Vec2X      ;       yRbeta     =    Vec1Y+Vec2Y
    xRgemma   =   -Vec1X            ;       yRgemma    =   -Vec1Y

    xRalpha = -xRalpha  ;   yRalpha  = -yRalpha
    xRbeta  = -xRbeta   ;   yRbeta   = -yRbeta
    xRgemma = -xRgemma  ;   yRgemma  = -yRgemma

    CS_KRalpha      =   COS(xk*xRalpha+yk*yRalpha)
    CS_KRbeta       =   COS(xk*xRbeta+yk*yRbeta)
    CS_KRgemma      =   COS(xk*xRgemma+yk*yRgemma)

!!%%%%%%%%%%%%%%        ST-----ST            %%%%%%%%%%%%%%%%

    Hamiltonian1(1,1)=2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP !+v23
    Hamiltonian1(1,2)=2.0_DP*( tpp(1,2,6)*CS_KRalpha+tpp(1,2,4)*CS_KRbeta+tpp(1,2,2)*CS_KRgemma )
    Hamiltonian1(1,3)=2.0_DP*( tpp(1,3,6)*CS_KRalpha+tpp(1,3,4)*CS_KRbeta+tpp(1,3,2)*CS_KRgemma )

    Hamiltonian1(2,2)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP !+v23
    Hamiltonian1(2,3)=2.0_DP*( tpp(2,3,6)*CS_KRalpha+tpp(2,3,4)*CS_KRbeta+tpp(2,3,2)*CS_KRgemma )
    Hamiltonian1(3,3)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ !+v23

!!%%%%%%%%%%%%%%        SB-----SB            %%%%%%%%%%%%%%%%

    Hamiltonian1(9,9)  =2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP !-v12
    Hamiltonian1(9,10)=Hamiltonian1(1,2)
    Hamiltonian1(9,11)=Hamiltonian1(1,3)

    Hamiltonian1(10,10)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP !-v12
    Hamiltonian1(10,11)=Hamiltonian1(2,3)
    Hamiltonian1(11,11)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ !-v12

!!%%%%%%%%%%%%%%        Mo-----Mo            %%%%%%%%%%%%%%%%

    Hamiltonian1(4,4)=2.0_DP*( tdd(1,1,6)*CS_KRalpha+tdd(1,1,4)*CS_KRbeta+tdd(1,1,2)*CS_KRgemma )+Delta2  !  v2=0 eV
    Hamiltonian1(4,5)=2.0_DP*( tdd(1,2,6)*CS_KRalpha+tdd(1,2,4)*CS_KRbeta+tdd(1,2,2)*CS_KRgemma )
    Hamiltonian1(4,6)=2.0_DP*( tdd(1,3,6)*CS_KRalpha+tdd(1,3,4)*CS_KRbeta+tdd(1,3,2)*CS_KRgemma )
    Hamiltonian1(4,7)=2.0_DP*( tdd(1,4,6)*CS_KRalpha+tdd(1,4,4)*CS_KRbeta+tdd(1,4,2)*CS_KRgemma )
    Hamiltonian1(4,8)=2.0_DP*( tdd(1,5,6)*CS_KRalpha+tdd(1,5,4)*CS_KRbeta+tdd(1,5,2)*CS_KRgemma )

    Hamiltonian1(5,5)=2.0_DP*( tdd(2,2,6)*CS_KRalpha+tdd(2,2,4)*CS_KRbeta+tdd(2,2,2)*CS_KRgemma )+Delta1  !  v2=0 eV
    Hamiltonian1(5,6)=2.0_DP*( tdd(2,3,6)*CS_KRalpha+tdd(2,3,4)*CS_KRbeta+tdd(2,3,2)*CS_KRgemma )
    Hamiltonian1(5,7)=2.0_DP*( tdd(2,4,6)*CS_KRalpha+tdd(2,4,4)*CS_KRbeta+tdd(2,4,2)*CS_KRgemma )
    Hamiltonian1(5,8)=2.0_DP*( tdd(2,5,6)*CS_KRalpha+tdd(2,5,4)*CS_KRbeta+tdd(2,5,2)*CS_KRgemma )

    Hamiltonian1(6,6)=2.0_DP*( tdd(3,3,6)*CS_KRalpha+tdd(3,3,4)*CS_KRbeta+tdd(3,3,2)*CS_KRgemma )+Delta1    !  v2=0 eV
    Hamiltonian1(6,7)=2.0_DP*( tdd(3,4,6)*CS_KRalpha+tdd(3,4,4)*CS_KRbeta+tdd(3,4,2)*CS_KRgemma )
    Hamiltonian1(6,8)=2.0_DP*( tdd(3,5,6)*CS_KRalpha+tdd(3,5,4)*CS_KRbeta+tdd(3,5,2)*CS_KRgemma )

    Hamiltonian1(7,7)=2.0_DP*( tdd(4,4,6)*CS_KRalpha+tdd(4,4,4)*CS_KRbeta+tdd(4,4,2)*CS_KRgemma )+Delta2    !  v2=0 eV
    Hamiltonian1(7,8)=2.0_DP*( tdd(4,5,6)*CS_KRalpha+tdd(4,5,4)*CS_KRbeta+tdd(4,5,2)*CS_KRgemma )
    Hamiltonian1(8,8)=2.0_DP*( tdd(5,5,6)*CS_KRalpha+tdd(5,5,4)*CS_KRbeta+tdd(5,5,2)*CS_KRgemma )+Delta0    !  v2=0 eV

!!%%%%%%%%%%%%%%        Mo-----> ST    dega+        %%%%%%%%%%%%%%%%

    xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
    xRbeta    =   0                                       ; yRbeta    =    carbondistance          
    xRgemma   =  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma   =   -0.5*carbondistance 

    xRalpha = -xRalpha  ;   yRalpha  = -yRalpha
    xRbeta  = -xRbeta   ;   yRbeta   = -yRbeta
    xRgemma = -xRgemma  ;   yRgemma  = -yRgemma

    Expon_KRalpha   = exp(-AI*(xk*xRalpha+yk*yRalpha ))
    Expon_KRbeta    = exp(-AI*(xk*xRbeta+yk*yRbeta ))
    Expon_KRgemma   = exp(-AI*(xk*xRgemma+yk*yRgemma))

    phase = exp(AI*(xk*xRbeta+yk*yRbeta ))     !!!!   Bloch theory in one Unit Cell

    Hamiltonian1(1,4)=tdp(1,1,3)*Expon_KRalpha+tdp(1,1,2)*Expon_KRgemma+tdp(1,1,1)*Expon_KRbeta
    Hamiltonian1(1,5)=tdp(1,2,3)*Expon_KRalpha+tdp(1,2,2)*Expon_KRgemma+tdp(1,2,1)*Expon_KRbeta
    Hamiltonian1(1,6)=tdp(1,3,3)*Expon_KRalpha+tdp(1,3,2)*Expon_KRgemma+tdp(1,3,1)*Expon_KRbeta 
    Hamiltonian1(1,7)=tdp(1,4,3)*Expon_KRalpha+tdp(1,4,2)*Expon_KRgemma+tdp(1,4,1)*Expon_KRbeta  
    Hamiltonian1(1,8)=tdp(1,5,3)*Expon_KRalpha+tdp(1,5,2)*Expon_KRgemma+tdp(1,5,1)*Expon_KRbeta 

    Hamiltonian1(2,4)=tdp(2,1,3)*Expon_KRalpha+tdp(2,1,2)*Expon_KRgemma+tdp(2,1,1)*Expon_KRbeta
    Hamiltonian1(2,5)=tdp(2,2,3)*Expon_KRalpha+tdp(2,2,2)*Expon_KRgemma+tdp(2,2,1)*Expon_KRbeta  
    Hamiltonian1(2,6)=tdp(2,3,3)*Expon_KRalpha+tdp(2,3,2)*Expon_KRgemma+tdp(2,3,1)*Expon_KRbeta
    Hamiltonian1(2,7)=tdp(2,4,3)*Expon_KRalpha+tdp(2,4,2)*Expon_KRgemma+tdp(2,4,1)*Expon_KRbeta 
    Hamiltonian1(2,8)=tdp(2,5,3)*Expon_KRalpha+tdp(2,5,2)*Expon_KRgemma+tdp(2,5,1)*Expon_KRbeta 

    Hamiltonian1(3,4)=tdp(3,1,3)*Expon_KRalpha+tdp(3,1,2)*Expon_KRgemma+tdp(3,1,1)*Expon_KRbeta  
    Hamiltonian1(3,5)=tdp(3,2,3)*Expon_KRalpha+tdp(3,2,2)*Expon_KRgemma+tdp(3,2,1)*Expon_KRbeta 
    Hamiltonian1(3,6)=tdp(3,3,3)*Expon_KRalpha+tdp(3,3,2)*Expon_KRgemma+tdp(3,3,1)*Expon_KRbeta 
    Hamiltonian1(3,7)=tdp(3,4,3)*Expon_KRalpha+tdp(3,4,2)*Expon_KRgemma+tdp(3,4,1)*Expon_KRbeta  
    Hamiltonian1(3,8)=tdp(3,5,3)*Expon_KRalpha+tdp(3,5,2)*Expon_KRgemma+tdp(3,5,1)*Expon_KRbeta
    Hamiltonian1(1:3 , 4:8 )=Hamiltonian1(1:3 , 4:8 )*phase

!!%%%%%%%%%%%%%%        SB -----> Mo            %%%%%%%%%%%%%%%%

    xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
    xRbeta    =   0                                       ; yRbeta    =    carbondistance          
    xRgemma   =  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma   =   -0.5*carbondistance 

    xRalpha = -xRalpha  ;   yRalpha  = -yRalpha
    xRbeta  = -xRbeta   ;   yRbeta   = -yRbeta
    xRgemma = -xRgemma  ;   yRgemma  = -yRgemma

    Expon_KRalpha   = exp(AI*(xk*xRalpha+yk*yRalpha ))
    Expon_KRbeta    = exp(AI*(xk*xRbeta+yk*yRbeta))
    Expon_KRgemma   = exp(AI*(xk*xRgemma+yk*yRgemma))

    phase = exp(-AI*(xk*xRbeta+yk*yRbeta ))    !!!!   Bloch theory in one Unit Cell

    Hamiltonian1(4,9)  =tdp(1,1,6)*Expon_KRalpha+tdp(1,1,5)*Expon_KRgemma+tdp(1,1,4)*Expon_KRbeta
    Hamiltonian1(5,9)  =tdp(1,2,6)*Expon_KRalpha+tdp(1,2,5)*Expon_KRgemma+tdp(1,2,4)*Expon_KRbeta 
    Hamiltonian1(6,9)  =tdp(1,3,6)*Expon_KRalpha+tdp(1,3,5)*Expon_KRgemma+tdp(1,3,4)*Expon_KRbeta 
    Hamiltonian1(7,9)  =tdp(1,4,6)*Expon_KRalpha+tdp(1,4,5)*Expon_KRgemma+tdp(1,4,4)*Expon_KRbeta 
    Hamiltonian1(8,9)  =tdp(1,5,6)*Expon_KRalpha+tdp(1,5,5)*Expon_KRgemma+tdp(1,5,4)*Expon_KRbeta 

    Hamiltonian1(4,10)  =tdp(2,1,6)*Expon_KRalpha+tdp(2,1,5)*Expon_KRgemma+tdp(2,1,4)*Expon_KRbeta 
    Hamiltonian1(5,10)  =tdp(2,2,6)*Expon_KRalpha+tdp(2,2,5)*Expon_KRgemma+tdp(2,2,4)*Expon_KRbeta 
    Hamiltonian1(6,10)  =tdp(2,3,6)*Expon_KRalpha+tdp(2,3,5)*Expon_KRgemma+tdp(2,3,4)*Expon_KRbeta 
    Hamiltonian1(7,10)  =tdp(2,4,6)*Expon_KRalpha+tdp(2,4,5)*Expon_KRgemma+tdp(2,4,4)*Expon_KRbeta 
    Hamiltonian1(8,10)  =tdp(2,5,6)*Expon_KRalpha+tdp(2,5,5)*Expon_KRgemma+tdp(2,5,4)*Expon_KRbeta 

    Hamiltonian1(4,11)  =tdp(3,1,6)*Expon_KRalpha+tdp(3,1,5)*Expon_KRgemma+tdp(3,1,4)*Expon_KRbeta 
    Hamiltonian1(5,11)  =tdp(3,2,6)*Expon_KRalpha+tdp(3,2,5)*Expon_KRgemma+tdp(3,2,4)*Expon_KRbeta 
    Hamiltonian1(6,11)  =tdp(3,3,6)*Expon_KRalpha+tdp(3,3,5)*Expon_KRgemma+tdp(3,3,4)*Expon_KRbeta 
    Hamiltonian1(7,11)  =tdp(3,4,6)*Expon_KRalpha+tdp(3,4,5)*Expon_KRgemma+tdp(3,4,4)*Expon_KRbeta 
    Hamiltonian1(8,11)  =tdp(3,5,6)*Expon_KRalpha+tdp(3,5,5)*Expon_KRgemma+tdp(3,5,4)*Expon_KRbeta 
    Hamiltonian1(4:8 , 9:11)  =Hamiltonian1(4:8 , 9:11)*phase 

!!%%%%%%%%%%%%%%        SB -----> ST           %%%%%%%%%%%%%%%%
!!    SB----> ST

    Hamiltonian1(1,9)=tpp(1,1,10)
    Hamiltonian1(1,10)=tpp(1,2,10)
    Hamiltonian1(1,11)=tpp(1,3,10)

    Hamiltonian1(2,9)=tpp(2,1,10)
    Hamiltonian1(2,10)=tpp(2,2,10)
    Hamiltonian1(2,11)=tpp(2,3,10)

    Hamiltonian1(3,9)=tpp(3,1,10)
    Hamiltonian1(3,10)=tpp(3,2,10)
    Hamiltonian1(3,11)=tpp(3,3,10)

!!%%%%%%%%%%%%%%        ST one layer -----> SB up layer    daga+       %%%%%%%%%%%%%%%%

    xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
    xRbeta    =  0.d0                                     ; yRbeta    =    carbondistance          
    xRgemma   =  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma   =   -0.5*carbondistance
 
    Expon_KRalpha      =exp(AI*(xk*xRalpha+yk*yRalpha ))
    Expon_KRbeta       =exp(AI*(xk*xRbeta+yk*yRbeta ))
    Expon_KRgemma      =exp(AI*(xk*xRgemma+yk*yRgemma))
    
    phase =exp(-AI*(xk*xRbeta+yk*yRbeta ))    !!!!   Bloch theory in one Unit Cell

    Hamiltonian1inter(9,12)=tpp(1,1,9)*Expon_KRalpha+tpp(1,1,8)*Expon_KRgemma+tpp(1,1,7)*Expon_KRbeta
    Hamiltonian1inter(9,13)=tpp(1,2,9)*Expon_KRalpha+tpp(1,2,8)*Expon_KRgemma+tpp(1,2,7)*Expon_KRbeta
    Hamiltonian1inter(9,14)=tpp(1,3,9)*Expon_KRalpha+tpp(1,3,8)*Expon_KRgemma+tpp(1,3,7)*Expon_KRbeta	

    Hamiltonian1inter(10,12)=tpp(2,1,9)*Expon_KRalpha+tpp(2,1,8)*Expon_KRgemma+tpp(2,1,7)*Expon_KRbeta
    Hamiltonian1inter(10,13)=tpp(2,2,9)*Expon_KRalpha+tpp(2,2,8)*Expon_KRgemma+tpp(2,2,7)*Expon_KRbeta
    Hamiltonian1inter(10,14)=tpp(2,3,9)*Expon_KRalpha+tpp(2,3,8)*Expon_KRgemma+tpp(2,3,7)*Expon_KRbeta	


    Hamiltonian1inter(11,12)=tpp(3,1,9)*Expon_KRalpha+tpp(3,1,8)*Expon_KRgemma+tpp(3,1,7)*Expon_KRbeta
    Hamiltonian1inter(11,13)=tpp(3,2,9)*Expon_KRalpha+tpp(3,2,8)*Expon_KRgemma+tpp(3,2,7)*Expon_KRbeta
    Hamiltonian1inter(11,14)=tpp(3,3,9)*Expon_KRalpha+tpp(3,3,8)*Expon_KRgemma+tpp(3,3,7)*Expon_KRbeta	

    Hamiltonian1inter( 9:11, 12:14)=phase*Hamiltonian1inter( 9:11, 12:14)

    Hamiltonian1inter(9,14)=-Hamiltonian1inter(9,14)
    Hamiltonian1inter(10,14)=-Hamiltonian1inter(10,14)
    Hamiltonian1inter(11,12)=-Hamiltonian1inter(11,12)
    Hamiltonian1inter(11,13)=-Hamiltonian1inter(11,13)


!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!           Fu ZHI
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Chu SHI Values
do i= 1, neigenvalues
        Hamiltonian(  i ,  1: neigenvalues  )= ( 0.d0, 0.d0 )
end do


If ( neigenvalues==11 ) then
        lneigenvalues= 11
        Hamiltonian(  lneigenvalues-10, lneigenvalues-10: lneigenvalues  )= Hamiltonian0( 1, 1:11  )
        Hamiltonian(  lneigenvalues-9, lneigenvalues-9: lneigenvalues  )= Hamiltonian0( 2, 2:11  )
        Hamiltonian(  lneigenvalues-8, lneigenvalues-8: lneigenvalues  )= Hamiltonian0( 3, 3:11  )
        Hamiltonian(  lneigenvalues-7, lneigenvalues-7: lneigenvalues  )= Hamiltonian0( 4, 4:11  )
        Hamiltonian(  lneigenvalues-6, lneigenvalues-6: lneigenvalues  )= Hamiltonian0( 5, 5:11  )
        Hamiltonian(  lneigenvalues-5, lneigenvalues-5: lneigenvalues  )= Hamiltonian0( 6, 6:11  )
        Hamiltonian(  lneigenvalues-4, lneigenvalues-4: lneigenvalues  )= Hamiltonian0( 7, 7:11  )
        Hamiltonian(  lneigenvalues-3, lneigenvalues-3: lneigenvalues  )= Hamiltonian0( 8, 8:11  )
        Hamiltonian(  lneigenvalues-2, lneigenvalues-2: lneigenvalues  )= Hamiltonian0( 9, 9:11  )
        Hamiltonian(  lneigenvalues-1, lneigenvalues-1: lneigenvalues  )= Hamiltonian0( 10, 10:11  )
        Hamiltonian(  lneigenvalues , lneigenvalues  )= Hamiltonian0( 11, 11  )

elseif (  neigenvalues==22      ) then

        lneigenvalues= 22 
        Hamiltonian(  lneigenvalues-21, lneigenvalues-21: lneigenvalues-11  )= Hamiltonian0( 1, 1:11  )
        Hamiltonian(  lneigenvalues-20, lneigenvalues-20: lneigenvalues-11  )= Hamiltonian0( 2, 2:11  )
        Hamiltonian(  lneigenvalues-19, lneigenvalues-19: lneigenvalues-11  )= Hamiltonian0( 3, 3:11  )
        Hamiltonian(  lneigenvalues-18, lneigenvalues-18: lneigenvalues-11  )= Hamiltonian0( 4, 4:11  )
        Hamiltonian(  lneigenvalues-17, lneigenvalues-17: lneigenvalues-11  )= Hamiltonian0( 5, 5:11  )
        Hamiltonian(  lneigenvalues-16, lneigenvalues-16: lneigenvalues-11  )= Hamiltonian0( 6, 6:11  )
        Hamiltonian(  lneigenvalues-15, lneigenvalues-15: lneigenvalues-11  )= Hamiltonian0( 7, 7:11  )
        Hamiltonian(  lneigenvalues-14, lneigenvalues-14: lneigenvalues-11  )= Hamiltonian0( 8, 8:11  )
        Hamiltonian(  lneigenvalues-13, lneigenvalues-13: lneigenvalues-11  )= Hamiltonian0( 9, 9:11  )
        Hamiltonian(  lneigenvalues-12, lneigenvalues-12: lneigenvalues-11  )= Hamiltonian0( 10, 10:11  )
        Hamiltonian(  lneigenvalues-11 , lneigenvalues-11  )= Hamiltonian0( 11, 11  )
    
        Hamiltonian(  lneigenvalues-10, lneigenvalues-10: lneigenvalues  )= Hamiltonian1( 1, 1:11  )
        Hamiltonian(  lneigenvalues-9, lneigenvalues-9: lneigenvalues  )= Hamiltonian1( 2, 2:11  )
        Hamiltonian(  lneigenvalues-8, lneigenvalues-8: lneigenvalues  )= Hamiltonian1( 3, 3:11  )
        Hamiltonian(  lneigenvalues-7, lneigenvalues-7: lneigenvalues  )= Hamiltonian1( 4, 4:11  )
        Hamiltonian(  lneigenvalues-6, lneigenvalues-6: lneigenvalues  )= Hamiltonian1( 5, 5:11  )
        Hamiltonian(  lneigenvalues-5, lneigenvalues-5: lneigenvalues  )= Hamiltonian1( 6, 6:11  )
        Hamiltonian(  lneigenvalues-4, lneigenvalues-4: lneigenvalues  )= Hamiltonian1( 7, 7:11  )
        Hamiltonian(  lneigenvalues-3, lneigenvalues-3: lneigenvalues  )= Hamiltonian1( 8, 8:11  )
        Hamiltonian(  lneigenvalues-2, lneigenvalues-2: lneigenvalues  )= Hamiltonian1( 9, 9:11  )
        Hamiltonian(  lneigenvalues-1, lneigenvalues-1: lneigenvalues  )= Hamiltonian1( 10, 10:11  )
        Hamiltonian(  lneigenvalues , lneigenvalues  )= Hamiltonian1( 11, 11  )

        Hamiltonian(  lneigenvalues-21 : lneigenvalues-19,    lneigenvalues-10: lneigenvalues-8  )=Hamiltonian0inter(9:11 , 12:14) 
        Hamiltonian(  lneigenvalues-13 : lneigenvalues-11,    lneigenvalues-2: lneigenvalues  )=Hamiltonian1inter(9:11 , 12:14) 

end if


!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!           Vonset
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    do i=1, neigenvalues
        m=mod( i, 11 )
        j = (i-1)/11+1
        if(  m>0 .and. m<4   ) then
            Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 3,  j )
            else if (  m>3 .and. m<9 ) then
            Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 2,  j )
            else if (  (m>8 .and. m<11) .or. (m==0)  ) then
            Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 1,  j )
        end if
    !~ print*
    !~ print*, Hamiltonian( i, i )
    end do

!!==================================================================
!                            End     Fu ZHI  
!!==================================================================
!write(*, '(2E16.8, 2E16.8, /, 2E16.8, 2E16.8)')   ((Hamiltonian(M,N), M=1,2),N=1,2)
!write(*, *)  Hamiltonian
!open (unit =UNIT4, file = FILENAME4,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)

K=0
do i=1,neigenvalues
   do j=1,i
   H(K)=Hamiltonian(j,i)
   K=K+1
   end do
end do   

        call ZHPEVX('V','A','U',neigenvalues,H,Vlower,Vupper,1,neigenvalues,abstol,Mfound,eig,vec,neigenvalues,cwork,aux,iwork,ifail,info)
        if(info.ne.0) then
        write(6,*) 'info: ',info
        write(6,*) 'Mfound ',Mfound
        do i=0,3
        write(6,*) 'i,ifail: ',i,ifail(i)
        enddo
        stop
        endif
!     write(6,*) eig
!	write(6,*)

end subroutine Hamiltonianofxs2VII

















!~ ==========================================================

    !~ Hamiltonian of XS2    2II version

!~ ==========================================================
subroutine Hamiltonianofxs2II( neigenvalues,  XK, Yk, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eig, vec) !  ,v1,v2,v3
use moduleone
Implicit real (kind=DP) (a-h,o-z) 
!----------------------------------------------------------------------------------------------------------------------
        integer(kind=4):: neigenvalues
        real(kind=DP),dimension(0:neigenvalues-1):: eig
        complex(kind=DP):: vec(0:neigenvalues-1,0:neigenvalues-1)
        complex(kind=DP):: H(0:(neigenvalues*(neigenvalues+1))/2-1)
        complex(kind=DP):: cwork(0:2*neigenvalues-1)
        real(kind=DP):: aux(0:7*neigenvalues-1)
        integer:: iwork(0:5*neigenvalues-1),ifail(0:neigenvalues-1)
        real(kind=DP),parameter:: abstol=2*tiny(abstol)  ! see ZHPEVX

        !!! =========  Define more parameters     =====
		integer(kind=4):: NdimensionVonset,  lneigenvalues
		real(kind=DP):: XK, Yk
		complex(kind=DP), dimension(neigenvalues,neigenvalues)::Hamiltonian, Hamiltonian0   
		complex(kind=DP)::Expon_KRalpha, Expon_KRbeta,Expon_KRgemma,  phase 
		complex(kind=DP):: tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)
		real(kind=DP)::Delta0,Delta1,Delta2,DeltaP,DeltaZ
		real(kind=DP):: v12, v23
		real(kind=DP), dimension ( 1:3, 1:(neigenvalues/11) ):: Vonset    		!~! Allocate(  Vonset  ( 1:3, Nlayer  )     )
		
!----------------------------------------------------------------------------------------------------------------------

!======================================================================================
!~ CALL MKL_SET_NUM_THREADS(1)
!~ !!! 

xRalpha   = -0.5*avalue          ; yRalpha   =    0.5*sqrt(3.0_DP)*avalue
xRbeta    =   avalue                 ; yRbeta    =    0.0
xRgemma= -0.5*avalue          ; yRgemma=   -0.5*sqrt(3.0_DP)*avalue

CS_KRalpha      =   COS(xk*xRalpha+yk*yRalpha)
CS_KRbeta       =   COS(xk*xRbeta+yk*yRbeta)
CS_KRgemma  =   COS(xk*xRgemma+yk*yRgemma)

!!%%%%%%%%%%%%%%        ST-----ST            %%%%%%%%%%%%%%%%

Hamiltonian0(1,1)=2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP !+v23
Hamiltonian0(1,2)=2.0_DP*( tpp(1,2,6)*CS_KRalpha+tpp(1,2,4)*CS_KRbeta+tpp(1,2,2)*CS_KRgemma )
Hamiltonian0(1,3)=2.0_DP*( tpp(1,3,6)*CS_KRalpha+tpp(1,3,4)*CS_KRbeta+tpp(1,3,2)*CS_KRgemma )

Hamiltonian0(2,2)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP !+v23
Hamiltonian0(2,3)=2.0_DP*( tpp(2,3,6)*CS_KRalpha+tpp(2,3,4)*CS_KRbeta+tpp(2,3,2)*CS_KRgemma )
Hamiltonian0(3,3)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ !+v23

!!%%%%%%%%%%%%%%        SB-----SB            %%%%%%%%%%%%%%%%

Hamiltonian0(9,9)  =2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP !-v12
Hamiltonian0(9,10)=Hamiltonian0(1,2)
Hamiltonian0(9,11)=Hamiltonian0(1,3)

Hamiltonian0(10,10)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP !-v12
Hamiltonian0(10,11)=Hamiltonian0(2,3)
Hamiltonian0(11,11)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ !-v12

!!%%%%%%%%%%%%%%        Mo-----Mo            %%%%%%%%%%%%%%%%

Hamiltonian0(4,4)=2.0_DP*( tdd(1,1,6)*CS_KRalpha+tdd(1,1,4)*CS_KRbeta+tdd(1,1,2)*CS_KRgemma )+Delta2  !  v2=0 eV
Hamiltonian0(4,5)=2.0_DP*( tdd(1,2,6)*CS_KRalpha+tdd(1,2,4)*CS_KRbeta+tdd(1,2,2)*CS_KRgemma )
Hamiltonian0(4,6)=2.0_DP*( tdd(1,3,6)*CS_KRalpha+tdd(1,3,4)*CS_KRbeta+tdd(1,3,2)*CS_KRgemma )
Hamiltonian0(4,7)=2.0_DP*( tdd(1,4,6)*CS_KRalpha+tdd(1,4,4)*CS_KRbeta+tdd(1,4,2)*CS_KRgemma )
Hamiltonian0(4,8)=2.0_DP*( tdd(1,5,6)*CS_KRalpha+tdd(1,5,4)*CS_KRbeta+tdd(1,5,2)*CS_KRgemma )

Hamiltonian0(5,5)=2.0_DP*( tdd(2,2,6)*CS_KRalpha+tdd(2,2,4)*CS_KRbeta+tdd(2,2,2)*CS_KRgemma )+Delta1  !  v2=0 eV
Hamiltonian0(5,6)=2.0_DP*( tdd(2,3,6)*CS_KRalpha+tdd(2,3,4)*CS_KRbeta+tdd(2,3,2)*CS_KRgemma )
Hamiltonian0(5,7)=2.0_DP*( tdd(2,4,6)*CS_KRalpha+tdd(2,4,4)*CS_KRbeta+tdd(2,4,2)*CS_KRgemma )
Hamiltonian0(5,8)=2.0_DP*( tdd(2,5,6)*CS_KRalpha+tdd(2,5,4)*CS_KRbeta+tdd(2,5,2)*CS_KRgemma )

Hamiltonian0(6,6)=2.0_DP*( tdd(3,3,6)*CS_KRalpha+tdd(3,3,4)*CS_KRbeta+tdd(3,3,2)*CS_KRgemma )+Delta1    !  v2=0 eV
Hamiltonian0(6,7)=2.0_DP*( tdd(3,4,6)*CS_KRalpha+tdd(3,4,4)*CS_KRbeta+tdd(3,4,2)*CS_KRgemma )
Hamiltonian0(6,8)=2.0_DP*( tdd(3,5,6)*CS_KRalpha+tdd(3,5,4)*CS_KRbeta+tdd(3,5,2)*CS_KRgemma )

Hamiltonian0(7,7)=2.0_DP*( tdd(4,4,6)*CS_KRalpha+tdd(4,4,4)*CS_KRbeta+tdd(4,4,2)*CS_KRgemma )+Delta2    !  v2=0 eV
Hamiltonian0(7,8)=2.0_DP*( tdd(4,5,6)*CS_KRalpha+tdd(4,5,4)*CS_KRbeta+tdd(4,5,2)*CS_KRgemma )
Hamiltonian0(8,8)=2.0_DP*( tdd(5,5,6)*CS_KRalpha+tdd(5,5,4)*CS_KRbeta+tdd(5,5,2)*CS_KRgemma )+Delta0    !  v2=0 eV

!!%%%%%%%%%%%%%%        Mo-----> ST    dega+        %%%%%%%%%%%%%%%%

xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
xRbeta    =   0                                                         ; yRbeta    =    carbondistance          
xRgemma=  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma=    -0.5*carbondistance 

Expon_KRalpha      =exp(-AI*(xk*xRalpha+yk*yRalpha ))
Expon_KRbeta       =exp(-AI*(xk*xRbeta+yk*yRbeta ))
Expon_KRgemma   =exp(-AI*(xk*xRgemma+yk*yRgemma))

phase =exp(AI*(xk*xRbeta+yk*yRbeta ))     !!!!   Bloch theory in one Unit Cell

Hamiltonian0(1,4)=tdp(1,1,3)*Expon_KRalpha+tdp(1,1,2)*Expon_KRgemma+tdp(1,1,1)*Expon_KRbeta
Hamiltonian0(1,5)=tdp(1,2,3)*Expon_KRalpha+tdp(1,2,2)*Expon_KRgemma+tdp(1,2,1)*Expon_KRbeta
Hamiltonian0(1,6)=tdp(1,3,3)*Expon_KRalpha+tdp(1,3,2)*Expon_KRgemma+tdp(1,3,1)*Expon_KRbeta 
Hamiltonian0(1,7)=tdp(1,4,3)*Expon_KRalpha+tdp(1,4,2)*Expon_KRgemma+tdp(1,4,1)*Expon_KRbeta  
Hamiltonian0(1,8)=tdp(1,5,3)*Expon_KRalpha+tdp(1,5,2)*Expon_KRgemma+tdp(1,5,1)*Expon_KRbeta 

Hamiltonian0(2,4)=tdp(2,1,3)*Expon_KRalpha+tdp(2,1,2)*Expon_KRgemma+tdp(2,1,1)*Expon_KRbeta
Hamiltonian0(2,5)=tdp(2,2,3)*Expon_KRalpha+tdp(2,2,2)*Expon_KRgemma+tdp(2,2,1)*Expon_KRbeta  
Hamiltonian0(2,6)=tdp(2,3,3)*Expon_KRalpha+tdp(2,3,2)*Expon_KRgemma+tdp(2,3,1)*Expon_KRbeta
Hamiltonian0(2,7)=tdp(2,4,3)*Expon_KRalpha+tdp(2,4,2)*Expon_KRgemma+tdp(2,4,1)*Expon_KRbeta 
Hamiltonian0(2,8)=tdp(2,5,3)*Expon_KRalpha+tdp(2,5,2)*Expon_KRgemma+tdp(2,5,1)*Expon_KRbeta 

Hamiltonian0(3,4)=tdp(3,1,3)*Expon_KRalpha+tdp(3,1,2)*Expon_KRgemma+tdp(3,1,1)*Expon_KRbeta  
Hamiltonian0(3,5)=tdp(3,2,3)*Expon_KRalpha+tdp(3,2,2)*Expon_KRgemma+tdp(3,2,1)*Expon_KRbeta 
Hamiltonian0(3,6)=tdp(3,3,3)*Expon_KRalpha+tdp(3,3,2)*Expon_KRgemma+tdp(3,3,1)*Expon_KRbeta 
Hamiltonian0(3,7)=tdp(3,4,3)*Expon_KRalpha+tdp(3,4,2)*Expon_KRgemma+tdp(3,4,1)*Expon_KRbeta  
Hamiltonian0(3,8)=tdp(3,5,3)*Expon_KRalpha+tdp(3,5,2)*Expon_KRgemma+tdp(3,5,1)*Expon_KRbeta
Hamiltonian0(1:3 , 4:8 )=Hamiltonian0(1:3 , 4:8 )*phase

!!%%%%%%%%%%%%%%        SB -----> Mo            %%%%%%%%%%%%%%%%

xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
xRbeta    =   0                                                         ; yRbeta    =    carbondistance          
xRgemma=  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma=    -0.5*carbondistance 

Expon_KRalpha      =exp(AI*(xk*xRalpha+yk*yRalpha ))
Expon_KRbeta       =exp(AI*(xk*xRbeta+yk*yRbeta))
Expon_KRgemma   =exp(AI*(xk*xRgemma+yk*yRgemma))

phase =exp(-AI*(xk*xRbeta+yk*yRbeta ))    !!!!   Bloch theory in one Unit Cell

Hamiltonian0(4,9)  =tdp(1,1,6)*Expon_KRalpha+tdp(1,1,5)*Expon_KRgemma+tdp(1,1,4)*Expon_KRbeta
Hamiltonian0(5,9)  =tdp(1,2,6)*Expon_KRalpha+tdp(1,2,5)*Expon_KRgemma+tdp(1,2,4)*Expon_KRbeta 
Hamiltonian0(6,9)  =tdp(1,3,6)*Expon_KRalpha+tdp(1,3,5)*Expon_KRgemma+tdp(1,3,4)*Expon_KRbeta 
Hamiltonian0(7,9)  =tdp(1,4,6)*Expon_KRalpha+tdp(1,4,5)*Expon_KRgemma+tdp(1,4,4)*Expon_KRbeta 
Hamiltonian0(8,9)  =tdp(1,5,6)*Expon_KRalpha+tdp(1,5,5)*Expon_KRgemma+tdp(1,5,4)*Expon_KRbeta 

Hamiltonian0(4,10)  =tdp(2,1,6)*Expon_KRalpha+tdp(2,1,5)*Expon_KRgemma+tdp(2,1,4)*Expon_KRbeta 
Hamiltonian0(5,10)  =tdp(2,2,6)*Expon_KRalpha+tdp(2,2,5)*Expon_KRgemma+tdp(2,2,4)*Expon_KRbeta 
Hamiltonian0(6,10)  =tdp(2,3,6)*Expon_KRalpha+tdp(2,3,5)*Expon_KRgemma+tdp(2,3,4)*Expon_KRbeta 
Hamiltonian0(7,10)  =tdp(2,4,6)*Expon_KRalpha+tdp(2,4,5)*Expon_KRgemma+tdp(2,4,4)*Expon_KRbeta 
Hamiltonian0(8,10)  =tdp(2,5,6)*Expon_KRalpha+tdp(2,5,5)*Expon_KRgemma+tdp(2,5,4)*Expon_KRbeta 

Hamiltonian0(4,11)  =tdp(3,1,6)*Expon_KRalpha+tdp(3,1,5)*Expon_KRgemma+tdp(3,1,4)*Expon_KRbeta 
Hamiltonian0(5,11)  =tdp(3,2,6)*Expon_KRalpha+tdp(3,2,5)*Expon_KRgemma+tdp(3,2,4)*Expon_KRbeta 
Hamiltonian0(6,11)  =tdp(3,3,6)*Expon_KRalpha+tdp(3,3,5)*Expon_KRgemma+tdp(3,3,4)*Expon_KRbeta 
Hamiltonian0(7,11)  =tdp(3,4,6)*Expon_KRalpha+tdp(3,4,5)*Expon_KRgemma+tdp(3,4,4)*Expon_KRbeta 
Hamiltonian0(8,11)  =tdp(3,5,6)*Expon_KRalpha+tdp(3,5,5)*Expon_KRgemma+tdp(3,5,4)*Expon_KRbeta 
Hamiltonian0(4:8 , 9:11)  =Hamiltonian0(4:8 , 9:11)*phase 

!!%%%%%%%%%%%%%%        SB -----> ST           %%%%%%%%%%%%%%%%
!!    SB----> ST

    Hamiltonian0(1,9)=tpp(1,1,10)
    Hamiltonian0(1,10)=tpp(1,2,10)
    Hamiltonian0(1,11)=tpp(1,3,10)

    Hamiltonian0(2,9)=tpp(2,1,10)
    Hamiltonian0(2,10)=tpp(2,2,10)
    Hamiltonian0(2,11)=tpp(2,3,10)

    Hamiltonian0(3,9)=tpp(3,1,10)
    Hamiltonian0(3,10)=tpp(3,2,10)
    Hamiltonian0(3,11)=tpp(3,3,10)

!!%%%%%%%%%%%%%%        ST one layer -----> SB up layer    daga+       %%%%%%%%%%%%%%%%

 
    xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
    xRbeta    =   0                                                         ; yRbeta    =    carbondistance          
    xRgemma=  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma=    -0.5*carbondistance
 
	Expon_KRalpha      =exp(-AI*(xk*xRalpha+yk*yRalpha ))
    Expon_KRbeta       =exp(-AI*(xk*xRbeta+yk*yRbeta ))
    Expon_KRgemma   =exp(-AI*(xk*xRgemma+yk*yRgemma))
    
	phase =exp(AI*(xk*xRalpha+yk*yRalpha ))     !!!!   Bloch theory in one Unit Cell
	
	Hamiltonian0(9,12)=tpp(1,1,9)*Expon_KRalpha+tpp(1,1,8)*Expon_KRgemma+tpp(1,1,7)*Expon_KRbeta
	Hamiltonian0(9,13)=tpp(1,2,9)*Expon_KRalpha+tpp(1,2,8)*Expon_KRgemma+tpp(1,2,7)*Expon_KRbeta
	Hamiltonian0(9,14)=tpp(1,3,9)*Expon_KRalpha+tpp(1,3,8)*Expon_KRgemma+tpp(1,3,7)*Expon_KRbeta	

	Hamiltonian0(10,12)=tpp(2,1,9)*Expon_KRalpha+tpp(2,1,8)*Expon_KRgemma+tpp(2,1,7)*Expon_KRbeta
	Hamiltonian0(10,13)=tpp(2,2,9)*Expon_KRalpha+tpp(2,2,8)*Expon_KRgemma+tpp(2,2,7)*Expon_KRbeta
	Hamiltonian0(10,14)=tpp(2,3,9)*Expon_KRalpha+tpp(2,3,8)*Expon_KRgemma+tpp(2,3,7)*Expon_KRbeta	


	Hamiltonian0(11,12)=tpp(3,1,9)*Expon_KRalpha+tpp(3,1,8)*Expon_KRgemma+tpp(3,1,7)*Expon_KRbeta
	Hamiltonian0(11,13)=tpp(3,2,9)*Expon_KRalpha+tpp(3,2,8)*Expon_KRgemma+tpp(3,2,7)*Expon_KRbeta
	Hamiltonian0(11,14)=tpp(3,3,9)*Expon_KRalpha+tpp(3,3,8)*Expon_KRgemma+tpp(3,3,7)*Expon_KRbeta	 
    Hamiltonian0(9:11 , 12:14)  =Hamiltonian0(9:11 , 12:14)*phase 
 
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                 Fu ZHI
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Chu SHI Values
do i= 1, neigenvalues-1
        Hamiltonian(  i ,  i: neigenvalues  )= ( 0.d0, 0.d0 )
end do

do lneigenvalues= neigenvalues,  11, -11
    Hamiltonian(  lneigenvalues-10, lneigenvalues-10: lneigenvalues  )= Hamiltonian0( 1, 1:11  )
    Hamiltonian(  lneigenvalues-9, lneigenvalues-9: lneigenvalues  )= Hamiltonian0( 2, 2:11  )
    Hamiltonian(  lneigenvalues-8, lneigenvalues-8: lneigenvalues  )= Hamiltonian0( 3, 3:11  )
    Hamiltonian(  lneigenvalues-7, lneigenvalues-7: lneigenvalues  )= Hamiltonian0( 4, 4:11  )
    Hamiltonian(  lneigenvalues-6, lneigenvalues-6: lneigenvalues  )= Hamiltonian0( 5, 5:11  )
    Hamiltonian(  lneigenvalues-5, lneigenvalues-5: lneigenvalues  )= Hamiltonian0( 6, 6:11  )
    Hamiltonian(  lneigenvalues-4, lneigenvalues-4: lneigenvalues  )= Hamiltonian0( 7, 7:11  )	
    Hamiltonian(  lneigenvalues-3, lneigenvalues-3: lneigenvalues  )= Hamiltonian0( 8, 8:11  )	
    Hamiltonian(  lneigenvalues-2, lneigenvalues-2: lneigenvalues  )= Hamiltonian0( 9, 9:11  )
    Hamiltonian(  lneigenvalues-1, lneigenvalues-1: lneigenvalues  )= Hamiltonian0( 10, 10:11  )	
    Hamiltonian(  lneigenvalues , lneigenvalues  )= Hamiltonian0( 11, 11  )	
	 if (   lneigenvalues>11 ) then
    Hamiltonian(  lneigenvalues-13 : lneigenvalues-11,    lneigenvalues-10: lneigenvalues-8  )=Hamiltonian0(9:11 , 12:14) 
	 end if
end do	

		!~! Allocate(  Vonset  ( 1:3, Nlayer  )     )
		
	do i=1, neigenvalues
		!~ !print*, (i-1)/11+1
        k=mod( i, 11)
		if(  0<k .and. k<=3   ) then
			Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 1,  (i-1)/11+1 )
			elseif (  3<k .and. k<=8 ) then
			Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 2,  (i-1)/11+1 )
			elseif (  (8<k .and. k<=10) .or. k==0  ) then
			Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 3,  (i-1)/11+1 )			
		end if
	end do

	!~ ==============     check            ====================
	!~ do i=1, Nlayer 
		!~ Vonset( 1, i  )=real(i)*10.0+1.0	
		!~ Vonset( 2, i  )=real(i)*10.0+2.0	
		!~ Vonset( 3, i  )=real(i)*10.0+3.0		
		!~ print*, Vonset( :, i  )
    !~ end do
	!~ Hamiltonian( 1:neigenvalues, 1:neigenvalues)=0.d0

	!~ do i=1, neigenvalues
		!print*, (i-1)/11+1 
		!~ if(  0< mod( i, 11 )<=3   ) then
			!~ Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 1,  (i-1)/11+1 )
			!~ elseif (  3< mod( i, 11 )<=8 ) then
			!~ Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 2,  (i-1)/11+1 )
			!~ elseif (  8< mod( i, 11 )<=10 .or. mod( i, 11 )==0  ) then
			!~ Hamiltonian( i, i )=Hamiltonian( i, i )+Vonset( 3,  (i-1)/11+1 )			
		!~ end if
	!~ end do
	
	!~ do  j=1, neigenvalues
		!~ do i=1, neigenvalues
		!~ write(8,  '( a,I2, I2,a, f20.8 )'  ) "H_( ", i, j, " ) = ",  Hamiltonian( i, j )
		!!write(8, '( a )'    ) 
		!~ end do
	!~ end do
	!~ ==============     End  check            ================

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                            End     Fu ZHI  
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!~ =============================================================
!write(*, '(2E16.8, 2E16.8, /, 2E16.8, 2E16.8)')   ((Hamiltonian(M,N), M=1,2),N=1,2)
!write(*, *)  Hamiltonian
!open (unit =UNIT4, file = FILENAME4,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)

K=0
do i=1,neigenvalues
   do j=1,i
   H(K)=Hamiltonian(j,i)
   K=K+1
   end do
end do   

!~ print*, K
!~ write(*,  '( 2(2xe20.10) )'  )  H(1)
!~ print*
!~ write(*,  '( 2(2xe20.10) )'  )  H(65)
!~ print*
!!!~ Should be same
!~ print*
!~ write(*,  '( 2(2xe20.10) )'  )  H(89)
!~ print*
!~ write(*,  '( 2(2xe20.10) )'  )  H(K-1)
! write(UNIT4, '(2E16.8)')  H(K)

        call ZHPEVX('V','A','U',neigenvalues,H,Vlower,Vupper,1,neigenvalues,abstol,Mfound,eig,vec,neigenvalues,cwork,aux,iwork,ifail,info)
        if(info.ne.0) then
        write(6,*) 'info: ',info
        write(6,*) 'Mfound ',Mfound
        do i=0,3
        write(6,*) 'i,ifail: ',i,ifail(i)
        enddo
        stop
        endif
!     write(6,*) eig
!	write(6,*)

end subroutine Hamiltonianofxs2II







!~ ==========================================================
                !~ Hamiltonian of Graphene
!~ ==========================================================
subroutine Hamiltonianof1layergraphene(tp, xk, yk, eig, vec) 
use moduleone
Implicit real (kind=DP) (a-h,o-z) 
!----------------------------------------------------------------------------------------------------------------------
        integer,parameter:: neigenvalues=2
        real(kind=DP),dimension(0:neigenvalues-1):: eig
        complex(kind=DP):: vec(0:neigenvalues-1,0:neigenvalues-1)
        complex(kind=DP):: H(0:(neigenvalues*(neigenvalues+1))/2-1)
        complex(kind=DP):: cwork(0:2*neigenvalues-1)
        real(kind=DP):: aux(0:7*neigenvalues-1)
        integer:: iwork(0:5*neigenvalues-1),ifail(0:neigenvalues-1)
        real(kind=DP),parameter:: abstol=2*tiny(abstol)  ! see ZHPEVX
		
        real(kind=DP),intent(IN):: xk, yk, tp
        complex(kind=DP), dimension(neigenvalues,neigenvalues)::Hamiltonian
        complex(kind=DP)::Expon_KRalpha, Expon_KRbeta,Expon_KRgemma 
        complex(kind=DP)::Expon_A1, 	Expon_A2	
	
!~ CALL MKL_SET_NUM_THREADS(1)
!!!!  avalue=sqrt(3)*carbondistance
!======================================================================================
xA1=0.5d0*sqrt(3.0_DP)*carbondistance ;    yA1=1.5d0*carbondistance 
xA2=-0.5d0*sqrt(3.0_DP)*carbondistance ;   yA2=1.5d0*carbondistance 

Expon_A1=exp(-AI*(xk*xA1+yk*yA1 ))
Expon_A2=exp(-AI*(xk*xA2+yk*yA2 ))

Hamiltonian(1,1)  =cmplx(0.0d0, 0.0d0)
Hamiltonian(1,2)  =-tp*Expon_A1-tp*Expon_A2-tp*1.0d0 
Hamiltonian(2,2)  =cmplx(0.0d0, 0.0d0)
!====================================================================================================
!write(*, '(2E16.8, 2E16.8, /, 2E16.8, 2E16.8)')   ((Hamiltonian(M,N), M=1,2),N=1,2)
!write(*, *)  Hamiltonian
!open (unit =UNIT4, file = FILENAME4,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)

K=0
do i=1,neigenvalues
   do j=1,i
   H(K)=Hamiltonian(j,i)
   K=K+1
   end do
end do   
  ! write(UNIT4, '(2E16.8)')  H(K)
        call ZHPEVX('V','A','U',neigenvalues,H,Vlower,Vupper,1,neigenvalues,abstol,Mfound,eig,vec,neigenvalues,cwork,aux,iwork,ifail,info)
        if(info.ne.0) then
        write(6,*) 'info: ',info
        write(6,*) 'Mfound ',Mfound
        do i=0,3
        write(6,*) 'i,ifail: ',i,ifail(i)
        enddo
        stop
        endif
end
!end subroutine 

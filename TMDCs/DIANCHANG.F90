    program DIAN_CHANG_QIANG_DU
    implicit real(kind=8) (a-h, o-z)
    
    character (len=50) :: material_name
    REAL(KIND=8) KAPPA, nprime
    real(kind=8), parameter :: pi=3.141592653589793238460d0,  sq3=dsqrt(3.0d0)
    !~ real(kind=8), parameter :: AMos2=3.160*10.0**(-10) 
    !~ real(kind=8), parameter :: AWs2=3.153*10.0**(-10) 
    !~ real(kind=8), parameter :: AMoSe2=3.288*10.0**(-10) 
    !~ real(kind=8), parameter :: AWSe2=3.260*10.0**(-10) 
    
    real(kind=8), parameter :: AMos2=3.157*10.0**(-10) 
    real(kind=8), parameter :: AWs2=3.160*10.0**(-10) 
    real(kind=8), parameter :: AMoSe2=3.290*10.0**(-10) 
    real(kind=8), parameter :: AWSe2=3.290*10.0**(-10)     

    real(kind=8), parameter ::  Sunitcellmos2=0.50d0*sq3*AMos2**2
    real(kind=8), parameter ::  SunitcellWs2=0.50d0*sq3*AWs2**2      
    real(kind=8), parameter ::  Sunitcellgraphene=2.619*2.0*10.0**(-20)
    real(kind=8), parameter ::  SunitcellmoSe2=0.50d0*sq3*AMoSe2**2 
    real(kind=8), parameter ::  SunitcellWSe2=0.50d0*sq3*AWSe2**2 

    character (len=*), parameter:: material_nameMoS2="MoS2", material_nameWS2="WS2"
    character (len=*), parameter:: material_nameMoSe2="MoSe2", material_nameWSe2="WSe2"
    character (len=*), parameter:: material_nameGraphene="Graphene"
    character (len=*), parameter:: material_nameTriGraphene="TriGraphene"

    PRINT*
    write(*, '( A )' ), "Please enter Material name ( MoS2 =1; WS2 =2; MoSe2 =3; WSe2 =4 )::"
    READ*, SetTransitionMetalDichalcogenides
    
    IF (SetTransitionMetalDichalcogenides==0 ) then
        material_name=material_nameGraphene
        L=8    ! lens of the material name
        !~ write(*, '(a)')
        write(*, '(a)')"The object of material is Graphene. "
    else if (SetTransitionMetalDichalcogenides==1 ) then
        
        Sunitcell=SunitcellmoS2               !! from module1.f90
        material_name=material_nameMoS2
        L=4    ! lens of the material name
        !~ write(*, '(a)')
        write(*, '(a)')"The object of material is MoS2. "
        print*
        write(*, '( a, 1xES23.15E3 )' ) "Lattice Constant = ", AMos2
    else if (SetTransitionMetalDichalcogenides==2 ) then

        Sunitcell=SunitcellWS2
        material_name=material_nameWS2
        L=3    ! lens of the material name
        !~ write(*, '(a)')
        write(*, '(a)')"The object of material is WS2. "
        print*
        write(*, '( a, 1xES23.15E3 )' ) "Lattice Constant = ",  AWs2
    else if (SetTransitionMetalDichalcogenides==3 ) then
        Sunitcell=SunitcellMoSe2
        material_name=material_nameMoSe2
        L=5    ! lens of the material name
        !~ write(*, '(a)')
        write(*, '(a)')"The object of material is MoSe2. "
        print*
        write(*, '( a, 1xES23.15E3 )' ) "Lattice Constant = ", AMoSe2        
    else if (SetTransitionMetalDichalcogenides==4 ) then
        Sunitcell=SunitcellWSe2
        material_name=material_nameWSe2
        L=4    ! lens of the material name
        !~ write(*, '(a)')
        write(*, '(a)')"The object of material is WSe2. "
        print*
        write(*, '( a, 1xES23.15E3 )' ) "Lattice Constant = ", AWSe2        
    else
        !~ write(*, '(a)')
        write(*, '(a)') "Error!  What the name of Material? Please have a check!!!"	
        stop
    end if

    ELECTRON = 1.60217653*10.0**(-19) !!   unit ( C )
    kappa = (3.90d0+1.0d0)/2.0d0 ;  
    epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V

   !! nprime  UNIT   (     10^{13} cm^{-2}   )
    !~ PRINT*
    WRITE(*, '( A )' ) "PLEASE ENTER THE NPRIME VALUE ( ~ 10^{13} cm^{-2} )::"
    READ*,RNprimereal
    write(*, '( a,  f14.8,  a )') "The given carrier density=",  RNprimereal , "*E13/cm^2"
    
    RNprime=RNprimereal*Sunitcell*10.0d0**(17) 
    write(*, '( a)') ""
    write(*, '( a, 1xES23.15E3, a )'  ) "The given carrier density=",  RNprime, ", which is £¨ n in electrons per unit £©."

    
    Efiled=electron*RNprimereal*10.d0**(8)/(2.0d0*epsilong*kappa)    !! UNIT  (    V/nm   )
    PRINT*
    WRITE(*, '(  A, F9.5  )')"THE EFILED IN UNIT ( V/NM ) IS ", Efiled



    end program DIAN_CHANG_QIANG_DU
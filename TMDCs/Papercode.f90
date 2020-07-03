program mainTMD
use moduleone
use moduletwo
use qsort_c_module
Implicit real (kind=DP) (a-h,o-z) 
integer(KINd=4)::  SetTransitionMetalDichalcogenides, layers
integer(kind=4)::  numberofcarrier, judge, KNumberatfermi,  Icalculatepolarizationfunction
integer(kind=4)::  neNumberofWs2,  neNumberofMos2, neNumberofGraphene, NE0
integer(kind=4)::  neNumberofWSe2,  neNumberofMoSe2, Doping_way
integer(kind=4), allocatable, dimension(  :, :, : ):: MBianhua

real(kind=DP)::  Omega, Deltasmall, kappa
real(kind=DP)::  Delta0,Delta1,Delta2,DeltaP,DeltaZ, Atmdc
real(kind=DP)::  probability, RNprime, distanceoflayer,  distanceoflayerW
real(kind=DP)::  EMiu, EMiu0, Eerror
real(kind=DP)::  RdistanceX,  RdistanceY 
real(kind=DP)::  corelationeff
real(kind=DP)::  deltax, prob0
real(kind=DP)::  tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)

real(kind=DP), dimension( 0 : Blengths-1 ) :: Bkx, Bky
real(KIND=DP), ALLOCATABLE, dimension( : ) :: Densityofstatecheck,  DensityofstateTMDC, DosTMDC
real(kind=DP), ALLOCATABLE, dimension( : ) :: RNprimelayer, Potentialonlayer, checkprimelayer
real(kind=DP), ALLOCATABLE, dimension( : ) :: carrier_density 
real(kind=DP), ALLOCATABLE, dimension( : ) :: rankx, ranky, rankxQ, rankyQ, abs_q
real(kind=DP), ALLOCATABLE, dimension( : ) :: rankFSx, rankFSy, rankFSxnew, rankFSynew
real(kind=DP), Allocatable, dimension( : ) :: eigX
real(kind=DP), Allocatable, dimension( :, : ) :: Vonset
real(kind=DP), ALLOCATABLE, dimension( :, : ) :: Dos,  Densityofstate
real(kind=DP), ALLOCATABLE, dimension( :, : ) :: DoslayerU, DoslayerD
real(kind=DP), ALLOCATABLE, dimension( :, : ) :: DensityofstatelayerU, DensityofstatelayerD, Densityofstatelayer
real(kind=DP), ALLOCATABLE, dimension( :, : ) :: prob, ckeig, ckeigQ  
real(kind=DP), allocatable, dimension( :, : ) :: FermiKS, QKS, oldQKS
real(kind=DP), allocatable, dimension( :, :, : ) ::  Fenzi
complex(kind=DP):: DynamicalpolarizationPI, complexcoeff,  phase,  Cmiu
complex(kind=DP), allocatable, dimension( : ) :: CdielectricFun   !! epsilong(q, w)
complex(kind=DP), allocatable, dimension( : ) :: CWcoulomb_interaction
complex(kind=DP), allocatable, dimension( :, : ) :: vecX, ckvecM
complex(kind=DP), Allocatable, dimension( :, :, : ) :: ckvec, ckvecQ

character (len=50) :: test_name, material_name

real(kind=DP), external ::   FermiDistribution
real(kind=DP), external ::   FouriercomponentCoulombinteractionFME,  FouriercomponentCoulombinteractionFZ

!~ CALL MKL_SET_NUM_THREADS(8)
!~ ===============================================================================
    open (unit =UNIT9, file = FILENAME9,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)         !!===   open informatjion.out
    open (unit =UNIT7, file = FILENAME7,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)         !!===   Bandstructure0.dat 
    open (unit =UNIT200, file = FILENAME200,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)     !!===   dos0.dat
    print*, "=================================================="
    call timestamp (6)
    call timestamp (UNIT9)
    write(UNIT9, '( 79("=")  )') 
    write(6,  '(1x,a)' ) "The program is running"
    write(6,  '(1x,a)' ) "Waiting ......."
    write(6,  '(1x,a)' ) ' '
!===============================================================================
!========================     Read input.txt  and  FuZHI      ==================
!===============================================================================
    !!~ Intent(OUT) ::  all arguments
    call read_input_txt_file ( SetTransitionMetalDichalcogenides,  NKnumberBZ, Ndeltax,  Nlayer,  Bandstructure0_Windows_Y0, &
            & Bandstructure0_Windows_Y1, Gap1_Windows_Y0, Gap1_Windows_Y1, RNprimereal, Eerror, kappa, Icalculatepolarizationfunction, &
            & deltaX, neNumberofGraphene, neNumberofMoS2,  neNumberofWS2,  neNumberofMoSe2, neNumberofWSe2, Doping_way, NE0,  EMiu  )
    !~ SetTransitionMetalDichalcogenides
    !~ NKnumberBZ
    !~ Ndeltax
    !~ Nlayer
    !~ Bandstructure0_Windows_Y0 
    !~ Bandstructure0_Windows_Y1
    !~ RNprimereal
    !~ Eerror  (not used so far)
    !~ kappa
    !~ Icalculatepolarizationfunction
    !~ deltaX
    !~ neNumberofGraphene
    !~ neNumberofMoS2
    !~ neNumberofWS2
    !~ neNumberofMoSe2
    !~ neNumberofWSe2
    !~ NE0
    !~ EMiu
    print*
    print*, 'Doping_way =', Doping_way
    !~ stop
    allocate(    carrier_density                ( 1:3*Nlayer )    )
    allocate(    RNprimelayer                   ( 1:3*Nlayer )    )
    allocate(    Potentialonlayer               ( 1:3*Nlayer )    )
    allocate(    checkprimelayer                ( 1:3*Nlayer )    )
    allocate(    Densityofstatecheck        ( 0:Ndeltax )    )
    allocate(    DensityofstateTMDC         ( 0:Ndeltax )    )
    allocate(    DosTMDC                    ( 0:Ndeltax )    )
    Allocate(    Vonset             ( 1:3,         1:Nlayer  )    )
    allocate(    Dos                ( 0:Ndeltax,  1:3*Nlayer )    )
    allocate(    Densityofstate     ( 0:Ndeltax,  1:3*Nlayer )    )    
!===============================================================================
!~ ==============          Make_data_folder   name      ========================
!===============================================================================
    call Making_data_folder ( SetTransitionMetalDichalcogenides, neNumberofWSe2,  neNumberofMoSe2, &
            & neNumberofWS2,  neNumberofMoS2, neNumberofGraphene,   Nlayer,  RNprimereal, material_name, &
            & numberofcarrier,  Atmdc, Sunitcell,  distanceoflayer, test_name   )
    !!~ Intent(OUT) ::  numberofcarrier, Atmdc, Sunitcell,  distanceoflayer, test_name
    !~ print*
    !~ print*,"Atmdc=", Atmdc
!===============================================================================
!~ ======================           Getting Ks         =========================
!===============================================================================
    ConstantMoverN=1.2d0              !!!  (2.0d0/sq3=1.1547)   Width/Height
    N = NKnumberBZ                         !   Y axis    N=500*4
    M = NINT(N*ConstantMoverN)
    write(UNIT9,'(a)')
    write(UNIT9,'(a,2I6)')"M and N are:: ", M, N 
    allocate(rankx( 1:(M+1)*(N+1)  ) )
    allocate(ranky( 1:(M+1)*(N+1)  ) )
    allocate(rankxQ( 1:(M+1)*(N+1)  ) )
    allocate(rankyQ( 1:(M+1)*(N+1)  ) )
    rankx(:)=0.d0
    ranky(:)=0.d0
    !! lens is the number of Ks in BZ
    !!~ integer(kind=4), intent(IN)     :: M, N
    !!~ real(Kind=DP),   intent(OUT) :: deltaK
    !!~ integer(kind=4), intent(OUT) :: Lens
    !!~ real(KIND=DP), intent(OUT) :: rankx, ranky
    call wavenumbersII( M, N, deltak, Lens, rankx, ranky )
    write(UNIT9,'(a)')
    write(UNIT9,'(a,1xES23.15E3, a, I6)')"deltak value is 2.0*Ybound/real(N,DP)= ", deltak, "   and Lens=", Lens 
    write(UNIT9,'(a)')   
    write(UNIT9,'(a)') "please have a check rankx(7) and ranky(7) in wavenumber.dat in 9th line"
    write(UNIT9,'(a)')
    write(UNIT9,'(a, 2(1xES23.15E3) )') "rankx(7), ranky(7) are:: ",  rankx(7), ranky(7)
    !=====================================
    !~        ^  N
    !~        |
    !~  ------|----------->   M
    !~        |
    ! wavenumbers in BZ
    !
    !  Kprime (2Pi/3sqrt(3)a ,  2pi/3a  ) ;  K  (4pi/3asqrt3, 0 )
    !                      -----
    !                   /          \
    !                   \          / 
    !                       ----
    !!!          BZ for graphene
    !
    !  peiliang write at Feb. 2018
!============================================================================
!~ ============              Bandstructure0.dat             =================
!============================================================================
    if (SetTransitionMetalDichalcogenides>0 .and. SetTransitionMetalDichalcogenides<5 ) then 
        Idimension=11*Nlayer
    else
        write(UNIT9,'(a)') "Error!  Please check the SetTransitionMetalDichalcogenides!!!"
    end if 
    
    write( UNIT9,'(a)'  )
    write( UNIT9,'(a, I4)' ) "Idimension = ",  Idimension
    
    allocate(    eigX                                    (1:Idimension)    )
    allocate(    vecX                   ( 1:Idimension,  1:Idimension )    )
    allocate(    ckeig                  ( 1:Idimension,        1:Lens )    )
    allocate(    ckvec             ( Idimension , Idimension , 1:Lens )    )
    allocate(    prob                   ( Idimension,      1:3*Nlayer )    )

    layers_tipsi=Nlayer
     !!===  get the hopping values ts 	   < TshoppingIV_SMO.f90 >
    call Hopping_amplitudes (SetTransitionMetalDichalcogenides, layers_tipsi,  tdp, tdd, tpp, Delta0, Delta1, Delta2, DeltaP, DeltaZ ) 
    !!~ ~ ~~~~~~ compare Hopping_Amplitudes ~~~~~~~~~~~~~~~~~~~
     !~ call Hoppingvalues(SetTransitionMetalDichalcogenides, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ)
     !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	 
    Vonset( 1:3,   1:Nlayer)=0.d0   ! eV 
    write( UNIT9, *  ) "Vonset:", Vonset
    
    call BandKs2(Blengths, Bkx, Bky)   !! Blenghts defined in moduleone  ! wavenumbers.f90
      do i=0, Blengths-1
            call Hamiltonianofxs2VII( Idimension, Bkx(i), Bky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vonset, eigX, vecX)  
            do j=1,Idimension    
                write(UNIT7, '(1x,2(1xE18.11) )') real(i),  eigX(j)      !!!===   Bandstructure0.dat 
            end do
      end do
      !!~ plot Bandstructure0.pdf/eps
    call showBandstructure0plot (SetTransitionMetalDichalcogenides, Bandstructure0_Windows_Y0, Bandstructure0_Windows_Y1)  
    write( UNIT9,  '( a )'   )
    write( UNIT9, '( a, 1xES23.15E3)') 'Bandstructure0.dat has already finished.'
    
    
    !~ =================================== check Hamiltonian_Data ======================================
    !~ Vonset( 3,   1)=  0.3d0
    !~ Vonset( 2,   1)=  0.d0
    !~ Vonset( 1,   1)=  -0.3d0
    !~ print*
    !~ print*, "Vonset =", Vonset
    !~ i=0          ! gamma point
    !~ i= 296       ! K point
    !~ call Hamiltonianofxs2VII_DATA ( Idimension, Bkx(i), Bky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vonset, eigX, vecX)


!============================================================================
!~ ===========     Density of State without doping         ==================
!============================================================================
    DensityofstateTMDC(:)=0.d0
    do i=1, Lens
        call Hamiltonianofxs2VII( Idimension, rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eigX, vecX)
            do j=1,Idimension       
                        prob0=1.0d0
                        call counter( eigX(j), prob0, deltax, Ndeltax+1, DosTMDC )
                        DensityofstateTMDC=DosTMDC+DensityofstateTMDC

                        ckeig( j, i )=eigX(j)
                        ckvec( :, j, i )=vecX(:, j)
            end do 
    end do
    call Rnorm(Ndeltax+1, DensityofstateTMDC, r0)                  !! LIBRARY.f90

    DensityofstateTMDC=DensityofstateTMDC*2.0_DP/r0
    write(  UNIT9,    '( a )'   ) 
    write(  UNIT9,    '( a, 1xES23.15E3 )'   ) "The separation of Density of state [deltax] is ", deltaX  !! valued in main.f90 file

    rdos=0.d0
    do i=0, NE0
        r=DensityofstateTMDC(i)
        rdos=rdos+r
    end do  
    write(  UNIT9,    '( a )'   ) ''
    write(  UNIT9,    '( a, F14.10 )'   )  'rdos is ', rdos

    rdos2=0.d0 
    do i=NE0+1, Ndeltax
            r=DensityofstateTMDC(i)
            rdos2=rdos2+r
    end do
    write(  UNIT9,    '( a )'   ) ''
    write(  UNIT9,    '( a, F14.10 )'   )  'rdos2 is ' ,rdos2
    write(  UNIT9,    '( a )'   ) ''
    write(  UNIT9,    '( a, F14.10 )'   )  'rdos+rods2 should be 2 here = ',  rdos+rdos2

    DensityofstateTMDC=DensityofstateTMDC/deltaX

    write( UNIT200,    '(5x,"Energy(E)" ,6x, "Densityofstate" )' )   !!=== dos0.out
    do i=0,  Ndeltax
        write(  UNIT200,   '(1x,E16.8,3x, E16.8 )' )   AL+real(i)*deltax ,  DensityofstateTMDC(i)       !!=== dos0.out 
    end do   
!============================================================================
!~ ===============      Carrier Density  Input             ==================
!============================================================================

    RNprime=RNprimereal*Sunitcell*10.0d0**(17)   !! R~~ 10^13 cm^-2
    write(*, '( a)') ""
    write(*, '( 1x,a,  f14.8,  a )') "The given carrier density=",  RNprimereal , "*E13/cm^2  "
    write(UNIT9, '( a)') ""
    write(UNIT9, '( a,  f14.8,  a )') "The given carrier density=",  RNprimereal , "*E13/cm^2  "
    write(UNIT9, '( a)') ""
    write(UNIT9, '( a, 1xES23.15E3, a )'  ) "The given carrier density=",  RNprime, ", which is £¨ n in electrons per unit £©."
    write(UNIT9, '( a)') ""
    RNprime=RNprime/real(Idimension) 
    R2K_zhuanhuanguanxi=Sunitcell*10.0d0**(17)/real(Idimension)
    write(UNIT9, '( a)') ""
    write(UNIT9, '( a, 1xES23.15E3 )') "R2K_zhuanhuanguanxi( ~10^13 cm^{-2} to ~ e/unit) =", R2K_zhuanhuanguanxi

!============================================================================
!~ ====================     SELF-CONSISTENT CALCULATION     =================
!============================================================================

    rdosprime0=rdos
    call Chargedistribution( RNprime, Doping_way, Nlayer, RNprimelayer )
    call Kappavalue ( Doping_way, kappa )
    !## open file  'screening_interaction.dat' 206
    open (unit =UNIT206, file = FILENAME206,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror) 
    write(UNIT206, '(a, 4(8xa) )' ) '#-interation',  'layer1',  'layer2',  'layer3', ' Total' 

    do  M=1,50     !! DO LOOP III
        !======================   Loops   ================================ 
        write(UNIT9, '(a, I2, a)') "---------------- LOOP ", M, " ---------------------------"
        
        if  ( Doping_way==1 ) then
            call Vfun( kappa, RNprimelayer(2)+RNprimelayer(3), Idimension, distanceoflayer,  Sunitcell, Vdeta12)  ! same to PRB79,035421 model
            call Vfun( kappa, RNprimelayer(3), Idimension, distanceoflayer, Sunitcell, Vdeta23)
            Vonset( 3,   1)= -Vdeta23  
            Vonset( 2,   1)=  0.d0
            Vonset( 1,   1)= +Vdeta12   
        else if ( Doping_way==2 ) then
            call Vfun( kappa, RNprimelayer(2), Idimension, distanceoflayer, Sunitcell, Vdeta)
            Vdeta = Vdeta/2.0d0
            Vonset( 3,   1)= Vdeta  
            Vonset( 2,   1)=  0.d0
            Vonset( 1,   1)= Vdeta   
        end if
        write( UNIT9,  '( a )'   )
        write( UNIT9,  '( a,  2x,E16.9 )'   ) "V1 = ", Vonset( 1,   1)
        write( UNIT9,  '( a )'   )
        write( UNIT9,  '( a,  2x,E16.9 )'   ) "V3 = ", Vonset( 3,   1)     !!  v2=0 eV 
        
        write( UNIT9,  '( a )'   )
        write( UNIT9,  '( a,  3(2xE16.9) )'  ) "Vonset = ", Vonset( 1:3, 1 )

        DensityofstateTMDC(:)=0.d0
        Densityofstate(:,:)=0.d0   !! Densityofstate( 0:Ndeltax, 1:3Nlayer )
        !~ do i=1, 2
        do i=1, Lens
            call Hamiltonianofxs2VII( Idimension, rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eigX, vecX)
            do j=1, Idimension
                    prob(j, 1)= abs(vecX(1,j))**2+abs(vecX(2,j))**2+abs(vecX(3,j))**2
                    prob(j, 2)= abs(vecX(4,j))**2+abs(vecX(5,j))**2+abs(vecX(6,j))**2+abs(vecX(7,j))**2+abs(vecX(8,j))**2       
                    prob(j, 3)= abs(vecX(9,j))**2+abs(vecX(10,j))**2+abs(vecX(11,j))**2
                    prob0=1.0d0
                    call counter( eigX(j), prob(j, 1), deltax, Ndeltax+1, DosTMDC )
                    Densityofstate(:, 3)=DosTMDC+Densityofstate(:, 3)
                    call counter( eigX(j), prob(j, 2), deltax, Ndeltax+1, DosTMDC )
                    Densityofstate(:, 2)=DosTMDC+Densityofstate(:, 2)
                    call counter( eigX(j), prob(j, 3), deltax, Ndeltax+1, DosTMDC )
                    Densityofstate(:, 1)=DosTMDC+Densityofstate(:, 1)
                    call counter( eigX(j), prob0, deltax, Ndeltax+1, DosTMDC )
                    DensityofstateTMDC=DosTMDC+DensityofstateTMDC
            end do     
        end do
        call Rnorm(Ndeltax+1, DensityofstateTMDC, r0)
        DensityofstateTMDC=DensityofstateTMDC*2.0_DP/(deltaX*r0*rdos)
        Densityofstate(:, 1)=Densityofstate(:, 1)*2.0_DP/(deltaX*r0*rdos) 
        Densityofstate(:, 2)=Densityofstate(:, 2)*2.0_DP/(deltaX*r0*rdos) 
        Densityofstate(:, 3)=Densityofstate(:, 3)*2.0_DP/(deltaX*r0*rdos) 

        Biggest_inter_Value=0.0_DP
        do i=1, Ndeltax
            B=  deltaX*abs( DensityofstateTMDC(i)-DensityofstateTMDC(i-1) )
                if( B>= Biggest_inter_Value ) then
                    Biggest_inter_Value=B
                else        
                    cycle
                end if 
        end do
        Serror=1.1_DP*Biggest_inter_Value
        write( UNIT9,  '(a)'  )
        write( UNIT9,  '(a, 1xES23.15E3)'  ) "The smallest Serror is equal to ",  Serror
        
        call paoguang( Ndeltax+1, DensityofstateTMDC )
        call paoguang( Ndeltax+1, Densityofstate(:, 1) )
        call paoguang( Ndeltax+1, Densityofstate(:, 2) )
        call paoguang( Ndeltax+1, Densityofstate(:, 3) )
        
        Densityofstatecheck=DensityofstateTMDC-Densityofstate(:, 1)- Densityofstate(:, 2)- Densityofstate(:, 3)
        call paoguang( Ndeltax+1, Densityofstatecheck )
        N0=0
        call  area( Ndeltax+1, Densityofstatecheck, N0, Ndeltax, S )
        S=S*deltaX
        write( UNIT9,  '(a)'  )
        write( UNIT9,  '(a, 1xES23.15E3)'  )    'Zero should be here = ', S 
    
        rdosprime0= 1.0d0
        Do kndelta=10, Ndeltax  
            N0=0
            S=0.0_DP 
            call  area( Ndeltax+1, DensityofstateTMDC, N0, kndelta, S )
            S=S*deltaX
            !!~ write(  8,   *  )  'S= ' ,S
        
            if ( (S-rdosprime0) > Serror) then
                write( UNIT9, '( a, F14.10 )'   )  'S= ' ,S
                write( UNIT9, '(  a  )'   )
                write( UNIT9, '(  a  )'   ) 'The program was stopped, please have a check Serror. '
                write( UNIT9, '(  a  )'   )
                stop  
            else if ( abs(S-rdosprime0) .LE. Serror ) then
                s2=0.0_DP
                call  area( Ndeltax+1, DensityofstateTMDC, N0, kndelta+1, s2 )
                s2=s2*deltaX
                
                if( abs(s2-rdosprime0)*1e5*1.d0  .GT. abs(S-rdosprime0)*1e5*1.d0   ) then
                    write( UNIT9, '(  a  )'   )
                    write( UNIT9, '(  a, 1xES23.15E3  )'   ) 'S should be close to one = ' ,S 
                                neutralitypoint=kndelta
                    write( UNIT9, '(  a  )'   )
                    write( UNIT9, '(  a, 6I  )'   )'The neutralitypoint is   ', neutralitypoint
                    
                    call  area( Ndeltax+1, DensityofstateTMDC, N0, neutralitypoint+1, S )
                            S=S*deltaX
                    !!~ write( 8,'( a )')''
                    !!~ write( 8,   *   )  'S(0:neutralitypoint+1)= ' ,S
                    
                    exit
                else 
                    cycle
                end if 
            end if
        end do
        
        Do kndelta=neutralitypoint+3, Ndeltax   
                S=0.0_DP 
                call  area(Ndeltax+1, DensityofstateTMDC, neutralitypoint, kndelta, S )
                S=S*deltaX
    
                if ( S-RNprime > Serror) then
                    write( UNIT9, '( a, F14.10 )'   )  'S= ' ,S
                    write( UNIT9, '(  a  )'   )
                    write( UNIT9, '(  a  )'   )'The program was stopped, please have a check Serror.'
                    write( UNIT9, '(  a  )'   )
                    stop  
                else if ( abs(S-RNprime) .LE. Serror ) then
                    s2=0.0_DP
                    call  area(Ndeltax+1, DensityofstateTMDC, neutralitypoint, kndelta+1, s2 )
                    s2=s2*deltaX
                    if( abs(s2-RNprime)  .GT. abs(S-RNprime)   ) then
                            write( UNIT9, '(  a  )'   )
                            write( UNIT9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu)-Sprime(neutralitypoint)= ", S
                                Kmiu=kndelta
                            write( UNIT9, '(  a  )'   )
                            write( UNIT9, '(  a, 6I  )'   )	 'The Kmiu is   ', Kmiu
                            exit
                    else 
                            cycle
                    end if 
                end if
        end do
        call  area( Ndeltax+1, DensityofstateTMDC, neutralitypoint, Kmiu+1, S )
                S=S*deltaX
        write( Unit9, '(  a  )'   )
        write( Unit9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu+1)-Sprime(neutralitypoint)= ", S
        call  area( Ndeltax+1, DensityofstateTMDC, neutralitypoint, Kmiu-1, S )
                S=S*deltaX
        write( Unit9, '(  a  )'   )
        write( Unit9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu-1)-Sprime(neutralitypoint)= ", S

        EMiu0=AL+real(neutralitypoint)*deltaX
        EMiu= AL+real(Kmiu)*deltaX  
        write (UNIT9, '(a)') ''
        write (UNIT9, '(a, F16.8)') 'EMiu0=  ', EMiu0    
        write (UNIT9, '(a)') ''      
        write (UNIT9, '(a, F16.8)') 'EMiu=  ', EMiu   

        A1=0.0_DP 
        call  area( Ndeltax+1, Densityofstate(:, 1), neutralitypoint, KMiu, A1)
        A1=A1*deltaX
        A2=0.0_DP 
        call  area( Ndeltax+1, Densityofstate(:, 2), neutralitypoint, KMiu, A2)
        A2=A2*deltaX
        A3=0.0_DP 
        call  area( Ndeltax+1, Densityofstate(:, 3), neutralitypoint, KMiu, A3)
        A3=A3*deltaX
        
        Errordelta=Biggest_inter_Value*0.25_DP
        
        If (   ( abs(A1-RNprimelayer(1))<= Errordelta  .and. abs(A2-RNprimelayer(2))<= Errordelta) .and.  (abs(A3-RNprimelayer(3))<= Errordelta)  ) then
            
            RNprimelayer(1)=(RNprimelayer(1)+A1)/2.0d0
            RNprimelayer(2)=(RNprimelayer(2)+A2)/2.0d0
            RNprimelayer(3)=(RNprimelayer(3)+A3)/2.0d0
            precisionofnprime= abs( RNprimelayer(1)+RNprimelayer(2)+RNprimelayer(3)-RNprime )*100.0d0/RNprime                                  ! percent
            write(UNIT9, '(a)')
            write(UNIT9, '(a, 2(1xES23.15E3) )')  'RNprimelayer(3) =  ', RNprimelayer(3), RNprimelayer(3)/R2K_zhuanhuanguanxi
            write(UNIT9, '(a, 2(1xES23.15E3) )')  'RNprimelayer(2) =  ', RNprimelayer(2), RNprimelayer(2)/R2K_zhuanhuanguanxi
            write(UNIT9, '(a, 2(1xES23.15E3) )')  'RNprimelayer(1) =  ', RNprimelayer(1), RNprimelayer(1)/R2K_zhuanhuanguanxi
            write(UNIT9, '(a)')
            write(UNIT9, '(a, 1xES23.15E3)')  'Nprime_new=  ', RNprimelayer(1)+RNprimelayer(2)+RNprimelayer(3)
            write(UNIT9, '(  a  )' )
            write(UNIT9, '(  a, f10.6, a   )' ) 'The precisionofnprime = ', precisionofnprime, '%.'
            write(UNIT9, '(  40("-") , a, 40("-")    )' ) "  Loop End  "
            
            write(UNIT206, '( 2xI3,4x,4ES26.15E3)' ) M,  A1,  A2, A3, A3+A2+A1 
            write(UNIT206, '( 2xI3,4x,4ES26.15E3)' ) M,  RNprimelayer(1),  RNprimelayer(2), RNprimelayer(3),  RNprimelayer(1)+RNprimelayer(2)+RNprimelayer(3) 
            exit 
        else  
            RNprimelayer(1)=(RNprimelayer(1)+A1)/2.0d0
            RNprimelayer(2)=(RNprimelayer(2)+A2)/2.0d0
            RNprimelayer(3)=(RNprimelayer(3)+A3)/2.0d0 
            write(UNIT9, '(a)')""
            write(UNIT9, '(a,  1xES23.15E3)')  'New_Nprime3=  ', RNprimelayer(3)
            write(UNIT9, '(a,  1xES23.15E3)')  'New_Nprime2=  ', RNprimelayer(2)
            write(UNIT9, '(a,  1xES23.15E3)')  'New_Nprime1=  ', RNprimelayer(1)            
            write(UNIT9, '(a)')""
            write(UNIT9, '(a)') 'Need to run again with Self-consistent.... '  
    
            write(UNIT206, '( 2xI3,4x,4ES26.15E3)' ) M,  RNprimelayer(1),  RNprimelayer(2), RNprimelayer(3),  RNprimelayer(1)+RNprimelayer(2)+RNprimelayer(3) 
        end if 
    end do   !! DO LOOP III
    close(UNIT206)    ! screening_interaction.dat


!============================================================================
!~ ====================             show gap1 & dos         =================
!============================================================================
    DosTMDC(:)=0.d0
    Mdeltax = Ndeltax+1    
    
    open (unit =UNIT8, file = FILENAME8,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)        ! gap1.dat
    open (unit =UNIT103, file = FILENAME103,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)    ! Sarea.dat
    if  ( Doping_way==1 ) then
        call Vfun( kappa, RNprimelayer(2)+RNprimelayer(3), Idimension, distanceoflayer,  Sunitcell, Vdeta12)  ! same to PRB79,035421 model
        call Vfun( kappa, RNprimelayer(3), Idimension, distanceoflayer, Sunitcell, Vdeta23)
        Vonset( 3,   1)= -Vdeta23   ! eV 
        Vonset( 2,   1)=  0.d0
        Vonset( 1,   1)= +Vdeta12   ! eV 
    else if ( Doping_way==2 ) then
        call Vfun( kappa, RNprimelayer(2), Idimension, distanceoflayer, Sunitcell, Vdeta)
        Vdeta = Vdeta/2.0d0
        Vonset( 3,   1)= Vdeta   ! eV 
        Vonset( 2,   1)=  0.d0
        Vonset( 1,   1)= Vdeta   ! eV 
    end if
    
    print*
    print*, "Vonset =", Vonset
    
    write( UNIT9,  '( a )'   )
    write( UNIT9,  '( a,  3(2x,E16.9) )'   ) "Vonset = ", Vonset( 1:3, 1 )
    
    do i=1, Lens
        call Hamiltonianofxs2VII( Idimension, rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eigX, vecX)
        do j=1, Idimension
                    prob(j, 1)= abs(vecX(1,j))**2+abs(vecX(2,j))**2+abs(vecX(3,j))**2
                    prob(j, 2)= abs(vecX(4,j))**2+abs(vecX(5,j))**2+abs(vecX(6,j))**2+abs(vecX(7,j))**2+abs(vecX(8,j))**2       
                    prob(j, 3)= abs(vecX(9,j))**2+abs(vecX(10,j))**2+abs(vecX(11,j))**2
                    prob0=1.0d0
                    call counter( eigX(j), prob(j, 1), deltax, Mdeltax , DosTMDC(:) )
                    Densityofstate(:, 3)=DosTMDC+Densityofstate(:, 3)
                    call counter( eigX(j), prob(j, 2), deltax, Mdeltax , DosTMDC )
                    Densityofstate(:, 2)=DosTMDC+Densityofstate(:, 2)
                    call counter( eigX(j), prob(j, 3), deltax, Mdeltax , DosTMDC )
                    Densityofstate(:, 1)=DosTMDC+Densityofstate(:, 1)
                    call counter( eigX(j), prob0, deltax, Mdeltax , DosTMDC )
                    DensityofstateTMDC=DosTMDC+DensityofstateTMDC
        end do     
    end do
    call Rnorm(Ndeltax+1, DensityofstateTMDC, r0)
    DensityofstateTMDC=DensityofstateTMDC*2.0_DP/(deltaX*r0*rdos)
    Densityofstate(:, 1)=Densityofstate(:, 1)*2.0_DP/(deltaX*r0*rdos) 
    Densityofstate(:, 2)=Densityofstate(:, 2)*2.0_DP/(deltaX*r0*rdos) 
    Densityofstate(:, 3)=Densityofstate(:, 3)*2.0_DP/(deltaX*r0*rdos) 

    Biggest_inter_Value=0.0_DP
    do i=1, Ndeltax
            B=  deltaX*abs( DensityofstateTMDC(i)-DensityofstateTMDC(i-1) )
                if( B>= Biggest_inter_Value ) then
                    Biggest_inter_Value=B
                else        
                    cycle
                end if 
    end do
    Serror=1.1_DP*Biggest_inter_Value
    write( UNIT9,  '(a)'  )
    write( UNIT9,  '(a, 1xES23.15E3)'  ) "The smallest Serror is equal to ",  Serror

    rdosprime0= 1.0d0
    Do kndelta=10, Ndeltax  
            N0=0
            S=0.0_DP 
            call  area( Ndeltax+1, DensityofstateTMDC, N0, kndelta, S )
            S=S*deltaX
            write( UNIT103, '(a, 1xES23.15E3)'  )  'S= ' ,S
        
            if ( (S-rdosprime0) > Serror) then
                write( UNIT9, '( a, F14.10 )'   )  'S= ' ,S
                write( UNIT9, '(  a  )'   )
                write( UNIT9, '(  a  )'   ) 'The program was stopped, please have a check Serror. '
                write( UNIT9, '(  a  )'   )
                stop  
            else if ( abs(S-rdosprime0) .LE. Serror ) then
                s2=0.0_DP
                call  area( Ndeltax+1, DensityofstateTMDC, N0, kndelta+1, s2 )
                s2=s2*deltaX
                
                if( abs(s2-rdosprime0)*1e5*1.d0  .GT. abs(S-rdosprime0)*1e5*1.d0   ) then
                    write( UNIT9, '(  a  )'   )
                    write( UNIT9, '(  a, 1xES23.15E3  )'   ) 'S should be close to one = ' ,S 
                                neutralitypoint=kndelta
                    write( UNIT9, '(  a  )'   )
                    write( UNIT9, '(  a, 6I  )'   )'The neutralitypoint is   ', neutralitypoint
                    
                    call  area( Ndeltax+1, DensityofstateTMDC, N0, neutralitypoint+1, S )
                            S=S*deltaX
                    write( UNIT103,'( a )')''
                    write( UNIT103,'(a, 1xES23.15E3)' )  'S(0:neutralitypoint+1)= ', S
                    
                    exit
                else 
                    cycle
                end if 
            end if
    end do
    
    Do kndelta=neutralitypoint+3, Ndeltax   
            S=0.0_DP 
            call  area(Ndeltax+1, DensityofstateTMDC, neutralitypoint, kndelta, S )
            S=S*deltaX
    
            if ( S-RNprime > Serror) then
                    write( UNIT9, '( a, F14.10 )'   )  'S= ' ,S
                    write( UNIT9, '(  a  )'   )
                    write( UNIT9, '(  a  )'   )'The program was stopped, please have a check Serror.'
                    write( UNIT9, '(  a  )'   )
                    stop  
            else if ( abs(S-RNprime) .LE. Serror ) then
                    s2=0.0_DP
                    call  area(Ndeltax+1, DensityofstateTMDC, neutralitypoint, kndelta+1, s2 )
                    s2=s2*deltaX
                    if( abs(s2-RNprime)  .GT. abs(S-RNprime)   ) then
                            write( UNIT9, '(  a  )'   )
                            write( UNIT9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu)-Sprime(neutralitypoint)= ", S
                                Kmiu=kndelta
                            write( UNIT9, '(  a  )'   )
                            write( UNIT9, '(  a, 6I  )'   )	 'The Kmiu is   ', Kmiu
                            exit
                    else 
                            cycle
                    end if 
            end if
    end do
    call  area( Ndeltax+1, DensityofstateTMDC, neutralitypoint, Kmiu+1, S )
                S=S*deltaX
    write( Unit9, '(  a  )'   )
    write( Unit9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu+1)-Sprime(neutralitypoint)= ", S
    call  area( Ndeltax+1, DensityofstateTMDC, neutralitypoint, Kmiu-1, S )
                S=S*deltaX
    write( Unit9, '(  a  )'   )
    write( Unit9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu-1)-Sprime(neutralitypoint)= ", S

    EMiu0=AL+real(neutralitypoint)*deltaX
    EMiu= AL+real(Kmiu)*deltaX  
    write (UNIT9, '(a)') ''
    write (UNIT9, '(a, F16.8)') 'EMiu0=  ', EMiu0    
    write (UNIT9, '(a)') ''      
    write (UNIT9, '(a, F16.8)') 'EMiu=  ',  EMiu 

                    !*********   End EMiu0 EMiu   *********
    !## open file  'wavenumber.dat'
    open (unit =UNIT106, file = FILENAME106,  status="UNKNOWN" ,  iostat=ierror)   
    write(UNIT106,'(69("=") )')
    write(UNIT106, '( 12(" "), a, 20(" "), a, 30(" "), a )' )  "Kx", "Ky", "Energy(1:11)"

    do i=0, Blengths-1
        call Hamiltonianofxs2VII( Idimension, Bkx(i), Bky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vonset, eigX, vecX)  
        do j=1,Idimension    
                write(UNIT8, '(1x,2(1xE18.11) )') real(i),  eigX(j)      !!!    !gap1.dat
        end do
    end do    

    call showgap1plot ( SetTransitionMetalDichalcogenides, Gap1_Windows_Y0, Gap1_Windows_Y1, RNprimereal, EMiu0, EMiu  )   ! gap1.dat.pdf
    close(UNIT8)    ! gap1.dat
    close(UNIT103)  ! Sarea.dat
    !*********   End gap1   *********

    do i=1, Lens
        call Hamiltonianofxs2VII( Idimension, rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eigX, vecX)
        do j=1, Idimension
                    prob(j, 1)= abs(vecX(1,j))**2+abs(vecX(2,j))**2+abs(vecX(3,j))**2
                    prob(j, 2)= abs(vecX(4,j))**2+abs(vecX(5,j))**2+abs(vecX(6,j))**2+abs(vecX(7,j))**2+abs(vecX(8,j))**2       
                    prob(j, 3)= abs(vecX(9,j))**2+abs(vecX(10,j))**2+abs(vecX(11,j))**2
                    prob0=1.0d0
                    call counter( eigX(j), prob(j, 1), deltax, Ndeltax+1, DosTMDC )
                    Densityofstate(:, 3)=DosTMDC+Densityofstate(:, 3)
                    call counter( eigX(j), prob(j, 2), deltax, Ndeltax+1, DosTMDC )
                    Densityofstate(:, 2)=DosTMDC+Densityofstate(:, 2)
                    call counter( eigX(j), prob(j, 3), deltax, Ndeltax+1, DosTMDC )
                    Densityofstate(:, 1)=DosTMDC+Densityofstate(:, 1)
                    call counter( eigX(j), prob0, deltax, Ndeltax+1, DosTMDC )
                    DensityofstateTMDC=DosTMDC+DensityofstateTMDC
        end do  
        write(UNIT106, '(1x, 2(1xES23.15E3), 11(1xES23.15E3)   )') rankx(i), ranky(i), eigX
    end do
    close( UNIT106 )  ! KxkyE_k.dat
    
    call Rnorm(Ndeltax+1, DensityofstateTMDC, r0)
    DensityofstateTMDC=DensityofstateTMDC*2.0_DP/(deltaX*r0)
    Densityofstate(:, 1)=Densityofstate(:, 1)*2.0_DP/(deltaX*r0) 
    Densityofstate(:, 2)=Densityofstate(:, 2)*2.0_DP/(deltaX*r0) 
    Densityofstate(:, 3)=Densityofstate(:, 3)*2.0_DP/(deltaX*r0) 
    N0=0
    call  area( Ndeltax+1, DensityofstateTMDC, N0, Ndeltax, S )
    S=S*deltaX
    write(  UNIT9,    '( a )'   ) ''
    write(  UNIT9,    '( a, F14.10 )'   )  'The area of Densityofstate should be 2 Here:: ' ,S
    
    open (unit =UNIT5, file = FILENAME5,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  ! dos.dat
    write(UNIT5,'(5x,"Energy(E)" ,6x, "Density of state",2x," Sum of DOS(123)", 4x," Densityofstate1 ", 2x,&
    & " Densityofstate2 ", 2x," Densityofstate3 "   )' )
    do i=0,Ndeltax
        write(UNIT5, '(1x,E16.8,3x, E16.8, 4(3x,E16.8) )' ) AL+real(i)*deltax , DensityofstateTMDC(i), Densityofstate(i, 1) +Densityofstate(i, 2)+Densityofstate(i, 3),  Densityofstate(i, 1), Densityofstate(i, 2), Densityofstate(i, 3)
    end do
    close(UNIT5)
                    !*********   End dos   *********

!============================================================================
!~ ====================        K points at fermi level      =================
!============================================================================

    open (unit =UNIT201, file = FILENAME201,  status="UNKNOWN" , ACTION="READWRITE", iostat=ierror)  !K_ferimilevel.dat  
    open (unit =UNIT12, file = FILENAME12,  status="UNKNOWN" , ACTION="READWRITE", iostat=ierror)  !  'qvectors.dat '  
    open (unit =UNIT15, file = FILENAME15,  status="UNKNOWN" , ACTION="READWRITE", iostat=ierror)  !  'qvectors1.dat '  

    !!~ Fermibound is defined in module1.f90
    write(UNIT201, '( a, f16.8, 4x, a, f16.8 )' ) "  Emiu =", Emiu,  "  Fermibound=", Fermibound  !K_ferimilevel.dat 
    write(UNIT201, '(  a )' ) ''

    KNumberoffermi=0
    allocate(  FermiKS ( 1:Lens , 1:2) )           
    FermiKS(:, :)=0.d0                                       

    NEf=0   !! The number of Efs in the fermi level
    do i=1, Lens
        call Hamiltonianofxs2VII( Idimension, rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eigX, vecX)
                judge=0
        do k=1,Idimension 
                If ( abs(eigX(k)-EMiu) <= Fermibound ) Then    ! Fermibound   defined in module1 file
                        judge=judge+1
                end if
                ckeig( k, i )=eigX(k)
                ckvec( :, k, i )=vecX(:, k)
        end do
        NEf=NEf+judge
        If (Judge.GE. 1) then 
                KNumberoffermi=KNumberoffermi+1
                write (UNIT201, '(1x,2f16.8)')  rankx(i), ranky(i)      !! FILENAME201='K_ferimilevel.dat'
                FermiKS( KNumberoffermi, 1 )=rankx(i)
                FermiKS( KNumberoffermi, 2 )=ranky(i)
        end if  
    end do   
    close( UNIT201 )      !K_ferimilevel.dat  
    call Kferimilevelplot (  )
    write( UNIT9,  '( a )' )
    write( UNIT9,  '( a, I10 )' ) "Nef = ", NEf

    r=0.d0
    do i=1,11
            r=r+abs(ckvec( i, 3, 9))**2
    end do
    write( UNIT9,  '(a)' )
    write( UNIT9,  '(a, 3f16.8 )' ) "The Vec_norm 1 should be here::",  r,  dot_product( ckvec( :, 3, 9),  ckvec( :, 3, 9) ) 
    write( UNIT9,  '(a)' )
    write( UNIT9,  '(a, I6)' ) "The number of the K points at Fermi level at the First time is ",  KNumberoffermi

    !~ ===      QK    ======
    LenQKS=KNumberoffermi*(KNumberoffermi-1)/2
    
    allocate(      abs_q    ( 1: LenQKS   )      )
    allocate(      QKS      ( LenQKS, 1:3 )      )  
    allocate(      oldQKS   ( LenQKS, 1:3 )      ) 
    allocate(      CdielectricFun        (  1: LenQKS )   )
    allocate(      CWcoulomb_interaction (  1: LenQKS )   )
    
    write ( UNIT9,  '(a)' )
    write ( UNIT9,  '( a, I10 )' ) "LenQKS = ", LenQKS
    write ( UNIT12, '(9x,a, 15x,a, 12x, a )') "qx", "qy", " || q || "
    write ( UNIT12, '( a )') 
    write ( UNIT15, '(9x,a, 15x,a, 12x, a )') "qx", "qy", " || q || "
    write ( UNIT15, '( a )') 

    Lqks=0
    do j=KNumberoffermi, 2, -1   
        do i=1,j-1
                Lqks=Lqks+1
                QKS(Lqks, 1:2)=FermiKS(j , :)-FermiKS(i , :)
                QKS(Lqks, 3)=sqrt( QKS(Lqks, 1)**2+QKS(Lqks, 2)**2  )
                write (UNIT12, '(1x,3f16.8)') FermiKS(j , :)-FermiKS(i , :), QKS(Lqks, 3)   !! 'qvectors.dat ' 
        end do 
    end do
    write( UNIT9,  '(a)' )
    write( UNIT9, '( a,  1x,3f16.8  )' ) "Random check QKS(13) in qvector.dat file at 15th line::",  QKS(13, :)

    oldQKS=QKS
    abs_q( 1:LenQKS )=QKS( 1:LenQKS, 3 )
    call QsortC(abs_q)
    iloop: do i=1, LenQKS
        Mloop:  do  M=1, LenQKS
                        if ( oldQKS( M , 3 )==abs_q(i)       ) then
                                    QKS(i, 1:3)=oldQKS(M, 1:3)
                                    oldQKS( M , 3 )=0.d0
                                    exit Mloop
                        end IF
                end do Mloop
                write (UNIT15, '(1x,3f16.8)')  QKS(i, 1:3)   !! 'qvectors1.dat ' 		
           end do  iloop
    
    close( UNIT12 )      !qvectors.dat
    close( UNIT15 )      !qvectors1.dat


If( 2.0 > 3.0 ) then 
    !!~======================    Coulomb_coupling_constant { Miu }  ========================
    !~ If(2>3) then  
    !~ ======      CMiu (q)    =========
    allocate(    ckvecM    ( Idimension, Idimension )   )
    allocate(    ckeigQ    ( Idimension, 1:Lens )       )
    allocate(    ckvecQ    ( Idimension, Idimension, 1:Lens )   )
    allocate(    Fenzi     ( Idimension, Idimension, 1:Lens )   )
    
    open (unit =UNIT19, file = FILENAME19,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  !    'WCoulomb_and_dielectric.dat'
    open (unit =UNIT106, file = FILENAME106,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  !   .dat
    open (unit =UNIT107, file = FILENAME107,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  !   CdielectricFun.dat

    write ( UNIT19,  '(9x,a, 15x,a, 9x, a, 10x, a )') "qx", "qy", "-Re_Dynal ",  "-Im_Dynal "
    write ( UNIT107, '(  3(20XA), 35xA )'  ) "FFFZ", "FCIAFM", "CdielectricFun(n)",  "CWcoulomb_interaction(n)"   !   CdielectricFun.dat


    RdistanceX=0.0d0;     RdistanceY=  -1.0d0    !! Mo-S distance a
    Omega=0.d0
    Fmiu = EMiu
    Deltasmall=0.005d0
    alpha=2.0/(4*pi**2)
    alpha=(8.0d0*pi**2/(3.0d0*sq3*Lens))*alpha
    
    Cmiu= cmplx(0._DP, 0._DP) 
   
    !~ do n=9000, 9160
    do n=1, LenQKS   
        Qx=QKS( n , 1 )
        Qy=QKS( n , 2 )	 
        phase=exp(-AI*(Qx*RdistanceX+Qy*RdistanceY ))   
        rankxQ=rankx+Qx
        rankyQ=ranky+Qy	 
        do i=1, Lens
        call Hamiltonianofxs2VII( Idimension, rankxQ(i), rankyQ(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ, Vonset, eigX, vecX )
            do j=1,Idimension       
                ckeigQ( j, i )=eigX(j)
                ckvecQ( :, j, i )=vecX(:, j)
            end do    
        end do
        ckvecQ( 4:8 , :, :)=ckvecQ( 4:8 , :, :)*phase
  
        DynamicalpolarizationPI=cmplx(0._DP, 0._DP)   
        DO Lk=1, Lens  
            do j=1, Idimension
                    Do i=1, Idimension
                        fermiIJ=FermiDistribution ( ckeig(i,Lk), Fmiu, BoltzmannBETA )-FermiDistribution( ckeigQ(j,Lk), Fmiu, BoltzmannBETA  )
                        complexcoeff= dot_product(ckvec( : , i, Lk) , ckvecQ( : , j, LK))     !!! ckvec( : , i, Lk)   should be at around  640 line
                        corelationeff=abs(complexcoeff)**2
                        Fenzi(i, j, Lk)=corelationeff*fermiIJ
                        DynamicalpolarizationPI=DynamicalpolarizationPI+ Fenzi(i, j, Lk)/(AI*Deltasmall+ckeig(i,Lk)-ckeigQ(j,Lk)) 
                    end do
            end do
        end do   ! LK
        DynamicalpolarizationPI= alpha*DynamicalpolarizationPI
        write(UNIT19, '( 1x, 2f16.8,   9xe20.10,  7xe20.10  )' )  QKS( n , 1:2 ),   -DynamicalpolarizationPI
        
        FCIAFM=FouriercomponentCoulombinteractionFME ( Atmdc, kappa, Qx,   Qy )
        CdielectricFun(n)=1.0d0- FCIAFM*DynamicalpolarizationPI
        
        

        FFFZ=FouriercomponentCoulombinteractionFZ ( Atmdc, Qx, Qy )
        !!~ FFFZ=FouriercomponentCoulombinteractionFZ ( Atmdc, kappa, Qx,   Qy )
        
        CWcoulomb_interaction(n)=FFFZ/( CdielectricFun(n) )
        
        write( Unit107, '( 6(2xES26.15E3) )' ) FFFZ, FCIAFM,  CdielectricFun(n), CWcoulomb_interaction(n)
        
        Cmiu = Cmiu + CWcoulomb_interaction(n)
    end do     !! LenQKS
    
    Cmiu = Cmiu/real(LenQKS) 
    Cmiu = Cmiu/real(NEf) 
    !~ Cmiu = Cmiu/DensityofstateTMDC( Kmiu ) 
    write( UNIT9,  '(a)' )
    write( UNIT9,  '( A, 2ES26.15E3 )' ) "(Complex number) Cmiu = ",  Cmiu
    close( UNIT19 )    !    'WCoulomb_and_dielectric.dat'
    close( UNIT107 )    !    CdielectricFun.dat
    close( Unit106 )
    !~ end if
    
    !~ open ( unit= UNIT104, file= FILENAME104,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror )  !    'fermiIJ.dat'
    !~ write (UNIT19, '(9x,a, 15x,a, 15x, a, 15x, a, 14x, a )') "qx", "qy",  "|| q || ",  "Re_W ",  "Im_W "   
        !~ write( UNIT106, '( I8, 2ES26.15E3 )'  ) n, CdielectricFun(n)
        !~ write(UNIT19, '( 1x, 2f16.8,   2(3xes20.10)  )' )  QKS( n ,1:2 ),  CWcoulomb_interaction(n)


end if
!=====================================================================
!=======================      Jie wei       ==========================
!=====================================================================
    write(UNIT9,  '(a)' )
    write(UNIT9, '( 79("=")  )') 
    call timestamp (UNIT9)
    PRINT*
    call timestamp (6)
    write(*, '( 20("=")  )') 
    !!~ ------------------------------------------------------------------
    close(UNIT7)       ! Bandstructure0.dat
    close(UNIT200 )   ! dos0.dat
    close(UNIT9)       ! information.out
!~ ----------------------------------------------------------------------
                  !~  Windows SYSTEM
    
    call system ("copy BAT\checkexist.bat ")
    call system ("copy BAT\move.bat ")
    call system ("checkexist.bat "//test_name)
    call system ("move.bat " // test_name)
    call system ('del *.bat')
    
!~ ----------------------------------------------------------------------
                  !~ LINUX SYSTEM

    !~ call system ("cp BAT/checkexist.sh ./ ")
    !~ call system ("cp BAT/move.sh ./ ")
    !~ call system ("./checkexist.sh "//test_name)
    !~ call system ("./move.sh " // test_name)
!~ ----------------------------------------------------------------------
end program mainTMD
!~ ===============================================
!~ ==============    END PROGRAM   ===============
!~ ===============================================




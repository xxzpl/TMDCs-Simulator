MODULE moduleone
implicit none
integer, parameter :: DP=SELECTED_REAL_KIND(15,307)
INTEGER (kind=4) :: IERROR
complex (kind=DP):: AI=(0._DP,1.0_DP)
complex (kind=DP):: Re=(1.0_DP,0._DP)
character(len=1024):: line,word,word2, lines

real(kind=DP), parameter :: pi=3.141592653589793238460_DP,  sq3=dsqrt(3.0_DP)
! Boltzmann constant =8.6173303(50)×10-5  ev/K
!!!~ real(kind=DP) :: kappa11 = ( (3.90d0+1.0d0)/2.0d0 )  !  Kappa*epsilong  unit  F/m or A^2s^4*Kg^{-1}m^{-3}

real(kind=DP), parameter :: Boltzmannconstant=8.617330350*10.0**(-5),   temperatureKelvin=300.0d0 

real(kind=DP) :: BoltzmannBETA=1.d0/(Boltzmannconstant* temperatureKelvin)

real(kind=DP) :: hbar=6.58211951440*10.0**(-16)   !  unit  eV.s
real(kind=DP) :: Planck_constant=4.135667662*10.0**(-15)   !  unit  eV.s
real(kind=DP) :: speedoflight=2.99792458*10.0**(8)   !  unit  m/s

real(kind=DP) :: electron=1.60217653*10.0**(-19) !!   unit ( C )
real(kind=DP) :: epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V   
real(kind=DP) :: Fine_structure_constant = 1.0d0/(137.035d0)  !! jingxi jiegou changshu


integer(KIND=4), parameter::  Blengths=514      !!514            86          !!  840

real(kind=DP), parameter :: carbondistance=1.0_DP                            ! 1.420_DP*(1.0_DP/10**10)
real(kind=DP), parameter :: avalue= carbondistance*dsqrt(3.0_DP)     ! 2.46_DP 10^{-10}
real(kind=DP), parameter ::     AU= 20.0_DP,   AL=-20.0_DP

!~ real(kind=DP), parameter ::     Fermibound=1.5_DP*deltaX          !0.01_DP
real(kind=DP), parameter ::     Fermibound=0.001_DP


real(kind=DP), parameter :: Vec1X=0.5_DP*avalue,   Vec1Y=0.5_DP*sq3*avalue
real(kind=DP), parameter :: Vec2X=0.5_DP*avalue,   Vec2Y=-0.5_DP*sq3*avalue
  
character (len=*), parameter:: material_nameMoS2="MoS2", material_nameWS2="WS2"
character (len=*), parameter:: material_nameMoSe2="MoSe2", material_nameWSe2="WSe2"
character (len=*), parameter:: material_nameGraphene="Graphene"
character (len=*), parameter:: material_nameTriGraphene="TriGraphene"

!~ ================   A  ====================
real(kind=DP), parameter :: AMos2=3.160*10.0**(-10) 
real(kind=DP), parameter :: AWs2=3.153*10.0**(-10) 
real(kind=DP), parameter :: AMoSe2=3.288*10.0**(-10) 
real(kind=DP), parameter :: AWSe2=3.260*10.0**(-10) 

real(kind=DP), parameter ::  Sunitcellmos2=0.50d0*sq3*AMos2**2
real(kind=DP), parameter ::  SunitcellWs2=0.50d0*sq3*AWs2**2      
real(kind=DP), parameter ::  Sunitcellgraphene=2.619*2.0*10.0**(-20)
real(kind=DP), parameter ::  SunitcellmoSe2=0.50d0*sq3*AMoSe2**2            
real(kind=DP), parameter ::  SunitcellWSe2=0.50d0*sq3*AWSe2**2             

!~ real(kind=DP), parameter ::  Sunitcellmos2=8.64778*10.0**(-20)              !!   S=SQRT(3)*a^2/2   a=3.16A
!~ real(kind=DP), parameter ::  SunitcellWs2=4.304756*2.0*10.0**(-20)       !!  a=3.153A
!~ real(kind=DP), parameter ::  Sunitcellgraphene=2.619*2.0*10.0**(-20)
!~ real(kind=DP), parameter ::  SunitcellmoSe2=9.36255*10.0**(-20)            !!   S=SQRT(3)*a^2/2   a=3.288A
!~ real(kind=DP), parameter ::  SunitcellWSe2=9.20377*10.0**(-20)              !!  a=3.260A

!~ ================   u  ====================
real(kind=DP), parameter :: distanceoflayermos2=1.586*10.0**(-10) 
real(kind=DP), parameter :: distanceoflayerWs2=1.571*10.0**(-10)
real(kind=DP), parameter :: distanceoflayerGraphene=3.350*10.0**(-10)
real(kind=DP), parameter :: distanceoflayermoSe2=1.664*10.0**(-10) 
real(kind=DP), parameter :: distanceoflayerWSe2=1.657*10.0**(-10)

!~ ================   W  ====================
real(kind=DP), parameter :: distanceoflayermos2W=2.968*10.0**(-10) 
real(kind=DP), parameter :: distanceoflayerws2W=3.018*10.0**(-10) 
real(kind=DP), parameter :: distanceoflayermoSe2W=3.123*10.0**(-10) 
real(kind=DP), parameter :: distanceoflayerwSe2W=3.108*10.0**(-10) 

!~ ================   Cprime  ====================
real(kind=DP), parameter :: CprimeMoS2=6.140*10.0**(-10) 
real(kind=DP), parameter :: CprimeWS2=6.160*10.0**(-10)
real(kind=DP), parameter :: CprimeMoSe2=6.451*10.0**(-10)
real(kind=DP), parameter :: CprimeWSe2=6.422*10.0**(-10) 

!==============================================================
INTEGER (kind=4) :: UNIT1=21
INTEGER (kind=4) :: UNIT2=22
INTEGER (kind=4) :: UNIT3=23
INTEGER (kind=4) :: UNIT4=24
INTEGER (kind=4) :: UNIT5=25
INTEGER (kind=4) :: UNIT6=26
INTEGER (kind=4) :: UNIT7=27
INTEGER (kind=4) :: UNIT8=28
INTEGER (kind=4) :: UNIT9=29
CHARACTER (LEN=56) :: FILENAME1='HoppingAmplitudes.py'
CHARACTER (LEN=56) :: FILENAME2='eigenvalues.dat'
CHARACTER (LEN=56) :: FILENAME3='CHECKHoppings.dat'
CHARACTER (LEN=56) :: FILENAME4='Hamitonilon_Re.dat'
CHARACTER (LEN=56) :: FILENAME5='dos.dat'
CHARACTER (LEN=56) :: FILENAME6='BandKs.dat'
CHARACTER (LEN=56) :: FILENAME7='Bandstructure0.dat'
CHARACTER (LEN=56) :: FILENAME8='gap1.dat'
CHARACTER (LEN=56) :: FILENAME9='informations.out'


INTEGER (kind=4) :: UNIT10=30
INTEGER (kind=4) :: UNIT11=31
INTEGER (kind=4) :: UNIT12=32
INTEGER (kind=4) :: UNIT13=33
INTEGER (kind=4) :: UNIT14=34
INTEGER (kind=4) :: UNIT15=35
INTEGER (kind=4) :: UNIT16=36
INTEGER (kind=4) :: UNIT17=37
INTEGER (kind=4) :: UNIT18=38
INTEGER (kind=4) :: UNIT19=39
INTEGER (kind=4) :: UNIT20=40
CHARACTER (LEN=56) :: FILENAME10='input.txt'
CHARACTER (LEN=56) :: FILENAME11='wavenumber.dat'
CHARACTER (LEN=56) :: FILENAME12='qvectors.dat'
CHARACTER (LEN=56) :: FILENAME13='showgap1.plt'
CHARACTER (LEN=56) :: FILENAME14='showBandstructure0.plt'
CHARACTER (LEN=56) :: FILENAME15='qvectors1.dat'
CHARACTER (LEN=56) :: FILENAME16='showdos0.plt '
CHARACTER (LEN=56) :: FILENAME17='Hamitonilon_AIMAG.dat'
CHARACTER (LEN=56) :: FILENAME18='DynamicPolyPI6VIII.dat'
CHARACTER (LEN=56) :: FILENAME19='WCoulomb_and_dielectric.dat'


INTEGER (kind=4) :: UNIT100=100
INTEGER (kind=4) :: UNIT101=101
INTEGER (kind=4) :: UNIT102=102
INTEGER (kind=4) :: UNIT103=103
INTEGER (kind=4) :: UNIT104=104
INTEGER (kind=4) :: UNIT105=105
INTEGER (kind=4) :: UNIT106=106
INTEGER (kind=4) :: UNIT107=107
INTEGER (kind=4) :: UNIT108=108
INTEGER (kind=4) :: UNIT109=109
CHARACTER (LEN=56) :: FILENAME100='checkexits.bat'
CHARACTER (LEN=56) :: FILENAME101='checkexits '
CHARACTER (LEN=56) :: FILENAME102='move.bat'
CHARACTER (LEN=56) :: FILENAME103='Sarea.dat '
CHARACTER (LEN=56) :: FILENAME104='fermiIJ.dat'
CHARACTER (LEN=56) :: FILENAME105=' '
CHARACTER (LEN=56) :: FILENAME106='KxkyE_k.dat '
CHARACTER (LEN=56) :: FILENAME107='CdielectricFun.dat'
CHARACTER (LEN=56) :: FILENAME108=' '
CHARACTER (LEN=56) :: FILENAME109=' '



INTEGER (kind=4) :: UNIT200=200
INTEGER (kind=4) :: UNIT201=201
INTEGER (kind=4) :: UNIT202=202
INTEGER (kind=4) :: UNIT203=203
INTEGER (kind=4) :: UNIT204=204
INTEGER (kind=4) :: UNIT205=205
INTEGER (kind=4) :: UNIT206=206
INTEGER (kind=4) :: UNIT207=207
INTEGER (kind=4) :: UNIT208=208
INTEGER (kind=4) :: UNIT209=209
CHARACTER (LEN=56) :: FILENAME200='dos0.dat'
CHARACTER (LEN=56) :: FILENAME201='K_ferimilevel.dat'
CHARACTER (LEN=56) :: FILENAME202='K_ferimilevelplot.plt'
CHARACTER (LEN=56) :: FILENAME203='K_ferimilevel_1.dat'
CHARACTER (LEN=56) :: FILENAME204='K_ferimilevel_2.dat'
CHARACTER (LEN=56) :: FILENAME205='dosU_D.dat'
CHARACTER (LEN=56) :: FILENAME206='screening_iteraction.dat'
CHARACTER (LEN=56) :: FILENAME207=' '
CHARACTER (LEN=56) :: FILENAME208=' '
CHARACTER (LEN=56) :: FILENAME209=' '



END MODULE moduleone
!==============================================================================



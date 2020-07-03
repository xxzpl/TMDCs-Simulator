subroutine read_input_txt_file (    SetTransitionMetalDichalcogenides,  NKnumberBZ, Ndeltax,  Nlayer,  Bandstructure0_Windows_Y0, &
    & Bandstructure0_Windows_Y1, Gap1_Windows_Y0, Gap1_Windows_Y1, RNprimereal,  Eerror, kappa, Icalculatepolarizationfunction,  &
    & deltaX,  neNumberofGraphene, neNumberofMoS2,  neNumberofWS2,  neNumberofMoSe2, neNumberofWSe2, Doping_way, NE0,  EMiu   )
    
use moduleone
use moduletwo
Implicit real (kind=DP) (a-h,o-z) 
integer(KINd=4)::  SetTransitionMetalDichalcogenides, NKnumberBZ, Ndeltax, Nlayer, Icalculatepolarizationfunction, Doping_way

real(kind=DP):: Bandstructure0_Windows_Y0,  Bandstructure0_Windows_Y1, Gap1_Windows_Y0, Gap1_Windows_Y1, RNprimereal,  Eerror, kappa

!~ ==============================================
!~ =============          Read input.txt        ================
!~ ==============================================

open (unit =UNIT10, file = FILENAME10,  status="OLD" , ACTION="READ", iostat=ios)            !!===  input.txt
    if(ios.ne.0) then
        stop ' I/O error, where is input file for reading'
     endif
    ! input-data-file interpreter
    ios=0
    do i=1, 20                                      !//Skip the first 20 lines
        read (  UNIT10, '( a )', iostat=ios  )
    end do
    
    read (  UNIT10, '( a )', iostat=ios  )  line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='SetTransitionMetalDichalcogenides='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (SetTransitionMetalDichalcogenides=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) SetTransitionMetalDichalcogenides
        endif
        endif
    end if
    
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='NKnumberBZ='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (NKnumberBZ=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) NKnumberBZ
        endif
        endif
    end if
    read (  UNIT10, '( a )', iostat=ios  )  line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='Ndeltax='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Ndeltax=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Ndeltax
        endif
        endif
    end if
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='Nlayer='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Nlayer=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Nlayer
        endif
        endif
    end if
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='Bandstructure0_Windows_Y0='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Bandstructure0_Windows_Y0=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Bandstructure0_Windows_Y0
        endif
        endif
    end if
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='Bandstructure0_Windows_Y1='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Bandstructure0_Windows_Y1=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Bandstructure0_Windows_Y1
        endif
        endif
    end if
    
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='Gap1_Windows_Y0='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Gap1_Windows_Y0=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Gap1_Windows_Y0
        endif
        endif
    end if
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then 
        word='Gap1_Windows_Y1='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Gap1_Windows_Y1=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Gap1_Windows_Y1
        endif
        endif
    end if
    
    
    
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then ! skip comment
        word='RNprime(unit/~E13.cm^{-2})='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (RNprime(unit/~E13.cm^{-2})=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) RNprimereal
        endif
        endif
        endif
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then ! skip comment
        word='Eerror of Precision(*deltax)='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Eerror of Precision(*deltax)=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Eerror
        endif
        endif
        endif
    
    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then ! skip comment
        word='Doping_way='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Doping_way=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Doping_way
        endif
        endif
        endif



    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then ! skip comment
        word='kappa='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (kappa=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) kappa
        endif
        endif
        endif



    read(UNIT10,'(a)',iostat=ios) line
    if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then ! skip comment
        word='Icalculatepolarizationfunction='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Icalculatepolarizationfunction=)'
        endif
    word=adjustl(word)
    read(word,*,iostat=ios) Icalculatepolarizationfunction
        endif
        endif
        endif
        
    close(UNIT10)       !! input.txt

!~ print*,   SetTransitionMetalDichalcogenides
!~ print*,   NKnumberBZ
!~ print*,   Ndeltax
!~ print*,   Nlayer
!~ print*,   Bandstructure0_Windows_Y0
!~ print*,   Bandstructure0_Windows_Y1
!~ print*,   RNprimereal
!~ print*,   Eerror
!~ print*,   Icalculatepolarizationfunction


    deltaX=(AU-AL)/( real(Ndeltax)*1.0d0)     ! module1.f90   for FENBIANLV of DOS dependent on Ndeltax
    neNumberofGraphene=2*Nlayer
    neNumberofMoS2=13*Nlayer
    neNumberofWS2=12*Nlayer	
    neNumberofMoSe2=13*Nlayer
    neNumberofWSe2=12*Nlayer

    NE0=Ndeltax/2

    EMiu=0.0_DP     
    Eerror=Eerror*deltaX 

    write( UNIT9,  '( a )'   )
    write( UNIT9, '( a, 1xES23.15E3)') 'Eerror of Precision=', Eerror
    write( UNIT9, '( a )' )
    write( UNIT9, '( a, I2 )' ) "The Charge Doping way is ", Doping_way



end subroutine read_input_txt_file



!~ =================================================================================================

!~ =================================================================================================



subroutine Making_data_folder ( SetTransitionMetalDichalcogenides, neNumberofWSe2,  neNumberofMoSe2, &
& neNumberofWS2,  neNumberofMoS2, neNumberofGraphene,   Nlayer,  RNprimereal, material_name, &
& numberofcarrier, Atmdc, Sunitcell,  distanceoflayer, test_name   )

use moduleone
use moduletwo
Implicit real (kind=DP) (a-h,o-z) 
character (len=50) :: material_name
integer(kind=4), intent(in) ::  SetTransitionMetalDichalcogenides
integer(kind=4), intent(in) :: 	neNumberofWSe2,  neNumberofMoSe2
integer(kind=4), intent(in) :: 	neNumberofWS2,  neNumberofMoS2, neNumberofGraphene

character (len=50), intent(out)  :: test_name
real(kind=DP), intent(out) :: numberofcarrier, Sunitcell, distanceoflayer


    IF (SetTransitionMetalDichalcogenides==0 ) then
        material_name=material_nameGraphene
        L=8    ! lens of the material name
        write(UNIT9, '(a)')
        write(UNIT9, '(a)')"The object of material is Graphene. "
        print*
        write(*, '(1x,a)')"The object of material is Graphene. "
    else if (SetTransitionMetalDichalcogenides==1 ) then
        Atmdc=AMos2
        numberofcarrier=neNumberofMoS2
        Sunitcell=SunitcellmoS2               !! from module1.f90
        distanceoflayer=distanceoflayermoS2
        material_name=material_nameMoS2
        L=4    ! lens of the material name
        write(UNIT9, '(a)')
        write(UNIT9, '(a)')"The object of material is MoS2. "
        print*
        write(*, '(1x,a)')"The object of material is MoS2. "
    else if (SetTransitionMetalDichalcogenides==2 ) then
        Atmdc=AWs2
        numberofcarrier=neNumberofWS2
        Sunitcell=SunitcellWS2
        distanceoflayer=distanceoflayerWS2
        material_name=material_nameWS2
        L=3    ! lens of the material name
        write(UNIT9, '(a)')
        write(UNIT9, '(a)')"The object of material is WS2. "
        print*
        write(*, '(1x,a)')"The object of material is WS2. "
    else if (SetTransitionMetalDichalcogenides==3 ) then
        Atmdc=AMoSe2
        numberofcarrier=neNumberofMoSe2
        Sunitcell=SunitcellMoSe2
        distanceoflayer=distanceoflayerMoSe2
        material_name=material_nameMoSe2
        L=5    ! lens of the material name
        write(UNIT9, '(a)')
        write(UNIT9, '(a)')"The object of material is MoSe2. "
        print*
        write(*, '(1x,a)')"The object of material is MoSe2. "
    else if (SetTransitionMetalDichalcogenides==4 ) then
        Atmdc=AWSe2
        numberofcarrier=neNumberofWSe2
        Sunitcell=SunitcellWSe2
        distanceoflayer=distanceoflayerWSe2
        material_name=material_nameWSe2
        L=4    ! lens of the material name
        write(UNIT9, '(a)')
        write(UNIT9, '(a)')"The object of material is WSe2. "
        print*
        write(*, '(1x,a)')"The object of material is WSe2. "
    
        else
            write(*, '(a)')
            write(*, '(a)') "Error!  What the name of Material? Please have a check!!!"	
            stop
        
    !~ else if (SetTransitionMetalDichalcogenides==5 ) then
        !~ material_name=material_nameGraphene
        !~ Lname=11    ! lens of the material name
        !~ write(UNIT9, '(a)')
        !~ write(UNIT9, '(a)')"The object of material is TriGraphene. "	
        !~ write(*, '(1x,a)')"The object of material is TriGraphene. "	
    end if
    material_name=adjustl(material_name)
    material_name=trim(material_name)

    !! mkdir a folder
    If (  Nlayer<10  ) then
        if( RNprimereal<0.00002d0) then
        write(test_name,  '( a , "Nprime0.00_N", I1 )' ) material_name(1:L),  Nlayer
        end if
        NYU=(RNprimereal-int(RNprimereal))*100 
        If ( int(RNprimereal) <10  .and.   NYU<10   ) then
        write(test_name,  '( a , "Nprime", I1,".",I1, a, "_N",I1 )' ) material_name(1:L),  int(RNprimereal), NYU, "0", Nlayer
        else if( int(RNprimereal) >=10 .and. NYU<10     ) then
        write(test_name,  '( a , "Nprime", I2,".",I1, a, "_N",I1 )' ) material_name(1:L),  int(RNprimereal), NYU, "0", Nlayer
        else if( int(RNprimereal) >=10.0  .and. NYU>=10  ) then
        write(test_name,  '( a , "Nprime", I2,".",I2, "_N",I1 )' ) material_name(1:L),  int(RNprimereal), NYU, Nlayer
        else if( int(RNprimereal) <10.0  .and. NYU>=10  ) then
        write(test_name,  '( a , "Nprime", I1,".",I2, "_N",I1 )' ) material_name(1:L),  int(RNprimereal), NYU, Nlayer
        end if
    else if  (  Nlayer>=10  ) then
        if( RNprimereal<0.00002d0) then
        write(test_name,  '( a , "Nprime0.00_N", I2 )' ) material_name(1:L),  Nlayer
        end if
        NYU=(RNprimereal-int(RNprimereal))*100 
        If ( int(RNprimereal) <10  .and.   NYU<10   ) then
        write(test_name,  '( a , "Nprime", I1,".",I1, a, "_N",I2 )' ) material_name(1:L),  int(RNprimereal), NYU, "0", Nlayer
        else if( int(RNprimereal) >=10 .and. NYU<10     ) then
        write(test_name,  '( a , "Nprime", I2,".",I1, a, "_N",I2 )' ) material_name(1:L),  int(RNprimereal), NYU, "0", Nlayer
        else if( int(RNprimereal) >=10.0  .and. NYU>=10  ) then
        write(test_name,  '( a , "Nprime", I2,".",I2, "_N",I2 )' ) material_name(1:L),  int(RNprimereal), NYU, Nlayer
        else if( int(RNprimereal) <10.0  .and. NYU>=10  ) then
        write(test_name,  '( a , "Nprime", I1,".",I2, "_N",I2 )' ) material_name(1:L),  int(RNprimereal), NYU, Nlayer
        end if
    end if


!~ write(*, '(1x,a)' ) ''
!~ write(*, '( 1x,a, a )') "test_name is ", test_name

write(UNIT9, '( a)') ""
write(UNIT9, '( a,  ES23.15E2  )') "Atmdc(nm) = ", Atmdc*10.0**(10)

write(UNIT9, '(a)' ) ''
write(UNIT9, '(a, a )') "test_name is ", test_name

write(UNIT9, '( a)') ""
write(UNIT9, '( a,  ES23.15E2  )') "The area of Unit Cell (/m**2) =",  Sunitcell

write(UNIT9, '( a)') ""
write(UNIT9, '( a,  ES23.15E2  )') "The distance of interlayers (/m) =",  distanceoflayer


end subroutine Making_data_folder

!~ ============================================================================

!~ ============================================================================




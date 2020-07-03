    program readata
    implicit none
    integer i,j,ierr, k
    character(len=128) str
    real*8,allocatable::data1(:,:)

    open(11,file='gap1.dat')


    i=0
    do               !//Obtain the rest lines of the file 
        read(11,*,iostat=ierr) str
        if(ierr/=0) exit
        i=i+1
    end do
    rewind(11)
    allocate(data1(2,i))

    do j=1,i !//Read the data from test.txt 
        read(11,*) data1(:,j)
    end do

    write(*, *) 'i=', i
    write(*,'(1x,2(1xE18.11) )')  data1(:,4)


    open(12,file='gap1_ZPL.dat')
    do j=8,i, 11 !// first Conductivity-band
        write(12,'(1x,2(1xE18.11) )')  data1(:,j)
    end do
    
    do j=9,i, 11 !// second Conductivity-band
        write(12,'(1x,2(1xE18.11) )')  data1(:,j)
    end do    

    do j=7,i, 11 !// Valens-band
        write(12,'(1x,2(1xE18.11) )')  data1(:,j)
    end do   
    

    close(12)
    close(11)

    end program readata
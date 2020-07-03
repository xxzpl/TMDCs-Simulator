!--------------------------------------------------------
!                             checkexit.bat
!--------------------------------------------------------
subroutine checkexitfolderwindows(  )
use moduleone

open(unit =UNIT100, file = FILENAME100, status="unknown", action="write" , DELIM='none', iostat=ierror)

write(UNIT100,'(a)')  " @echo off " 
write(UNIT100,'(a)')  " set g=%1%   "
write(UNIT100,'(a)')  " if exist %g%\NUL (    " 
write(UNIT100,'(a)')  " RMDIR /S /Q %g% " 
write(UNIT100,'(a)')  " )   " 
write(UNIT100,'(a)')  " @mkdir %g%  " 
write(UNIT100,'(a)')  "  " 
close(UNIT100)
end subroutine


!--------------------------------------------------------
!                             copy
!--------------------------------------------------------
subroutine copytools(  )
character (len=5) :: Foldername="Tools"
call system ("cd  "// Foldername )
call system ("copy *.bat  .. " )
call system ("cd .. " )

end subroutine



subroutine move_files_windows( test_name   )
use moduleone
character (len=50) :: test_name

open(unit =UNIT102, file = FILENAME102, status="unknown", action="write" , DELIM='none', iostat=ierror)

write(UNIT102,'(a)')  " @echo off " 
write(UNIT102,'(a)')  " @del *.obj   "
write(UNIT102,'(a)')  " @del *.mod   " 
write(UNIT102,'(a)')  " set g=%1 " 
write(UNIT102,'(a)')  " @move *.plt %g%   " 
write(UNIT102,'(a)')  " @move *.out %g%  " 
write(UNIT102,'(a)')  " @move *.pdf %g%" 
write(UNIT102,'(a)')  " @move *.dat %g%" 
write(UNIT102,'(a)')  " @move *.eps %g%" 
write(UNIT102,'(a)')  " @copy *.txt %g%" 
write(UNIT102,'(a)')  " @copy *.bat %g%" 
close(UNIT102)

call system ("move.bat " // test_name)
!~ call system ("delete move.bat"   )
!~ call system ( 'delete *.bat'   )

end subroutine




!--------------------------------------------------------
!                             Kferimilevelplot
!--------------------------------------------------------
subroutine Kferimilevelplot (  )
use moduleone
!============  start   K_ferimilevelplot.plt'=======

open(unit =UNIT202, file = FILENAME202, status="unknown", action="write" , DELIM='none', iostat=ierror)

write(UNIT202,'(a)')  " set terminal pdf " 
write(UNIT202,'(a)')  " set output 'Ks_points0.pdf'   "
write(UNIT202,'(a)')  " set border linewidth 1.5    " 
write(UNIT202,'(a)')  " set pointsize 0.20  " 
write(UNIT202,'(a)')  " set style line 1 lc rgb ""red"" pt 7 ps 0.20   # circle    " 
write(UNIT202,'(a)')  " unset key  " 
write(UNIT202,'(a)')  " set tics scale 0.75  " 
write(UNIT202,'(a)')  " set xtics 1 " 
write(UNIT202,'(a)')  " set ytics 1  "
write(UNIT202,'(a)')  " set xlabel 'Kx'  "
write(UNIT202,'(a)')  " set ylabel 'Ky'  "
write(UNIT202,'(a)')  " m = ""./K_ferimilevel.dat""  "
write(UNIT202,'(a)')  "  plot m using 1:2 notitle w p ls 1 "
write(UNIT202,'(a)')  " quit " 
close(UNIT202)

call system('gnuplot -persist  K_ferimilevelplot.plt')
!~ call system('del *.plt')

end subroutine






!--------------------------------------------------------
!                             showgap1plot
!--------------------------------------------------------
subroutine showgap1plot (Tswitch, Y0, Y1, np, ua, uz  )
use moduleone
use moduletwo
implicit real (Kind=DP) (a-h, o-z)
integer(KIND=4)::Tswitch
real (KInd=DP), intent(in) :: Y0, Y1, np, ua, uz
character (len=8) :: material_name

If( Tswitch==1 ) then
material_name="MoS_{2}"
else If( Tswitch==2 ) then
material_name="WS_{2}"
else If( Tswitch==3 ) then
material_name="MoSe_{2}"
else If( Tswitch==4 ) then
material_name="WSe_{2}"
end if
material_name=adjustl(material_name)
material_name=trim(material_name)



open(unit =UNIT13, file = FILENAME13, status="unknown", action="write" , DELIM='none', iostat=ierror) !!  showgap1.plt

!~ write(UNIT13,'(a)')  "set term postscript eps enhanced color 24  " 
!~ write(UNIT13,'(a)')  "set output 'gap1.dat.eps'   "

write(UNIT13,'(a)')  "set terminal pdfcairo enhanced  " 
write(UNIT13,'(a)')  "set output 'gap1_dat.pdf'   "
write(UNIT13,'(a)')  "set style data linespoints   " 
write(UNIT13,'(a)')  "set style line 1 lc rgb ""blue"" pt 7 ps 0.5   # circle    " 
write(UNIT13,'(a)')  "set style line 2 lc rgb ""red"" pt 7 ps 0.5   # circle    " 
write(UNIT13,'(a)')  "set style line 3 lc rgb '#0060ad' pt 9 ps 0.5   # triangle    " 

write(UNIT13,  '(  3a, f6.2, a)'   ) "set key title"" ", material_name, "---- n'=", np, "{/Symbol \264} 10^{13} cm^{-2}"" at screen 0.72, 0.45  "

write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 1 from  187,", Y0, ", 0  to  187,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 2 from  296,", Y0, ", 0  to  296,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 3 from  513,", Y0, ", 0  to  513,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
!~ write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 4 from  730,", Y0, ", 0  to  730,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
!~ write(UNIT13,'(a, f8.4, a, f8.4, a)')  "set arrow 5 from  0,", ua, ", 0  to 839,", ua, ", 0 nohead back nofilled linecolor rgbcolor ""violet""  lw 1.5 dashtype 2  " 
write(UNIT13,'(a, f8.4, a, f8.4, a)')  "set arrow 6 from  0,", uz, ", 0  to 513,", uz, ", 0 nohead back nofilled linecolor rgbcolor ""violet""  lw 1.5  " 

write(UNIT13,'(a)')  "set autoscale x   "
write(UNIT13,'(a)')  "set autoscale y   "
write(UNIT13,'(  2(a, f6.2), a )')  "set yrange [", y0, ":", y1, "]" 

write(UNIT13,'(a)')  "set xtics ( ""{/Symbol G}"" 0, ""M"" 187, ""K"" 296, ""{/Symbol G}"" 513)   font "",26""    "
write(UNIT13,'(a)')  "set ylabel ""{/Symbol=32 e }  (eV) ""   font "",28""  "
write(UNIT13,'(a)')  "plot ""gap1.dat""  using 1:2 notitle w p pt 8 ps 0.5 "
write(UNIT13,'(a)')  "quit " 
close(UNIT13)

call system('gnuplot -persist  showgap1.plt')
!~ call system('del *.plt')

end subroutine

!--------------------------------------------------------
!                             plot Bandstructure0.dat
!--------------------------------------------------------
subroutine showBandstructure0plot ( Tswitch, Y0, Y1)
use moduleone
use moduletwo
implicit real (Kind=DP) (a-h, o-z)

integer(KIND=4)::Tswitch
real (KInd=DP), intent(in) :: Y0, Y1
character (len=8) :: material_name
!~ y0=0.2
!~ y1=1.2

If( Tswitch==1 ) then
material_name="MoS_{2}"
else If( Tswitch==2 ) then
material_name="WS_{2}"
else If( Tswitch==3 ) then
material_name="MoSe_{2}"
else If( Tswitch==4 ) then
material_name="WSe_{2}"
end if
material_name=adjustl(material_name)
material_name=trim(material_name)


open(unit =UNIT14, file = FILENAME14, status="unknown", action="write" , DELIM='none', iostat=ierror) !!  'showBandstructure0.plt'

write(UNIT14,'(a)')  "set term postscript eps enhanced color 24  " 
write(UNIT14,'(a)')  "set output  'Bandstructure0.eps'    "
!~ write(UNIT14,'(a)')  "set terminal pdfcairo enhanced  " 
!~ write(UNIT14,'(a)')  "set output 'Bandstructure0.pdf'   "
write(UNIT14,'(a)')  "set style data linespoints   " 
write(UNIT14,'(a)')  "set style line 1 lc rgb ""blue"" pt 7 ps 0.5   # circle    " 
write(UNIT14,'(a)')  "set style line 2 lc rgb ""red"" pt 7 ps 0.5   # circle    " 
write(UNIT14,'(a)')  "set style line 3 lc rgb '#0060ad' pt 9 ps 0.5   # triangle    " 

write(UNIT14,  '(  3a  )'   ) "set key title"" ", material_name, " ""  "

write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 1 from  187,", Y0, ", 0  to  187,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 2 from  296,", Y0, ", 0  to  296,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 3 from  513,", Y0, ", 0  to  513,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
!~ write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 4 from  730,", Y0, ", 0  to  730,", Y1, ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
!~ write(UNIT14,'(   a    )'                )  "set arrow 5 from  0, 0, 0  to 839, 0, 0 nohead back nofilled linecolor rgbcolor ""violet""  lw 1.5 dashtype 2  " 
write(UNIT14,'(a)')  "set autoscale x   "
write(UNIT14,'(a)')  "set autoscale y   "
write(UNIT14,'(  2(a, f6.2), a )')  "set yrange [", y0, ":", y1, "]" 

!~ write(UNIT14,'(a)')  "set xtics ( ""{/Symbol G}"" 0, ""M"" 187, ""K"" 296, ""{/Symbol G}"" 513, ""K'"" 730, ""M"" 839)   font "",26""    "
write(UNIT14,'(a)')  "set xtics ( ""{/Symbol G}"" 0, ""M"" 187, ""K"" 296, ""{/Symbol G}"" 513)   font "",26""    "
write(UNIT14,'(a)')  "set ylabel ""{/Symbol=32 e }  (eV) ""   font "",28""  "
write(UNIT14,'(a)')  "plot ""Bandstructure0.dat""  using 1:2 notitle w p pt 8 ps 0.5"    !! w p pt 8 ps 0.5
write(UNIT14,'(a)')  "quit " 
close(UNIT14)

call system('gnuplot -persist  showBandstructure0.plt')
!~ call system('del *.plt')

end subroutine



!--------------------------------------------------------
!                             plot dos0.dat
!--------------------------------------------------------
subroutine showdos0plot ( Tswitch, X0, X1)
use moduleone
use moduletwo
implicit real (Kind=DP) (a-h, o-z)

integer(KIND=4)::Tswitch
real (KInd=DP), intent(in) :: X0, X1
character (len=8) :: material_name
!~ y0=0.2
!~ y1=1.2

If( Tswitch==1 ) then
material_name="MoS_{2}"
else If( Tswitch==2 ) then
material_name="WS_{2}"
else If( Tswitch==3 ) then
material_name="MoSe_{2}"
else If( Tswitch==4 ) then
material_name="WSe_{2}"
end if
material_name=adjustl(material_name)
material_name=trim(material_name)

open(unit =UNIT16, file = FILENAME16, status="unknown", action="write" , DELIM='none', iostat=ierror) !!  'showdos0.plt'

write(UNIT16,'(a)')  "set term postscript eps enhanced color 24  " 
write(UNIT16,'(a)')  "set output  'Dos0.eps'    "
!~ write(UNIT16,'(a)')  "set terminal pdfcairo enhanced  " 
!~ write(UNIT16,'(a)')  "set output 'Dos0.pdf'   "
write(UNIT16,'(a)')  "set time"
write(UNIT16,'(a)')  "set notimestamp "
write(UNIT16,'(a)')  "set title ""DOS as a function of energy"" "
write(UNIT16,'(a)')  "set style data line  " 
write(UNIT16,'(a)')  "set autoscale x"
write(UNIT16,'(a)')  "set autoscale y"
write(UNIT16,  '(  3a  )'   ) "set key title"" ", material_name, " ""  "
write(UNIT16,'(a)')  "set xzeroaxis  lw 1.5  "
write(UNIT16,'(  2(a, f6.2), a )')  "set xrange [", x0, ":", x1, "]" 
write(UNIT16,'(a)')  "set xlabel ""E"" "
write(UNIT16,'(a)')  "set ylabel ""Dos(E)"" "
write(UNIT16,'(a)')  "plot ""dos0.dat""  u 1:2 notit  w l "
write(UNIT16,'(a)')  "quit " 
close(UNIT16)
call system('gnuplot -persist  showdos0.plt')

end subroutine



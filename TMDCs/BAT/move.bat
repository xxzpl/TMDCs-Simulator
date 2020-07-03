 @echo off 
 
 set g=%1 
 @move *.plt %g%   
 @copy *.out %g%  
 @move *.pdf %g%
 @move *.dat %g%
 @move *.py %g%
 @move *.eps %g%
 @copy *.txt %g%
 @copy *.bat %g% 

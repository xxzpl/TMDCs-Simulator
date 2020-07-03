    @echo off

    set g=%1%

	if exist %g%\NUL (
	RMDIR /S /Q %g%
	)

    @mkdir %g%

	
	 

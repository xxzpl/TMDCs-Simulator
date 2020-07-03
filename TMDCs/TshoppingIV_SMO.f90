!~ ============================================
!~ program InterAtomshopping
!~ ============================================	
subroutine Hopping_amplitudes (Tswitch, layers, tdp, tdd, tpp, Delta0, Delta1, Delta2, DeltaP, DeltaZ )
    use moduleone
    use moduletwo
    implicit real (kind=DP) (a-h, o-z)

    integer (kind=4) :: Tswitch, layers

    real (kind=DP) :: Delta0, Delta1, Delta2, DeltaP, DeltaZ
    real (kind=DP) :: V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi
    real (kind=DP) :: V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp
    real (kind=DP) :: Vx, Vy, Vz, Vn, a, Vx2, Vy2, Vz2

    real (kind=DP) :: tdp(1:3,1:5,1:6),  tdd(1:5,1:5,1:6),  tpp(1:3,1:3,1:10)

    real(kind=DP) :: pp(1:3, 1:3)

    tdp(1:3, 1:5, 1:6) = 0.d0
    tdd(1:5, 1:5, 1:6) = 0.d0 
    tpp(1:3, 1:3, 1:10) = 0.d0

    open (unit =UNIT1, file = FILENAME1,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  !!===   HoppingAmplitudes.py		
    write(  UNIT1,  '( a )'  ) "import matplotlib  "
    write(  UNIT1,  '( a )'  ) "matplotlib.use('Agg')  "
    write(  UNIT1,  '( a )'  ) "import matplotlib.pyplot as plt  "
    write(  UNIT1,  '( a )'  ) "import numpy as np  "
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) "import sys  "
    write(  UNIT1,  '( a )'  ) "sys.path.append(""../.."") "
    write(  UNIT1,  '( a )'  ) "import tipsi  "

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) "def lattice(a = 0.316, d = 0.1586, z = 0.2968):  "
    write(  UNIT1,  '( 4x, a )'  ) """"""" "
    write(  UNIT1,  '( 4x, a )'  ) 'MoS2 lattice.'
    write(  UNIT1,  '( 4x, a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'Parameters '
    write(  UNIT1,  '( 4x, a )'  ) '------------- '
    write(  UNIT1,  '( 4x, a )'  ) 'a : float '
    write(  UNIT1,  '( 4x, a )'  ) '    lattice constant '
    write(  UNIT1,  '( 4x, a )'  ) 'd : float '
    write(  UNIT1,  '( 4x, a )'  ) '    z-distance between Mo and S atoms '
    write(  UNIT1,  '( 4x, a )'  ) 'z : float '
    write(  UNIT1,  '( 4x, a )'  ) '    interlayer distance '
    write(  UNIT1,  '( 4x, a )'  ) '  '
    write(  UNIT1,  '( 4x, a )'  ) 'Returns'
    write(  UNIT1,  '( 4x, a )'  ) '------------- '
    write(  UNIT1,  '( 4x, a )'  ) 'tipsi.Lattice object '
    write(  UNIT1,  '( 4x, a )'  ) '     XS2 lattice. '
    write(  UNIT1,  '( 4x, a )'  ) """"""""

    write(  UNIT1,  '( 4x, a )'  ) 'Z=z+2.0*d'
    write(  UNIT1,  '( 4x, a )'  ) 'b = a / np.sqrt(3.)'
    write(  UNIT1,  '( 4x, a )'  ) '# lattice vectors'
    write(  UNIT1,  '( 4x, a )'  ) 'vectors        = [[0.5 * a, 1.5 * b, 0],  '
    write(  UNIT1,  '( 4x, a )'  ) '                  [0.5 * a, -1.5 * b, 0], '
    write(  UNIT1,  '( 4x, a )'  ) '                  [0,  0,  - 2.0 * Z]] '

    write(  UNIT1,  '( 4x,a )'  ) '# first the three top S orbitals, then three bottom S orbitals, '
    write(  UNIT1,  '( 4x,a )'  ) '# finally five Mo orbitals'
    write(  UNIT1,  '( 4x,a )'  ) 'orbital_coords = [[0, 0, d],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, d],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, d],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, 0],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, 0],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, 0],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, 0],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, 0],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, -d],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, -d],  '
if (layers==2) then	
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, -d],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, d-Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, d-Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, d-Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, -Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, -Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, -Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, -Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, b, -Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, -d-Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, -d-Z],  '
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, -d-Z]]  '
else if (layers==1) then
    write(  UNIT1,  '( 4x,a )'  ) '                  [0, 0, -d]] '
end if
    write(  UNIT1,  '( 4x,a )'  ) 'return tipsi.Lattice(vectors, orbital_coords)  '
!~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    write(  UNIT1,  '( a )'  ) 	
    write(  UNIT1,  '( a )'  ) "def sheet_rectangle(W, H):  "
    write(  UNIT1,  '( 4x,a )'  ) """"""" "	
    write(  UNIT1,  '( 4x,a )'  ) 'XS2 SiteSet, rectangular.'
    write(  UNIT1,  '( 4x,a )'  ) 
    write(  UNIT1,  '( 4x,a )'  ) 'Parameters '
    write(  UNIT1,  '( 4x,a )'  ) '------------- '
    write(  UNIT1,  '( 4x,a )'  ) 'W : integer '
    write(  UNIT1,  '( 4x,a )'  ) '    width of the SiteSet, in unit cells '
    write(  UNIT1,  '( 4x,a )'  ) 'H : integer '
    write(  UNIT1,  '( 4x,a )'  ) '    height of the SiteSet, in unit cells '
    write(  UNIT1,  '( 4x,a )'  ) '  '
    write(  UNIT1,  '( 4x,a )'  ) 'Returns'
    write(  UNIT1,  '( 4x,a )'  ) '------------- '
    write(  UNIT1,  '( 4x,a )'  ) 'site_set : tipsi.SiteSet object '
    write(  UNIT1,  '( 4x,a )'  ) '     Rectangular XS2 SiteSet.'
    write(  UNIT1,  '( 4x,a )'  ) """"""""

    write(  UNIT1,  '( 4x,a )'  ) 'site_set = tipsi.SiteSet()'
    write(  UNIT1,  '( 4x,a )'  ) 'for x in range(W):'
    write(  UNIT1,  '( 4x,a )'  ) '    for y in range(int(H / 2)):'
    write(  UNIT1,  '( 4x,a )'  ) '        i, j = x + y, x - y  '
    write(  UNIT1,  '( 4x,a )'  ) '        for k in range(11): '
    write(  UNIT1,  '( 4x,a )'  ) '            unit_cell_coords = (i, j, 0) '
    write(  UNIT1,  '( 4x,a )'  ) '            site_set.add_site(unit_cell_coords, k)'
    write(  UNIT1,  '( 4x,a )'  ) 'return site_set'
   
!~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) "def pbc_rectangle(W, H, unit_cell_coords, orbital):  "
    write(  UNIT1,  '( 4x,a )'  ) """"""" "
    write(  UNIT1,  '( 4x,a )'  ) 'PBC for a rectangular XS2 sample.'
    write(  UNIT1,  '( 4x,a )'  ) 
    write(  UNIT1,  '( 4x,a )'  ) 'Parameters '
    write(  UNIT1,  '( 4x,a )'  ) '------------- '
    write(  UNIT1,  '( 4x,a )'  ) 'W : integer '
    write(  UNIT1,  '( 4x,a )'  ) '    width of the sample, in unit cells '
    write(  UNIT1,  '( 4x,a )'  ) 'H : integer '
    write(  UNIT1,  '( 4x,a )'  ) '    height of the SiteSet, in unit cells '
    write(  UNIT1,  '( 4x,a )'  ) 'unit_cell_coords : 3-tuple of integers '
    write(  UNIT1,  '( 4x,a )'  ) '    unit cell coordinates  '
    write(  UNIT1,  '( 4x,a )'  ) 'orbital : integer  '
    write(  UNIT1,  '( 4x,a )'  ) '    orbital index  '
    write(  UNIT1,  '( 4x,a )'  ) 
    write(  UNIT1,  '( 4x,a )'  ) 'Returns'
    write(  UNIT1,  '( 4x,a )'  ) '------------- '
    write(  UNIT1,  '( 4x,a )'  ) 'unit_cell_coords : 3-tuple of integers '
    write(  UNIT1,  '( 4x,a )'  ) '    unit cell coordinates  '
    write(  UNIT1,  '( 4x,a )'  ) 'orbital : integer  '
    write(  UNIT1,  '( 4x,a )'  ) '    orbital index  '
    write(  UNIT1,  '( 4x,a )'  ) """"""""

    write(  UNIT1,  '( 4x,a )'  ) '# get input'
    write(  UNIT1,  '( 4x,a )'  ) 'x, y, z = unit_cell_coords'
    write(  UNIT1,  '( 4x,a )'  ) '# transform to rectangular coords (xloc, yloc)'
    write(  UNIT1,  '( 4x,a )'  ) 'xloc = (x + y) / 2.  '
    write(  UNIT1,  '( 4x,a )'  ) 'yloc = (x - y) / 2. '
    write(  UNIT1,  '( 4x,a )'  ) '# use standard pbc'
    write(  UNIT1,  '( 4x,a )'  ) 'xloc = xloc % W'
    write(  UNIT1,  '( 4x,a )'  ) 'yloc = yloc % (H / 2)'

    write(  UNIT1,  '( 4x,a )'  ) '# transform back  '
    write(  UNIT1,  '( 4x,a )'  ) 'x = int(xloc + yloc) '
    write(  UNIT1,  '( 4x,a )'  ) 'y = int(xloc - yloc)'
    write(  UNIT1,  '( 4x,a )'  ) '# done'
    write(  UNIT1,  '( 4x,a )'  ) 'return (x, y, z), orbital'

!~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write(  UNIT1,  '( a )'  ) 
write(  UNIT1,  '( a )'  ) "def hop_dict(X=""S"") :  "
write(  UNIT1,  '( 4x,a )'  ) """"""" "	
write(  UNIT1,  '( 4x,a )'  ) 'MoS2/WS2 hopping dictionary'
write(  UNIT1,  '( 4x,a )'  ) 'Parameters '
write(  UNIT1,  '( 4x,a )'  ) '------------- '
write(  UNIT1,  '( 4x,a )'  ) 'X : string '
write(  UNIT1,  '( 4x,a )'  ) "    Either ""S"" or ""Se"" "
write(  UNIT1,  '( 4x,a )'  ) '  '
write(  UNIT1,  '( 4x,a )'  ) 'Returns'             
write(  UNIT1,  '( 4x,a )'  ) '------------- '
write(  UNIT1,  '( 4x,a )'  ) 'hops : tipsi.HopDict object '
write(  UNIT1,  '( 4x,a )'  ) '     XS2 HopDict. '
write(  UNIT1,  '( 4x,a )'  ) """"""""
    !===========================================

    !	MoS2
    if (Tswitch==1) then
        a=AMos2*10.0**(10)                        !  structure constant    a
        u=distanceoflayermos2*10.0**(10)
        cprime=CprimeMoS2*10.0**(10)
        w=distanceoflayermos2W*10.0**(10)

    !	WS2
    else if (Tswitch==2) then
        a=AWs2*10.0**(10) 
        u=distanceoflayerWs2*10.0**(10)
        cprime=CprimeWS2*10.0**(10)
        w=distanceoflayerWs2W*10.0**(10)
    
    !	MoSe2
    else if (Tswitch==3) then
        a=3.288d0
        u=1.664d0
        cprime=6.4510d0
        w=3.123d0

    !	WSe2
    else if (Tswitch==4) then
        a=3.260d0
        u=1.657d0
        cprime=6.422d0
        w=3.108d0
    end if


    !~ =============================================	
    !   Next Nearest Hopping
    !
    !           Mo-Mo  / W-W
    !
    !       6(alpha)    5(-gamma)
    !   1(-beta)            4(beta)
    !       2(gamma)    3(-alpha)
    !

    !~ ===============  d - d family  ==============================
    
    !~                       V-alpha
    !~ =============================================================

    Vx=  -0.5d0  ;       Vy=   0.5d0*sq3   ;      Vz= 0.d0   
    
    call    DD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
            & E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2, &
            & E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2, &
            & E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2, &
            & E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2, &
            & E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2 )
    !~ =============================================

    tdd( 1, 1:5 , 6 )=[  E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2  ]
    tdd( 2, 1:5 , 6 )=[  E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2  ]
    tdd( 3, 1:5 , 6 )=[  E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2  ]
    tdd( 4, 1:5 , 6 )=[  E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2 ]
    tdd( 5, 1:5 , 6 )=[  E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2  ]


    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '###################################### '
    write(  UNIT1,  '( 4x, a )'  ) '# define site-to-site hopping matrices	  '
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '# Mo-Mo, next-nearest neighbours ' 
    write(  UNIT1,  '( 4x, a )'  ) '# in three directions (alpha, beta, gamma)   '
    write(  UNIT1,  '( 4x, a )'  ) '#  '
    write(  UNIT1,  '( 4x, a )'  ) '# Mo_alpha   ' 
    write(  UNIT1,  '( 4x, a )'  ) '#       Mo  Mo_beta '
    write(  UNIT1,  '( 4x, a )'  ) '# Mo_gamma  '
    write(  UNIT1,  '( 4x, a )'  ) '#   '
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'A_Mo_Mo_nnb_alpha = np.zeros((5, 5), dtype = float)    ' 
    write(  UNIT1,  '( 4x, a )'  ) 'A_Mo_Mo_nnb_beta = np.zeros((5, 5), dtype = float)      '
    write(  UNIT1,  '( 4x, a )'  ) 'A_Mo_Mo_nnb_gamma = np.zeros((5, 5), dtype = float)  '
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 

    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[0,0] = ',  tdd(1, 1, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[0,1] = ',  tdd(1, 2, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[0,2] = ',  tdd(1, 3, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[0,3] = ',  tdd(1, 4, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[0,4] = ',  tdd(1, 5, 6)
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[1,0] = ',  tdd(2, 1, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[1,1] = ',  tdd(2, 2, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[1,2] = ',  tdd(2, 3, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[1,3] = ',  tdd(2, 4, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[1,4] = ',  tdd(2, 5, 6)
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[2,0] = ',  tdd(3, 1, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[2,1] = ',  tdd(3, 2, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[2,2] = ',  tdd(3, 3, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[2,3] = ',  tdd(3, 4, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[2,4] = ',  tdd(3, 5, 6)
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[3,0] = ',  tdd(4, 1, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[3,1] = ',  tdd(4, 2, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[3,2] = ',  tdd(4, 3, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[3,3] = ',  tdd(4, 4, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[3,4] = ',  tdd(4, 5, 6)
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[4,0] = ',  tdd(5, 1, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[4,1] = ',  tdd(5, 2, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[4,2] = ',  tdd(5, 3, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[4,3] = ',  tdd(5, 4, 6)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_alpha[4,4] = ',  tdd(5, 5, 6)
    
    !~ =======================================================================
    !~                    DD            V-Beta
    !~ =======================================================================

    Vx=  1.0d0*a  ;       Vy=   0.d0   ;      Vz= 0.d0   

	call	 DD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2, &
			& E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2, &
			& E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2, &
			& E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2, &
			& E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2 )
	!~ =============================================

	tdd( 1, 1:5 , 4 )=[  E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2    ]
	tdd( 2, 1:5 , 4 )=[  E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2    ]
	tdd( 3, 1:5 , 4 )=[  E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2    ]
	tdd( 4, 1:5 , 4 )=[  E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2 ]
	tdd( 5, 1:5 , 4 )=[  E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2    ]

	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[0,0] = ',  tdd(1, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[0,1] = ',  tdd(1, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[0,2] = ',  tdd(1, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[0,3] = ',  tdd(1, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[0,4] = ',  tdd(1, 5, 4)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[1,0] = ',  tdd(2, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[1,1] = ',  tdd(2, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[1,2] = ',  tdd(2, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[1,3] = ',  tdd(2, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[1,4] = ',  tdd(2, 5, 4)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[2,0] = ',  tdd(3, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[2,1] = ',  tdd(3, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[2,2] = ',  tdd(3, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[2,3] = ',  tdd(3, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[2,4] = ',  tdd(3, 5, 4)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[3,0] = ',  tdd(4, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[3,1] = ',  tdd(4, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[3,2] = ',  tdd(4, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[3,3] = ',  tdd(4, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[3,4] = ',  tdd(4, 5, 4)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[4,0] = ',  tdd(5, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[4,1] = ',  tdd(5, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[4,2] = ',  tdd(5, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[4,3] = ',  tdd(5, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_beta[4,4] = ',  tdd(5, 5, 4)
	
	!~ =========================================================================
	!~                   DD            V-Gamma   2
	!~ =========================================================================

	Vx=  -0.5d0*a  ;       Vy=   -0.5d0*a*sq3   ;      Vz= 0.d0   
	

	call	 DD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2, &
			& E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2, &
			& E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2, &
			& E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2, &
			& E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2 )

	!~ =============================================

	tdd( 1, 1:5 , 2 )=[  E_xy_xy, E_xy_yz, E_xy_xz, E_xy_x2y2, E_xy_z2    ]
	tdd( 2, 1:5 , 2 )=[  E_yz_xy,  E_yz_yz , E_yz_xz , E_yz_x2y2 , E_yz_z2    ]
	tdd( 3, 1:5 , 2 )=[  E_xz_xy,  E_xz_yz , E_xz_xz , E_xz_x2y2 , E_xz_z2    ]
	tdd( 4, 1:5 , 2 )=[  E_x2y2_xy,  E_x2y2_yz , E_x2y2_xz , E_x2y2_x2y2 , E_x2y2_z2 ]
	tdd( 5, 1:5 , 2 )=[  E_z2_xy,  E_z2_yz , E_z2_xz , E_z2_x2y2 , E_z2_z2    ]

	tdd(:,:,1)=tdd(:,:,4)
	tdd(:,:,3)=tdd(:,:,6)
	tdd(:,:,5)=tdd(:,:,2)
	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[0,0] = ',  tdd(1, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[0,1] = ',  tdd(1, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[0,2] = ',  tdd(1, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[0,3] = ',  tdd(1, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[0,4] = ',  tdd(1, 5, 2)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[1,0] = ',  tdd(2, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[1,1] = ',  tdd(2, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[1,2] = ',  tdd(2, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[1,3] = ',  tdd(2, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[1,4] = ',  tdd(2, 5, 2)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[2,0] = ',  tdd(3, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[2,1] = ',  tdd(3, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[2,2] = ',  tdd(3, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[2,3] = ',  tdd(3, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[2,4] = ',  tdd(3, 5, 2)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[3,0] = ',  tdd(4, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[3,1] = ',  tdd(4, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[3,2] = ',  tdd(4, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[3,3] = ',  tdd(4, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[3,4] = ',  tdd(4, 5, 2)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[4,0] = ',  tdd(5, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[4,1] = ',  tdd(5, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[4,2] = ',  tdd(5, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[4,3] = ',  tdd(5, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_nnb_gamma[4,4] = ',  tdd(5, 5, 2)

	!===========================================
	!	Nearest Hoping
	!~ W-S
	!~ W-Se
	!~ M0-S

				!~     Mo
				!~ 1(Beta)
				!~ S  Top
		!~ 3 (Alpha)      2(Gamma)

	!~ ==============================================================================
	!~                        PD            V-Beta S--top  1
	!~ ==============================================================================

	Vx=  0.d0  ;       Vy=   1.d0/sq3   ;      Vz= -0.5d0 

	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)
	!~ =============================================

	tdp(1, 1:5, 1)= [ E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2  ]
	tdp(2, 1:5, 1)= [ E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2  ]
    tdp(3, 1:5, 1)=[ E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2 ]
	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) '# S_top-Mo, nearest neighbours ' 
	write(  UNIT1,  '( 4x, a )'  ) '# in three directions (alpha, beta, gamma)   '
	write(  UNIT1,  '( 4x, a )'  ) '#  '
	write(  UNIT1,  '( 4x, a )'  ) '#      Mo_beta    ' 
	write(  UNIT1,  '( 4x, a )'  ) '#         S   '
	write(  UNIT1,  '( 4x, a )'  ) '# Mo_alpha  Mo_gamma  '
	write(  UNIT1,  '( 4x, a )'  ) '# HE TMDC XIANG FAN BOTTOM  '
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) 'A_Stop_Mo_nb_alpha = np.zeros((3, 5), dtype = float)  ' 
	write(  UNIT1,  '( 4x, a )'  ) 'A_Stop_Mo_nb_beta = np.zeros((3, 5), dtype = float)      '
	write(  UNIT1,  '( 4x, a )'  ) 'A_Stop_Mo_nb_gamma = np.zeros((3, 5), dtype = float)  '
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[0,0] = ',  tdp(1, 1, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[0,1] = ',  tdp(1, 2, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[0,2] = ',  tdp(1, 3, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[0,3] = ',  tdp(1, 4, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[0,4] = ',  tdp(1, 5, 1)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[1,0] = ',  tdp(2, 1, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[1,1] = ',  tdp(2, 2, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[1,2] = ',  tdp(2, 3, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[1,3] = ',  tdp(2, 4, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[1,4] = ',  tdp(2, 5, 1)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[2,0] = ',  tdp(3, 1, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[2,1] = ',  tdp(3, 2, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[2,2] = ',  tdp(3, 3, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[2,3] = ',  tdp(3, 4, 1)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_beta[2,4] = ',  tdp(3, 5, 1)
	
	!~ =============================================	
	!~                                  PD            V-Gamma S--top  2
	!~ =============================================	

	Vx=  0.5d0*a  ;       Vy=   -0.5d0/sq3*a   ;      Vz= -0.5d0*a   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)

	tdp(1, 1:5, 2)=[ E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2  ]
	tdp(2, 1:5, 2)=[ E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2  ]
	tdp(3, 1:5, 2)=[ E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2  ]

	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[0,0] = ',  tdp(1, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[0,1] = ',  tdp(1, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[0,2] = ',  tdp(1, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[0,3] = ',  tdp(1, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[0,4] = ',  tdp(1, 5, 2)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[1,0] = ',  tdp(2, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[1,1] = ',  tdp(2, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[1,2] = ',  tdp(2, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[1,3] = ',  tdp(2, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[1,4] = ',  tdp(2, 5, 2)	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[2,0] = ',  tdp(3, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[2,1] = ',  tdp(3, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[2,2] = ',  tdp(3, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[2,3] = ',  tdp(3, 4, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_gamma[2,4] = ',  tdp(3, 5, 2)	

	!~ =============================================	
	!~                                  PD            V-alpha S--top   3
	!~ =============================================	

	Vx=  -0.5d0*a  ;       Vy=   -0.5d0/sq3*a   ;      Vz= -0.5d0*a   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)


	tdp(1, 1:5, 3)=[ E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2  ]
	tdp(2, 1:5, 3)=[ E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2  ]
	tdp(3, 1:5, 3)=[ E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2  ]	

	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[0,0] = ',  tdp(1, 1, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[0,1] = ',  tdp(1, 2, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[0,2] = ',  tdp(1, 3, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[0,3] = ',  tdp(1, 4, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[0,4] = ',  tdp(1, 5, 3)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[1,0] = ',  tdp(2, 1, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[1,1] = ',  tdp(2, 2, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[1,2] = ',  tdp(2, 3, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[1,3] = ',  tdp(2, 4, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[1,4] = ',  tdp(2, 5, 3)	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[2,0] = ',  tdp(3, 1, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[2,1] = ',  tdp(3, 2, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[2,2] = ',  tdp(3, 3, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[2,3] = ',  tdp(3, 4, 3)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Mo_nb_alpha[2,4] = ',  tdp(3, 5, 3)	

	!!~ Bottom
	!	Nearest Hoping
	!~ W-S
	!~ W-Se
	!~ M0-S
	!	Mo-Se	Bottom layer

				!~ Mo
			!~ 4(Beta)
				!~ S  Bottom
	!~ 6 (Alpha)      5(Gamma)


	!~ =============================================	
	!~                                  PD            V-Beta S bottom  4
	!~ =============================================	

	Vx=  0.d0  ;       Vy=   1.d0/sq3*a   ;      Vz= 0.5d0*a   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)

	tdp(1, 1:5, 4)=[ E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2  ]
	tdp(2, 1:5, 4)=[ E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2  ]
	tdp(3, 1:5, 4)=[ E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2  ]	

	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) '# S_bottom-Mo, nearest neighbours ' 	
	write(  UNIT1,  '( 4x, a )'  ) '# in three directions (alpha, beta, gamma)   '
	write(  UNIT1,  '( 4x, a )'  ) '#	  '	
	write(  UNIT1,  '( 4x, a )'  ) '#      Mo_beta    ' 	
	write(  UNIT1,  '( 4x, a )'  ) '#         S   '
	write(  UNIT1,  '( 4x, a )'  ) '# Mo_alpha  Mo_gamma  '
	write(  UNIT1,  '( 4x, a )'  ) '# HE TMDC XIANG FAN BOTTOM	  '	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) 'A_Sbot_Mo_nb_alpha = np.zeros((3, 5), dtype = float)  ' 	
	write(  UNIT1,  '( 4x, a )'  ) 'A_Sbot_Mo_nb_beta = np.zeros((3, 5), dtype = float)      '
	write(  UNIT1,  '( 4x, a )'  ) 'A_Sbot_Mo_nb_gamma = np.zeros((3, 5), dtype = float)	  '	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[0,0] = ',  tdp(1, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[0,1] = ',  tdp(1, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[0,2] = ',  tdp(1, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[0,3] = ',  tdp(1, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[0,4] = ',  tdp(1, 5, 4)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[1,0] = ',  tdp(2, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[1,1] = ',  tdp(2, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[1,2] = ',  tdp(2, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[1,3] = ',  tdp(2, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[1,4] = ',  tdp(2, 5, 4)	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[2,0] = ',  tdp(3, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[2,1] = ',  tdp(3, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[2,2] = ',  tdp(3, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[2,3] = ',  tdp(3, 4, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_beta[2,4] = ',  tdp(3, 5, 4)	
	!~ =============================================	
	!~                                  PD            V-Gamma S bottom  5
	!~ =============================================	

	Vx=  0.5d0*a  ;       Vy=   -0.5d0/sq3*a   ;      Vz= 0.5d0*a   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)

	
	tdp(1, 1:5, 5)=[ E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2  ]
	tdp(2, 1:5, 5)=[ E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2  ]
	tdp(3, 1:5, 5)=[ E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2  ]		

	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[0,0] = ',  tdp(1, 1, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[0,1] = ',  tdp(1, 2, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[0,2] = ',  tdp(1, 3, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[0,3] = ',  tdp(1, 4, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[0,4] = ',  tdp(1, 5, 5)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[1,0] = ',  tdp(2, 1, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[1,1] = ',  tdp(2, 2, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[1,2] = ',  tdp(2, 3, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[1,3] = ',  tdp(2, 4, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[1,4] = ',  tdp(2, 5, 5)	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[2,0] = ',  tdp(3, 1, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[2,1] = ',  tdp(3, 2, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[2,2] = ',  tdp(3, 3, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[2,3] = ',  tdp(3, 4, 5)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_gamma[2,4] = ',  tdp(3, 5, 5)	

	!~ =============================================	
	!~                                  PD            V-alpha S bottom  6
	!~ =============================================	

	Vx=  -0.5d0*a  ;       Vy=   -0.5d0/sq3*a   ;      Vz= 0.5d0*a   


	call	 PD_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
			& E_x_xy,  E_x_yz,  E_x_xz,  E_x_x2y2, E_x_z2, &
			& E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2, &
			& E_z_xy,  E_z_yz ,  E_z_xz,  E_z_x2y2 , E_z_z2)


	tdp(1, 1:5, 6)=[ E_x_xy, E_x_yz, E_x_xz, E_x_x2y2, E_x_z2  ]
	tdp(2, 1:5, 6)=[ E_y_xy,  E_y_yz , E_y_xz , E_y_x2y2 , E_y_z2  ]
	tdp(3, 1:5, 6)=[ E_z_xy,  E_z_yz , E_z_xz , E_z_x2y2 , E_z_z2  ]	
	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[0,0] = ',  tdp(1, 1, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[0,1] = ',  tdp(1, 2, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[0,2] = ',  tdp(1, 3, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[0,3] = ',  tdp(1, 4, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[0,4] = ',  tdp(1, 5, 6)
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[1,0] = ',  tdp(2, 1, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[1,1] = ',  tdp(2, 2, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[1,2] = ',  tdp(2, 3, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[1,3] = ',  tdp(2, 4, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[1,4] = ',  tdp(2, 5, 6)	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[2,0] = ',  tdp(3, 1, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[2,1] = ',  tdp(3, 2, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[2,2] = ',  tdp(3, 3, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[2,3] = ',  tdp(3, 4, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Sbot_Mo_nb_alpha[2,4] = ',  tdp(3, 5, 6)		

	
					!	S-S   / Se-Se
					!	x-1, y-2, z-3
					!		6(alpha)	5(-gamma)
					!	1(-beta)			4(beta)
					!		2(gamma)	3(-alpha)


	!~ =============================================	
	!~                                  PP           V-alpha  6
	!~ =============================================	

	Vx=   -0.5d0*a  ;       Vy=   0.5d0*sq3*a  ;      Vz= 0.0d0   


	call	 PP_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
				& E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
				& E_z_x,  E_z_y , E_z_z)
	
	tpp(1, 1:3, 6)= [ E_x_x, E_x_y,  E_x_z ]	
	tpp(2, 1:3, 6)= [ E_y_x, E_y_y,  E_y_z ]	
	tpp(3, 1:3, 6)= [ E_z_x, E_z_y,  E_z_z  ]	


	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) '# S_top-S_top / S_bottom-S_bottom, next-nearest neighbours  ' 	
	write(  UNIT1,  '( 4x, a )'  ) '# in three directions (alpha, beta, gamma)     '
	write(  UNIT1,  '( 4x, a )'  ) '#	  '	
	write(  UNIT1,  '( 4x, a )'  ) '# S_alpha    ' 	
	write(  UNIT1,  '( 4x, a )'  ) '#      S  S_beta    '
	write(  UNIT1,  '( 4x, a )'  ) '# S_gamma	  '
	write(  UNIT1,  '( 4x, a )'  ) '#	  '	
	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) 'A_S_S_nnb_alpha = np.zeros((3, 3), dtype = float)   ' 	
	write(  UNIT1,  '( 4x, a )'  ) 'A_S_S_nnb_beta = np.zeros((3, 3), dtype = float)      '
	write(  UNIT1,  '( 4x, a )'  ) 'A_S_S_nnb_gamma = np.zeros((3, 3), dtype = float)	  '	
	write(  UNIT1,  '( a )'  ) 


	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[0,0] = ',  tpp(1, 1, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[0,1] = ',  tpp(1, 2, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[0,2] = ',  tpp(1, 3, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[1,0] = ',  tpp(2, 1, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[1,1] = ',  tpp(2, 2, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[1,2] = ',  tpp(2, 3, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[2,0] = ',  tpp(3, 1, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[2,1] = ',  tpp(3, 2, 6)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_alpha[2,2] = ',  tpp(3, 3, 6)

	!~ =========================================================================================
	!~                                  PP            V-Beta   4
	!~ =========================================================================================

	Vx=  1.0  ;       Vy=   0.d0   ;      Vz= 0.0d0   


	call	 PP_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
				& E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
				& E_z_x,  E_z_y , E_z_z)

	tpp(1, 1:3, 4)= [ E_x_x, E_x_y,  E_x_z ]
	tpp(2, 1:3, 4)= [ E_y_x, E_y_y,  E_y_z ]
	tpp(3, 1:3, 4)= [ E_z_x, E_z_y,  E_z_z ]
	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[0,0] = ',  tpp(1, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[0,1] = ',  tpp(1, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[0,2] = ',  tpp(1, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[1,0] = ',  tpp(2, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[1,1] = ',  tpp(2, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[1,2] = ',  tpp(2, 3, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[2,0] = ',  tpp(3, 1, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[2,1] = ',  tpp(3, 2, 4)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_beta[2,2] = ',  tpp(3, 3, 4)

	
	!~ =============================================
	!~               PP            V-Gamma   2
	!~ =============================================

	Vx=   -0.5d0*a  ;       Vy=   -0.5d0*sq3*a  ;      Vz= 0.0d0   


	call	 PP_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
				& E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
				& E_z_x,  E_z_y , E_z_z)

	tpp(1, 1:3, 2)= [ E_x_x, E_x_y,  E_x_z ]
	tpp(2, 1:3, 2)= [ E_y_x, E_y_y,  E_y_z ]
	tpp(3, 1:3, 2)= [ E_z_x, E_z_y,  E_z_z ]


	tpp(:,:,1)=tpp(:,:,4)
	tpp(:,:,3)=tpp(:,:,6)
	tpp(:,:,5)=tpp(:,:,2)
	
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[0,0] = ',  tpp(1, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[0,1] = ',  tpp(1, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[0,2] = ',  tpp(1, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[1,0] = ',  tpp(2, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[1,1] = ',  tpp(2, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[1,2] = ',  tpp(2, 3, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[2,0] = ',  tpp(3, 1, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[2,1] = ',  tpp(3, 2, 2)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_nnb_gamma[2,2] = ',  tpp(3, 3, 2)


	
			!	Interlayer Hoping
			!	S-S	
			!    8 (u/ gamma)		9 (v/ alpha)
			!		S
			!		7 (lambda /beta) 
			!	(Mo,S,Direction)  
	
	!~ =============================================
	!~     PP            V-    lambda /beta   7
	!~ =============================================

	Vx=   0.d0  ;       Vy=   a/sq3  ;      Vz= -w


	call	 PPINTERLAYER_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
				& E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
				& E_z_x,  E_z_y , E_z_z)

	tpp(1, 1:3, 7)= [ E_x_x, E_x_y,  E_x_z ]
	tpp(2, 1:3, 7)= [ E_y_x, E_y_y,  E_y_z ]
	tpp(3, 1:3, 7)= [ E_z_x, E_z_y,  E_z_z ]

	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) '# S-S interlayer  ' 
	write(  UNIT1,  '( 4x, a )'  ) '# in three directions (alpha, beta, gamma)    '
	write(  UNIT1,  '( 4x, a )'  ) '#  '
	write(  UNIT1,  '( 4x, a )'  ) '# Stop_gamma   Stop_alpha    ' 
	write(  UNIT1,  '( 4x, a )'  ) '#             S    '
	write(  UNIT1,  '( 4x, a )'  ) '#         Stop_beta  '
	write(  UNIT1,  '( 4x, a )'  ) '#  '
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a )'  ) 'A_S_S_inter_alpha = np.zeros((3, 3), dtype = float)     ' 
	write(  UNIT1,  '( 4x, a )'  ) 'A_S_S_inter_beta = np.zeros((3, 3), dtype = float)       '
	write(  UNIT1,  '( 4x, a )'  ) 'A_S_S_inter_gamma = np.zeros((3, 3), dtype = float)	  '
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( a )'  ) 
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[0,0] = ',  tpp(1, 1, 7)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[0,1] = ',  tpp(1, 2, 7)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[0,2] = ',  tpp(1, 3, 7)
        
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[1,0] = ',  tpp(2, 1, 7)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[1,1] = ',  tpp(2, 2, 7)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[1,2] = ',  tpp(2, 3, 7)

	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[2,0] = ',  tpp(3, 1, 7)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[2,1] = ',  tpp(3, 2, 7)
	write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_beta[2,2] = ',  tpp(3, 3, 7)

	!~ =============================================	
	!~                                  PP            V-   u/ gamma   8
	!~ =============================================	

	Vx=   0.5d0*a  ;       Vy=   -0.5d0*a/sq3  ;      Vz= -w


	call	 PPINTERLAYER_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
				& E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
				& E_z_x,  E_z_y , E_z_z)

    tpp(1, 1:3, 8)= [ E_x_x, E_x_y,  E_x_z ]
    tpp(2, 1:3, 8)= [ E_y_x, E_y_y,  E_y_z ]
    tpp(3, 1:3, 8)= [ E_z_x, E_z_y,  E_z_z ]

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 	
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[0,0] = ',  tpp(1, 1, 8)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[0,1] = ',  tpp(1, 2, 8)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[0,2] = ',  tpp(1, 3, 8)
        
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[1,0] = ',  tpp(2, 1, 8)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[1,1] = ',  tpp(2, 2, 8)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[1,2] = ',  tpp(2, 3, 8)

    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[2,0] = ',  tpp(3, 1, 8)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[2,1] = ',  tpp(3, 2, 8)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_gamma[2,2] = ',  tpp(3, 3, 8)

    !~ =============================================
    !~                                  PP            V-    v/ alpha 9
    !~ =============================================

    Vx=   -0.5d0*a  ;       Vy=   -0.5d0*a/sq3  ;      Vz= -w


    call    PPINTERLAYER_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
                & E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
                & E_z_x,  E_z_y , E_z_z)

    tpp(1, 1:3, 9)= [ E_x_x, E_x_y,  E_x_z ]
    tpp(2, 1:3, 9)= [ E_y_x, E_y_y,  E_y_z ]
    tpp(3, 1:3, 9)= [ E_z_x, E_z_y,  E_z_z ]
    
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 	
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[0,0] = ',  tpp(1, 1, 9)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[0,1] = ',  tpp(1, 2, 9)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[0,2] = ',  tpp(1, 3, 9)
        
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[1,0] = ',  tpp(2, 1, 9)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[1,1] = ',  tpp(2, 2, 9)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[1,2] = ',  tpp(2, 3, 9)

    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[2,0] = ',  tpp(3, 1, 9)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[2,1] = ',  tpp(3, 2, 9)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_inter_alpha[2,2] = ',  tpp(3, 3, 9)



    !~ tpp(1:3, 1:3, 7:9)= tpp(1:3, 1:3, 7:9)*1.0
    !~ pp= tpp(1:3,1:3, 7)
    !~ tpp(1:3, 1:3, 7:8)=tpp(1:3, 1:3, 8:9)
    !~ tpp(1:3, 1:3, 9)=pp

    !~ =============================================	
    !~                                  PP            V-    Vertical   10
    !~ =============================================	

    Vx=   0.d0  ;       Vy=   0.d0  ;      Vz= 1.0d0   


    call     PP_Bonds_Family( Tswitch, Vx, Vy, Vz, &  
                & E_x_x, E_x_y, E_x_z, E_y_x,  E_y_y , E_y_z, &
                & E_z_x,  E_z_y , E_z_z)

    tpp(1, 1:3, 10)= [ E_x_x, E_x_y,  E_x_z ]
    tpp(2, 1:3, 10)= [ E_y_x, E_y_y,  E_y_z ]
    tpp(3, 1:3, 10)= [ E_z_x, E_z_y,  E_z_z  ]

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '# S_top-S-bot, within the unit cell	  '
    write(  UNIT1,  '( 4x, a )'  ) 'A_Stop_Sbot_uc = np.zeros((3, 3), dtype = float)   ' 

    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[0,0] = ',  tpp(1, 1, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[0,1] = ',  tpp(1, 2, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[0,2] = ',  tpp(1, 3, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[1,0] = ',  tpp(2, 1, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[1,1] = ',  tpp(2, 2, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[1,2] = ',  tpp(2, 3, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[2,0] = ',  tpp(3, 1, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[2,1] = ',  tpp(3, 2, 10)
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Stop_Sbot_uc[2,2] = ',  tpp(3, 3, 10)


!!~ Get bond integerals
call Bond_integrals ( Tswitch, Delta0, Delta1, Delta2, DeltaP, DeltaZ, &  
        & V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, &
        & V_pp_sigma_interlayer,   V_pp_pi_interlayer,   SOC_dd,  SOC_pp)
    
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '# S_top-S-top / S_bottom-S_bottom, onsite  '
    write(  UNIT1,  '( 4x, a )'  ) 'A_S_S_uc = np.zeros((3, 3), dtype = float)   ' 
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_uc[0,0] = ', DeltaP
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_uc[1,1] = ', DeltaP
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_S_S_uc[2,2] = ', DeltaZ

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '# Mo-Mo, onsite  '
    write(  UNIT1,  '( 4x, a )'  ) 'A_Mo_Mo_uc = np.zeros((5, 5), dtype = float)  ' 
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_uc[0,0] = ', Delta2
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_uc[1,1] = ', Delta1
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_uc[2,2] = ', Delta1
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_uc[3,3] = ', Delta2
    write(  UNIT1,  '( 4x, a, F20.13 )'  ) 'A_Mo_Mo_uc[4,4] = ', Delta0

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '###################################### '
    write(  UNIT1,  '( 4x, a )'  ) '###################################### ' 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '# set empty hopping matrices'
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict = tipsi.HopDict()   '

if (layers==1) then	

    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((0,0,0), (11,11))  # uc; nb: beta  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((-1,0,0), (11,11)) # nb: alpha; nnb: gamma  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((0,1,0), (11,11)) # nb: gamma  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((0,-1,0), (11,11)) # nnb: alpha  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((1,1,0), (11,11))  # nnb: beta  '

else if( layers==2 ) then

    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((0,0,0), (22,22))  # uc; nb: beta   '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((-1,0,0), (22,22)) # nb: alpha; nnb: gamma  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((0,1,0), (22,22))  # nb: gamma    N2- nnb: alpha '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((0,-1,0), (22,22)) # nnb: alpha   N2- nb: gamma  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((1,1,0), (22,22))  # nnb: beta   '
    write(  UNIT1,  '( a )'  ) 	
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((1,0,0), (22,22))  # N2- nb: alpha; nnb: gamma '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.empty((-1,-1,0), (22,22))  # N2- nnb: beta   '
end if

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) '# fill hopping matrices  '
    write(  UNIT1,  '( 4x, a )'  ) '# conjugates are added automatically  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][0:3,0:3] = A_S_S_uc[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][8:11,8:11] = A_S_S_uc[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][3:8,3:8] = A_Mo_Mo_uc[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][0:3,8:11] = A_Stop_Sbot_uc[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][0:3,3:8] = A_Stop_Mo_nb_beta[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][8:11,3:8] = A_Sbot_Mo_nb_beta[:,:]  '

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,0,0)][0:3,3:8] = A_Stop_Mo_nb_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,0,0)][8:11,3:8] = A_Sbot_Mo_nb_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,0,0)][0:3,0:3] = A_S_S_nnb_gamma[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,0,0)][8:11,8:11] = A_S_S_nnb_gamma[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,0,0)][3:8,3:8] = A_Mo_Mo_nnb_gamma[:,:]  '

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,1,0)][0:3,3:8] = A_Stop_Mo_nb_gamma[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,1,0)][8:11,3:8] = A_Sbot_Mo_nb_gamma[:,:] '

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,-1,0)][0:3,0:3] = A_S_S_nnb_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,-1,0)][8:11,8:11] = A_S_S_nnb_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,-1,0)][3:8,3:8] = A_Mo_Mo_nnb_alpha[:,:]  '

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,1,0)][0:3,0:3] = A_S_S_nnb_beta[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,1,0)][8:11,8:11] = A_S_S_nnb_beta[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,1,0)][3:8,3:8] = A_Mo_Mo_nnb_beta[:,:]  '

if (layers==2) then
    write(  UNIT1,  '( 4x, a )'  ) '# N2 layer'
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][11:14,11:14] = A_S_S_uc[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][19:22,19:22] = A_S_S_uc[:,:]    '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][14:19,14:19] = A_Mo_Mo_uc[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][11:14,19:22] = A_Stop_Sbot_uc[:,:]   '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][11:14,14:19] = A_Stop_Mo_nb_beta[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][19:22,14:19] = A_Sbot_Mo_nb_beta[:,:] '

    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,0,0)][11:14,14:19] = A_Stop_Mo_nb_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,0,0)][19:22,14:19] = A_Sbot_Mo_nb_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,0,0)][11:14,11:14] = A_S_S_nnb_gamma[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,0,0)][19:22,19:22] = A_S_S_nnb_gamma[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(1,0,0)][14:19,14:19] = A_Mo_Mo_nnb_gamma[:,:] '

    write(  UNIT1,  '( a )'  ) 	
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,-1,0)][11:14,14:19] = A_Stop_Mo_nb_gamma[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,-1,0)][19:22,14:19] = A_Sbot_Mo_nb_gamma[:,:] '

    write(  UNIT1,  '( a )'  ) 	
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,1,0)][11:14,11:14] = A_S_S_nnb_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,1,0)][19:22,19:22] = A_S_S_nnb_alpha[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,1,0)][14:19,14:19] = A_Mo_Mo_nnb_alpha[:,:] '

    write(  UNIT1,  '( a )'  ) 	
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,-1,0)][11:14,11:14] = A_S_S_nnb_beta[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,-1,0)][19:22,19:22] = A_S_S_nnb_beta[:,:]   '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,-1,0)][14:19,14:19] = A_Mo_Mo_nnb_beta[:,:] '

    write(  UNIT1,  '( 4x, a )'  ) 	"# interaction"
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,0,0)][8:11,11:14] = A_S_S_inter_beta[:,:] '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(-1,0,0)][8:11,11:14] = A_S_S_inter_alpha[:,:]  '
    write(  UNIT1,  '( 4x, a )'  ) 'hop_dict.dict[(0,1,0)][8:11,11:14] = A_S_S_inter_gamma[:,:] '

end if
    write(  UNIT1,  '( a )'  ) 
    write(  UNIT1,  '( 4x, a )'  ) 'return hop_dict '

    close(UNIT1)
!~ ==
!~ ==
END





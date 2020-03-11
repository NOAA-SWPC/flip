C.................... RSMNSD.FOR ..................................... 
C      PROGRAM FLIP
C.... this main module for the Field Line Interhemispheric Plasma (FLIP)
C.... model. It reads initial conditions and sets up the
C.... coordinate system and neutral atmosphere, note s(4500) for d300
C
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      INTEGER ION,NEQ,WUNIT8,IRIHED,ICAUTN,IVERT,IVPERP,NCOL_EXB,JDAY
     > ,UTORLT,NOCON,AUR_FILE,APFILE
      INTEGER NUMUNIT,UNITNUMS(22)     !.. Unit numbers for writing rate warning
      INTEGER EFLAG(11,11)             !.. solution procedure error flags
      INTEGER INNO  !.. switch to turn on FLIP NO calculation if <0
      REAL SEC,DEC,ETRAN,F107,F107A,AP,BLON,D,T,CHI,SATN,SATF,SAT,
     > GLON,GLONN,GLONF,CHIN,CHIF,GLATN,GLATF,GLATJ,AP3HR(33001)
     > ,GRD,SW(25),UTHRS,UHED(2),DECMAG,SATEQ,L_STAG,L_SHORT
     > ,RE,PLAT,PLON,PI,ANGVEL,TIMEON,TIMEOF,AURALT,ZTOPTE,L_PHASE,SX
     > ,TD_PHASE,DELTIM,PHASE,SECWF8,GLOND,GLATD,BCOMPU,KPREAL,ZFLUX(37)
     > ,UVFAC(59)
      DIMENSION S(50000),D(19),T(2),FD(9),COSDIP(401),Q(401),RESULT(6)
     > ,COOL(22),EDEN(401),EHTJ(401),GL(401),N4S(401),GRAV(401),NNO(401)
      DOUBLE PRECISION XIONN(9,401),XIONV(9,401),TINF(401)
      DOUBLE PRECISION NHEAT(401)    !.. Neutral heating rate

      COMMON/V91/RVBT(401)
      COMMON/PALT/ZLO,ZHI,JWR,JTI
      COMMON/RK/VC(401)
      !.. N= [O+], [H+], sum minor, [He+] densities, Ti= H+,O+,e temps 
      !.. F= cty eqn functions
      COMMON/FJS/N(4,401),TI(3,401),F(20)
      !.. Ap, solar declination, ephemeris transit, magnetic longitude,
      !.. solar indices, UT in secs, magnetic latitude
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      !.. O+,H+ velocities, field factor (obsolete), magnetic field, 
      !.. GR=gravity, R=geocentric distance,SL= distance along B
      COMMON/VN/U(2,401),BG(401),BM(401),GR(2,401),R(2,401),SL(401)
      !.. [O],[H],[N2],[O2], O+ photoionization rate, Tn
      COMMON/ND/ON(401),HN(401),N2N(401),O2N(401),PHION(401),TN(401)
      !.. used in minor ions [O+], V(O+), [He], ???, He+ prod, ???
      COMMON/A/OPLS(401),VOPLS(401),HE(401),UE(2,401),HEPRD(401),OP(401)
      !.. Altitude(km), time step, ???, ???, ???, ITF= change solution
      !.. mode, JMAX=# grid points, JMAX1=JMAX-1,# Newton iterations, # ions
      COMMON/ALT/Z(401),DT,DH,THF,EPS,TF,ITF,JMAX,JMAX1,ITER,ION
      !.. NSAVE= saved [O+] and [H+], O+, H+ flux, EHT=thermal e heating rate
      COMMON/SAV/NSAVE(2,401),TISAV(3,401),FY(2,401),UN(401),EHT(3,401)
      !.. solar zenith angle, ??? ??? ??? ???, max # time steps ??? ???
      COMMON/LPS/SZA(401),EPSN,DC,ISPC,I1,I2,IPV,JTIMAX,IPP,ISKP
      !.. Stores vibrational N2 and odd N densities for printing. 
      !.. RCON= effect of vibrational N2 on O+ + N2 reaction
      COMMON/STAW/STR(12,401),RCON(401)
      !.. RE= Earth radius, (PLAT,PLON) = geographic location of north
      !.. magnetic pole, Pi, angular velocity of Earth 
      COMMON/GPAR/RE,PLAT,PLON,PI,ANGVEL
      !.. Common for O+ - O collision frequencies
      COMMON/CFACT/DHDU_FAC,COLFAC,CFAC(9),ISCH(9)
      !.. Switches for wind type and equations solved
      COMMON/CONFAC/IWIND,IVIBN2,IODDN,INPLUS,IHEP,ITEMP,IHPOP

      !.. JB= grid point of lower boundary, 3 hour Aps, winds flag, SW=MSIS
      !.. switches, ZPROD = upper boundary for EUV ion production
      DATA JB/1/, AP3HR/-1,33000*0/, IWTIM/-1/,SW/15*1,10*1/,ZPROD/700./
      !.. ITSAV ITAU track elapsed time, ISOL= type of solution, HFAC=
      !.. e heating factor, MAXIT= max Newton iters, JTIDT used to set time
      !.. ZLBDY=lower boundary altitude for solving equations, IRIHED IRI/HWM
      !.. wind switch, max solar zenith angle for ion production
      DATA ITSAV,ITAU, ISOL, HFAC, MAXIT, JTIDT, ZLBDY,  IRIHED, SETCHI
     >  / -10000, 0  , -1  ,  0.0 ,  17  ,   5   , 120. ,  0   ,  1.7/
      !.. north/south lat and long 250 km above station, Z0=absolute lower
      !.. boudary altitude, SCAL, GRD scale grid B and vertical grid points
      DATA GLONN, GLONF, GLATN, GLATF, Z0, SCAL, GRD / 7*0.0/,APFILE/-1/
      DATA JHFLX/0/,WUNIT8/1/,NCOL_EXB/0/, JJJ/0/, ICAUTN/0/
      DATA IVPERP/-1/,SCALK/-1.0/,PCO_MAX/0.0/, ISYM/0/, SECWF8/-999/
      DATA HMF2N,HMF2S/250.0,250.0/, JK,JC/10,10/, NOCON/0/
      DATA PCO_LIM/1.03/   !.. minimum allowable L-value
      DATA PCO_WIND/1.08/  !.. minimum allowable L-value for vertical grid
      DATA AUR_FILE/0/     !.. Flag for aurora data file

      !.. Unit numbers for writing reaction rate warning
      DATA NUMUNIT/9/,UNITNUMS/3,4,6,9,12,13,16,56,66,13*0/ 
 888  FORMAT(I4,I6,I6,I5,1P,2E9.2,1P,5E9.1,1P,12E9.1)
      END


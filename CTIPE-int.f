C............................ CTIPE-int.for .........................
C.. This file contains the basic interface routines for the CTIP interface
C::::::::::::::::::::::::::: Module FIELD_LINE_GRID :::::::::::::::::
C.... Contains parameters associated with the field line grid
C.... Written by P. Richards June 2010.
      MODULE FIELD_LINE_GRID
	  IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER JMIN,JMAX   !.. first and last indices on field line grid
        INTEGER, PARAMETER :: FLDIM = 1115     !.. Field line grid dimension
        DOUBLE PRECISION Z(FLDIM),SZA(FLDIM)  !.. altitude, Solar zenith angle
        DOUBLE PRECISION BM(FLDIM)            !.. magnetic field strength
        DOUBLE PRECISION SL(FLDIM)            !.. distance along the field line
        DOUBLE PRECISION GR(FLDIM),GL(FLDIM)  !.. Gravity, magnetic latitude
      END MODULE FIELD_LINE_GRID
C::::::::::::::::::::::::::: ION_DEN_VEL :::::::::::::::::
C.... Ion densities and velocities
C.... Written by P. Richards June 2010.
C.... O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      MODULE ION_DEN_VEL
	  IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: IDIM = 1115  !.. Field line grid dimension
        INTEGER, PARAMETER :: ISPEC = 9   !.. Species dimension
      !.. Ion densities & Velocities
      DOUBLE PRECISION XIONN(ISPEC,IDIM),XIONV(ISPEC,IDIM)
      END MODULE ION_DEN_VEL
C:::::::::::::::::::::::::::Module SOLVAR :::::::::::::::::
C.... This module is used to conserve space used by the Newton solver.
C.... The arrays are workspace for the solver
C.... Written by P. Richards June 2010.
      MODULE SOLVARR
	  IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: SDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION DELTA(10*SDIM),RHS(10*SDIM)
        DOUBLE PRECISION WORK(50*SDIM),S(50*SDIM)
      END MODULE SOLVARR
C::::::::::::::::::::::::::: Module THERMOSPHERE :::::::::::::::::
C.... Contains the neutral densities and temperatures, O+ production, electron 
C.... heating
C.... Written by P. Richards June 2010.
      MODULE THERMOSPHERE
	  IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: TDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION ON(TDIM),HN(TDIM),N2N(TDIM),O2N(TDIM),HE(TDIM)
        DOUBLE PRECISION TN(TDIM),UN(TDIM),TINF(TDIM)
        DOUBLE PRECISION EHT(3,TDIM)
        DOUBLE PRECISION COLFAC  !.. O+ - O collision frequency Burnside factor 1-1.7
      END MODULE THERMOSPHERE
C::::::::::::::::::::::::::: AVE_PARAMS :::::::::::::::::
C.... Computes midpoint values associated with the field line grid
C.... Written by P. Richards June 2010.
      MODULE AVE_PARAMS
	  IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: AVDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION TEJ(AVDIM),TIJ(AVDIM),NUX(2,AVDIM),RBM(AVDIM)
        DOUBLE PRECISION UNJ(AVDIM),DS(AVDIM)
        DOUBLE PRECISION GRADTE(AVDIM),GRADTI(AVDIM),GRAV(AVDIM)
        DOUBLE PRECISION OLOSS(AVDIM),HLOSS(AVDIM),HPROD(AVDIM)
      END MODULE AVE_PARAMS
C::::::::::::::::::::::::::: PRODUCTION :::::::::::::::::
C.... Holds the EUV, photoelectron, and auroral production
C.... Written by P. Richards August 2010.
      !.. Total ionization and excitation rates (EUV + PE + Aurora)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      MODULE PRODUCTION
	  IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: PDIM = 1115  !.. Field line grid dimension
      !.. EUV and photoelectron production rates 
      REAL EUVION(3,12,PDIM),PEXCIT(3,12,PDIM),PEPION(3,12,PDIM)
      DOUBLE PRECISION OTHPR1(6,PDIM),OTHPR2(6,PDIM)
      !.. Total ionization and excitation rates (EUV + PE + Aurora)
      REAL SUMION(3,12,PDIM),SUMEXC(3,12,PDIM)
      !.. Auroral ionization and excitation rates
      REAL PAUION(3,12,PDIM),PAUEXC(3,12,PDIM)
      DOUBLE PRECISION NPLSPRD(PDIM)  !.. N+ production
      DOUBLE PRECISION PHION(PDIM)    !.. Total O+ production
      END MODULE PRODUCTION
C::::::::::::::::::::::::::: Module MINORNEUT :::::::::::::::::
C.... Contains the minor neutral densities 
C.... Written by P. Richards September 2010.
      !.. USE MINORNEUT N4S N2D NNO N2P N2A O1D O1S
      MODULE MINORNEUT
	  IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: NDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION N4S(NDIM),N2D(NDIM),NNO(NDIM),N2P(NDIM),
     >    N2A(NDIM),O1D(NDIM),O1S(NDIM),EQN2D(NDIM)
      END MODULE MINORNEUT
C::::::::::::::::::::::::::::::::: PHOTOIONIZATION_DATA :::::::::::::::
C..... This module contains data for the PRIMPR photoionization routine
C..... Written by P. Richards in August 2010
      MODULE PHOTOIONIZATION_DATA  !.. NPOT LMAX PROB ZLAM SIGION SIGABS TPOT 
	IMPLICIT NONE
      INTEGER, PARAMETER :: LMAX =37  !.. # of EUV wavelengths
      INTEGER NPOT(3)                 !.. # of ionization potentials
      REAL PROB(3,6,LMAX)             !.. Ionization probability,
      REAL ZLAM(LMAX)                 !.. EUV wavelengths,
      REAL SIGION(3,LMAX)             !.. Ionization cross sections
      REAL SIGABS(3,LMAX)             !.. Absorption cross sections
      REAL TPOT(3,10)                 !.. Ionization potentials

      DATA NPOT/5,5,6/   !.. # of O, O2, N2 ionization potentials 
      !.. Ionization potentials
      DATA TPOT/13.6,12.1,15.6,16.9,16.5,16.7,18.6,18.2,18.8,
     >   28.5,20.0,25.0,40.0,25.0,29.0,0.0,0.0,37.0, 12*0.0/
      !.... wavelength data. average is taken for bands
      DATA ZLAM/1025.,1031.91,1025.72,975.,977.02,925.,875.,825.,775.,
     > 789.36,770.41,765.15,725.,703.36,675.,625.,629.73,609.76,575.,
     > 584.33,554.31,525.,475.,465.22,425.,375.,368.07,325.,303.78,
     > 303.31,275.,284.15,256.3,225.,175.,125.,75./

      DATA SIGION/0.0,2.59000E-19,0.00000E+00,0.00000E+00,0.0,
     > 0.00000E+00,0.00000E+00,1.05000E-18,0.00000E+00,0.00000E+00,
     > 1.39400E-17,0.00000E+00,0.00000E+00,1.55400E-17,0.00000E+00,
     > 1.31500E-18,9.37400E-18,0.00000E+00,4.55400E-18,5.49400E-18,
     > 0.00000E+00,3.49800E-18,6.41300E-18,0.00000E+00,5.09100E-18,
     > 1.05970E-17,1.42740E-17,3.74900E-18,1.01910E-17,8.86000E-18,
     > 3.89000E-18,8.47000E-18,8.50000E-18,4.00000E-18,1.17200E-17,
     > 6.58000E-17,1.07360E-17,2.38050E-17,1.50600E-17,1.14600E-17,
     > 2.37500E-17,2.54800E-17,1.72450E-17,2.13060E-17,2.92350E-17,
     > 1.33650E-17,2.49370E-17,2.33390E-17,1.34000E-17,3.11000E-17,
     > 2.33700E-17,1.34000E-17,2.63900E-17,2.27900E-17,1.30240E-17,
     > 2.66100E-17,2.27870E-17,1.30900E-17,2.27900E-17,2.24000E-17,
     > 1.25900E-17,2.60400E-17,2.41300E-17,1.20590E-17,2.46060E-17,
     > 2.45010E-17,1.21270E-17,2.31010E-17,2.34710E-17,1.19300E-17,
     > 2.19100E-17,2.31600E-17,1.14960E-17,2.03100E-17,2.16750E-17,
     > 9.68700E-18,1.81180E-17,1.63950E-17,9.84000E-18,1.83200E-17,
     > 1.69100E-17,8.69300E-18,1.74380E-17,1.38570E-17,7.70000E-18,
     > 1.68100E-17,1.17000E-17,7.68000E-18,1.68000E-17,1.16700E-17,
     > 6.46100E-18,1.43870E-17,1.04930E-17,7.08000E-18,1.57900E-17,
     > 1.09000E-17,6.05000E-18,1.33700E-17,1.02100E-17,5.20200E-18,
     > 1.09000E-17,8.39200E-18,3.73200E-18,7.50900E-18,4.95800E-18,
     > 1.83900E-18,3.80600E-18,2.26100E-18,7.30000E-19,1.31600E-18,
     > 7.20000E-19/
      DATA SIGABS/0.0000,1.34600E-18,0.000,0.00000E+00,1.00000E-18,
     > 0.00000E+00,0.00000E+00,1.63000E-18,0.00000E+00,0.00000E+00,
     > 2.11080E-17,5.09880E-17,0.00000E+00,1.87300E-17,2.24000E-18,
     > 1.31500E-18,1.28170E-17,9.68000E-18,4.55400E-18,8.56200E-18,
     > 2.02490E-17,3.49800E-18,1.66310E-17,1.69920E-17,5.09100E-18,
     > 2.21450E-17,3.35780E-17,3.74900E-18,2.66680E-17,1.64870E-17,
     > 3.89000E-18,1.89100E-17,1.41800E-17,4.00000E-18,2.08000E-17,
     > 1.20490E-16,1.07360E-17,2.85350E-17,2.46620E-17,1.14600E-17,
     > 2.74400E-17,2.65400E-17,1.72450E-17,2.19190E-17,3.17550E-17,
     > 1.33650E-17,2.60170E-17,2.33390E-17,1.34000E-17,3.20600E-17,
     > 2.33700E-17,1.34000E-17,2.80700E-17,2.27900E-17,1.30240E-17,
     > 2.66100E-17,2.27870E-17,1.30900E-17,2.27900E-17,2.24000E-17,
     > 1.25900E-17,2.60400E-17,2.41300E-17,1.20590E-17,2.46060E-17,
     > 2.45010E-17,1.21270E-17,2.31010E-17,2.34710E-17,1.19300E-17,
     > 2.19100E-17,2.31600E-17,1.14960E-17,2.03100E-17,2.16750E-17,
     > 9.68700E-18,1.81180E-17,1.63950E-17,9.84000E-18,1.83200E-17,
     > 1.69100E-17,8.69300E-18,1.74380E-17,1.38570E-17,7.70000E-18,
     > 1.68100E-17,1.17000E-17,7.68000E-18,1.68000E-17,1.16700E-17,
     > 6.46100E-18,1.43870E-17,1.04930E-17,7.08000E-18,1.57900E-17,
     > 1.09000E-17,6.05000E-18,1.33700E-17,1.02100E-17,5.20200E-18,
     > 1.09000E-17,8.39200E-18,3.73200E-18,7.50900E-18,4.95800E-18,
     > 1.83900E-18,3.80600E-18,2.26100E-18,7.30000E-19,1.31600E-18,
     > 7.20000E-19/
       DATA PROB/0.0,1.0,1.0,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,0.957794,1.000000,0.000000,0.042206,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,0.895954,1.000000,0.000000,0.104046,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.627701,0.565000,0.527407,0.372299,0.435000,0.472593,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.426702,0.668824,0.388933,0.573298,0.331176,0.611067,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.297362,0.395273,0.335986,0.664135,0.463591,0.664014,
     > 0.038504,0.141136,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.317995,0.240422,0.308000,0.460831,0.365353,0.589000,
     > 0.221175,0.340922,0.103000,0.000000,0.053304,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.289552,0.242130,0.308000,0.469403,0.359162,0.589000,
     > 0.241045,0.352241,0.103000,0.000000,0.046467,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.269403,0.234886,0.308000,0.460448,0.383177,0.589000,
     > 0.270149,0.306381,0.103000,0.000000,0.075556,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.277948,0.354208,0.332000,0.453778,0.280858,0.568250,
     > 0.268274,0.214870,0.099750,0.000000,0.150064,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.299465,0.303890,0.323043,0.459893,0.328831,0.575994,
     > 0.240642,0.213666,0.100963,0.000000,0.153613,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.289913,0.355053,0.351862,0.450357,0.235358,0.551077,
     > 0.259730,0.216569,0.097060,0.000000,0.193020,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.283689,0.359411,0.380000,0.450037,0.282907,0.526750,
     > 0.266274,0.128627,0.093250,0.000000,0.219594,0.000000,
     > 0.000000,0.009461,0.000000,0.000000,0.000000,0.000000,
     > 0.279871,0.265374,0.405764,0.450153,0.160558,0.469567,
     > 0.269976,0.077351,0.082915,0.000000,0.084977,0.039289,
     > 0.000000,0.411741,0.002465,0.000000,0.000000,0.000000,
     > 0.279966,0.242785,0.414980,0.450126,0.139718,0.465742,
     > 0.269908,0.068087,0.081994,0.000000,0.064217,0.034222,
     > 0.000000,0.485192,0.003062,0.000000,0.000000,0.000000,
     > 0.267310,0.419160,0.439757,0.425887,0.235695,0.446402,
     > 0.259743,0.120885,0.078528,0.047060,0.122618,0.028560,
     > 0.000000,0.101641,0.006687,0.000000,0.000000,0.000066,
     > 0.260452,0.397023,0.361839,0.400846,0.223099,0.480085,
     > 0.250026,0.122412,0.099384,0.088676,0.144374,0.026681,
     > 0.000000,0.113092,0.030068,0.000000,0.000000,0.001944,
     > 0.260061,0.393955,0.346750,0.400000,0.221354,0.478991,
     > 0.259959,0.122624,0.101074,0.079980,0.147389,0.029676,
     > 0.000000,0.114679,0.040767,0.000000,0.000000,0.002741,
     > 0.259864,0.374885,0.267385,0.396411,0.210504,0.462571,
     > 0.249971,0.123939,0.107338,0.093754,0.166130,0.022626,
     > 0.000000,0.124542,0.114657,0.000000,0.000000,0.025423,
     > 0.250000,0.365000,0.256053,0.370000,0.205000,0.440509,
     > 0.250000,0.125000,0.103843,0.090000,0.055006,0.021138,
     > 0.040000,0.249993,0.089033,0.000000,0.000000,0.089423,
     > 0.250000,0.365000,0.255808,0.370052,0.205000,0.440035,
     > 0.250000,0.125000,0.103766,0.089974,0.055006,0.021076,
     > 0.039974,0.249993,0.088240,0.000000,0.000000,0.091075,
     > 0.251973,0.365000,0.235293,0.359851,0.205000,0.404410,
     > 0.230305,0.125000,0.095588,0.098592,0.055006,0.011581,
     > 0.059279,0.249993,0.069486,0.000000,0.000000,0.183642,
     > 0.250000,0.365000,0.242093,0.370056,0.205000,0.416098,
     > 0.230226,0.125000,0.098350,0.100282,0.055006,0.014699,
     > 0.049435,0.249993,0.078089,0.000000,0.000000,0.150670,
     > 0.254380,0.365000,0.221395,0.353388,0.205000,0.380522,
     > 0.225289,0.125000,0.089942,0.098843,0.055006,0.006278,
     > 0.068099,0.249993,0.037670,0.000000,0.000000,0.264193,
     > 0.257513,0.365000,0.207040,0.346101,0.205000,0.355850,
     > 0.224710,0.125000,0.084110,0.099784,0.055006,0.000000,
     > 0.071892,0.249993,0.000000,0.000000,0.000000,0.353000,
     > 0.270685,0.365000,0.204800,0.332954,0.205000,0.352000,
     > 0.218368,0.125000,0.083200,0.098948,0.055006,0.000000,
     > 0.079045,0.249993,0.000000,0.000000,0.000000,0.360000,
     > 0.294011,0.365000,0.204800,0.320024,0.205000,0.352000,
     > 0.208711,0.125000,0.083200,0.098609,0.055006,0.000000,
     > 0.078645,0.249993,0.000000,0.000000,0.000000,0.360000,
     > 0.296412,0.365000,0.204800,0.321373,0.205000,0.352000,
     > 0.209048,0.125000,0.083200,0.096724,0.055006,0.000000,
     > 0.076443,0.249993,0.000000,0.000000,0.000000,0.360000/ 
      END MODULE PHOTOIONIZATION_DATA
C:::::::::::::::::::::::::::::::::: CTIPINT :::::::::::::::::::::::::::::
C.... Interface for CTIP model. This routine uploads the CTIPe variables to 
C.... modules and calls the different FLIP routines.
C.... The parameters ending in X are dummy variables to be transferred to 
C.... the FLIP module variables. 
C.... Dummy variables are also used because it is assumed that the values 
C.... of FLIP variables are not preserved from call to call
C---- Additional mods Change TI to (TI(2,DIM)
C---- Change SCOLUMN for grazing incidence
C.... Written by P. Richards June-September 2010.
      SUBROUTINE CTIPINT(
     >             JMINX,  !.. index of the first point on the field line
     >             JMAXX,  !.. index of the last point on the field line
     >           CTIPDIM,  !.. CTIPe array dimension, must equal to FLDIM
     >                ZX,  !.. array, altitude (km)
     >               PCO,  !.. p coordinate (L-shell)
     >               SLX,  !.. array, distance of point from northern hemisphere
     >               GLX,  !.. array, magnetic latitude (radians)
     >               BMX,  !.. array, magnetic field strength, (Tesla)
     >               GRX,  !.. array, gravity, cm2 s-1
     >                OX,  !.. array, O density (cm-3)
     >                HX,  !.. array, H density (cm-3)
     >               N2X,  !.. array, N2 density (cm-3)
     >               O2X,  !.. array, O2 density (cm-3)
     >               HEX,  !.. array, He density (cm-3)
     >              N4SX,  !.. array, N(4S) density (cm-3)
     >              INNO,  !.. switch to turn on FLIP NO calculation if <0
     >              NNOX,  !.. array, NO density (cm-3)
     >               TNX,  !.. array, Neutral temperature (K)
     >             TINFX,  !.. array, Exospheric Neutral temperature (K)
     >               UNX,  !.. array, Neutral wind (m/s)
     >                DT,  !.. CTIPe time step (secs)
     >             DTMIN,  !.. Minimum time step allowed (>=10 secs?)
     >              F107,  !.. Daily F10.7
     >             F107A,  !.. 81 day average F10.7
     >              SZAX,  !.. Solar Zenith angle (radians)
     >              FPAS,  !.. Pitch angle scattering fraction
     >              HPEQ,  !.. Sets initial equatorial H+ density. See declaration below
     >            HEPRAT,  !.. Intial He+/H+ ratio (.01 to 1.0)
     >           COLFACX,  !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
     >            IHEPLS,  !.. switches He+ diffusive solution on if > 0
     >             INPLS,  !.. switches N+ diffusive solution on if > 0
     >              EHTX,  !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
     >            TE_TIX,  !.. OUT: 2D array, Electron and ion temperatures (K) (see below)
     >     XIONNX,XIONVX,  !.. OUT: 2D array, Storage for ion densities and velocities 
     >             NHEAT,  !.. OUT: array, Neutral heating rate (eV/cm^3/s) 
     >             EFLAG,  !.. OUT: 2D array, Error Flags
     &                mp, 
     &                lp,
     &                utime, !dbg20141209
     &         hrate_cgs, mlt ) !.. OUTPUT: (eV/cm^3/s) !nm20121020

      USE THERMOSPHERE       !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE MINORNEUT          !.. N4S N2D NNO N2P N2A O1D O1S
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE ION_DEN_VEL        !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION         !.. EUV, photoelectron, and auroral production

      USE module_input_parameters,ONLY: sw_TEI,sw_OHPLS,sw_PE2S
     &, sw_DEBUG_flip,sw_debug,sw_output_fort167
     &,mpfort167,lpfort167,mype,peFort167
     &,sw_optw_flip
     &,start_time,ip_freq_output,sw_aurora,LPI,LevPI,GWatts
      USE module_IO,ONLY: LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4
!tmp20151112      USE module_FIELD_LINE_GRID_MKS,ONLY: mlon_rad
!tmp20151112      USE module_physical_constants,ONLY: zero !,pi
!tmp20151112      USE module_magfield,ONLY:sunlons
      IMPLICIT NONE
      include "gptl.inc"
      INTEGER CTIPDIM         !.. CTIPe array dimension, must equal to FLDIM
!nm20110923      INTEGER JTI             !.. Dummy variable to count the number of calls to this routine
      INTEGER I,J,JMINX,JMAXX !.. lcv + spatial grid indices
      INTEGER EFLAG(11,11)    !.. error flags, check =0 on return from FLIP
      INTEGER, INTENT(IN):: mp,lp,utime !dbg20141209
      INTEGER INNO            !.. switch to turn on FLIP NO calculation if <0
!nm20110810      INTEGER DEBUG           !.. switch to turn on debug writes 0=off, 1=on
      !.. IHEPLS,INPLS turn on diffusive solutions if > 0. no solution if 0, 
      !.. chemical equilibrium if < 0
      INTEGER IHEPLS,INPLS  !.. switches He+ and N+ diffusive solutions on
      REAL F107,F107A         !.. daily and 81 day average F10.7      !.. 
      REAL EUVFLUX(37),UVFAC(59)  !.. Solar EUV fluxes and mult. factors
      DOUBLE PRECISION M_to_CM,M3_to_CM3  !.. Unit conversion factors
      !.. Dummy variables ending in X for transferring CTIP to FLIP modules, 
      !.. see above for documentation
      DOUBLE PRECISION N4SX(FLDIM),NNOX(FLDIM),ZX(FLDIM)
      DOUBLE PRECISION UNX(FLDIM),COLFACX,BMX(FLDIM)
      DOUBLE PRECISION GRX(FLDIM),SLX(FLDIM),GLX(FLDIM)
      DOUBLE PRECISION OX(FLDIM),HX(FLDIM),N2X(FLDIM),O2X(FLDIM)
      DOUBLE PRECISION HEX(FLDIM),TNX(FLDIM),SZAX(FLDIM)
      DOUBLE PRECISION XIONNX(9,FLDIM),XIONVX(9,FLDIM)
      !.. TE_TI(3,J) = Te, TE_TIX(2,J) = Ti = TE_TIX(2,J)
      DOUBLE PRECISION TE_TIX(3,FLDIM)
      !.. EHTX(3,J) = e heating rate, EHTX(1,J) = ion heating rate, EHTX(2,J) unused
      DOUBLE PRECISION EHTX(3,FLDIM)
      !.. TINFX has to be an array for grazing incidence column densities
      DOUBLE PRECISION TINFX(FLDIM)  !.. Exospheric temperature
      !.. End dummy variable declarations 
      !.. HPEQ is equatorial H+ density = HPEQ * density of a full flux tube.
      !.. If zero the densities from the previous time step are used.
      !.. If positive, initial densities and temperatures are set. 
      !.. If negative, H+ and He+ densities are reset for flux tube depletion
      DOUBLE PRECISION HPEQ            !.. HPEQ is equatorial H+ density
      DOUBLE PRECISION DT,DTMIN,FD(9),BCKPRD,FPAS,HEPRAT,PCO
      DOUBLE PRECISION COLUM(3,FLDIM)  !.. Neutral column densities for PEPRIM
      DOUBLE PRECISION N(4,FLDIM)      !.. FLIP variable for O+ H+ & total ions
      DOUBLE PRECISION TI(3,FLDIM)     !.. FLIP variable for Te and Ti
      DOUBLE PRECISION NHEAT(FLDIM)    !.. Neutral heating rate
      DOUBLE PRECISION RTS(99)         !.. Reaction rates
      DOUBLE PRECISION EDEN(FLDIM)     !.. electron density
      DOUBLE PRECISION O2DISF(FLDIM)   !.. O2 dissociation frequency
!nm20121020
      DOUBLE PRECISION hrate_cgs(22,FLDIM)   !.. heating rates

      DATA M_TO_CM,M3_TO_CM3/1.0E+2,1.0E-6/    !.. Unit conversion factors
!dbg20110120:      DATA DEBUG/1/  !.. turn on debug writes if DEBUG=1
      INTEGER :: midpoint !nm20110312
      integer :: ret
!nm20150322
      integer :: j802 
!nm20151030 aurora
      integer :: tirosdim
      REAL :: gm_lat 
      REAL,intent(IN) :: mlt
      INTEGER :: tiros_activity_level
      REAL :: gw
      REAL*8,dimension(3,fldim) :: qiont !1:O;2:O2;3:N2 !units number/m3/s

      ret = gptlstart ('CTIPINT init_params')
      CALL initialize_module_parameters ( )
      ret = gptlstop  ('CTIPINT init_params')
!nm20110923      JTI=JTI+1

      !.. Load debug flag into EFLAG for sending to subroutines
      EFLAG(11,11) = sw_DEBUG_flip 

      !.. Check that dimensions agree
      IF(CTIPDIM.GT.FLDIM) THEN
        EFLAG(11,1) =-1
        RETURN
      ENDIF
      ret = gptlstart ('CTIPINT upload')
      !.. Upload field line grid parameters to FIELD_LINE_GRID module
      JMIN=JMINX
      JMAX=JMAXX
      DO J=JMIN,JMAX
        Z(J)=ZX(J)
        SZA(J)=SZAX(J)
        BM(J)=BMX(J)*1.0E+4    !Tesla to gauss 
        SL(J)=SLX(J)*M_to_CM  !nm110210
        GL(J)=GLX(J)
        GR(J)=GRX(J)*M_to_CM  !nm110210
      ENDDO

      !.. Upload densities and velocities to ION_DEN_VEL module
      DO J=JMIN,JMAX
      DO I=1,ISPEC
        XIONN(I,J)=XIONNX(I,J)*M3_to_CM3
        XIONV(I,J)=XIONVX(I,J)*M_to_CM
      ENDDO
      ENDDO
!
      !.. Transfer Te and Ti to FLIP variable TI
      DO J=JMIN,JMAX
      DO I=1,3
        TI(I,J)=TE_TIX(I,J)
!dbg20110712: debug check which output is essential?
        EHT(I,J)=EHTX(I,J)
      ENDDO
      ENDDO

      !.. Upload thermosphere parameters to THERMOSPHERE module
      COLFAC=COLFACX  !.. O+ - O collision frequency Burnside factor 1-1.7 
      DO J=JMIN,JMAX
        ON(J)=OX(J)*M3_to_CM3
        HN(J)=HX(J)*M3_to_CM3
        N2N(J)=N2X(J)*M3_to_CM3
        O2N(J)=O2X(J)*M3_to_CM3
        HE(J)=HEX(J)*M3_to_CM3
        TN(J)=TNX(J)
        TINF(J)=TINFX(J)
        N4S(J)=N4SX(J)*M3_to_CM3
!dbg20110712: debug check which output is essential?
        NNO(J)=NNOX(J)*M3_to_CM3
        UN(J)=UNX(J)*M_to_CM
        !.. transfer densities from storage to FLIP solution variable N
        N(1,J)=XIONN(1,J)
        N(2,J)=XIONN(2,J)
        N(3,J)=(XIONN(4,J)+XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J))
        N(4,J)=XIONN(3,J)
        NHEAT(J)=0.0
        O2DISF(J)=0.0
!nm20121020
        hrate_cgs(1:22,J)=0.0
      ENDDO
      ret = gptlstop ('CTIPINT upload')

!nm20151029 aurora
!nm20151102      if ( mp==17.and.lp==22 ) then !UT= 120.000
      ret = gptlstart ('ionize_ipe')
      if(sw_debug) then
!sms$ignore begin
        print *,'starting ionize_ipe',mype
!sms$ignore end
      endif
      gm_lat = GL(JMIN)
! call tiros only >50deg mlat
      if( abs(gm_lat*57.295779513)>=50.0 .AND. sw_aurora==1 ) then
!tmp20151112      mlt = mlon_rad(mp)*180./PI/15.0D0-sunlons(1)*12.0D0/PI+12.0 !;[hr]
!      mlt = 0.168238 ![hr]
         midpoint = (JMAX/2)+1
         tirosdim = midpoint-JMIN+1
         if(sw_debug) then
!sms$ignore begin
           print *,jmax,jmin,tirosdim, midpoint,' GL1=',gm_lat
     &        ,mlt
     &  , maxval(qiont(1,JMIN:midpoint)), minval(qiont(1,JMIN:midpoint))
!sms$ignore end
         endif
         tiros_activity_level = LevPI(LPI)
         gw = GWatts(LPI)
         if ( mp==1.and.lp==1 ) write(unit=1003,FMT='(I8,I3,f7.2)')
     &utime,tiros_activity_level,gw
         call IONIZE_IPE (
     &     tirosdim,z(JMIN:midpoint),gr(JMIN:midpoint),on(JMIN:midpoint)
     &        ,o2n(JMIN:midpoint)
     &        ,n2n(JMIN:midpoint),hn(JMIN:midpoint),he(JMIN:midpoint)
     &        ,tn(JMIN:midpoint)
     &        ,gm_lat,mlt
     &        ,tiros_activity_level,gw
     &        ,qiont(1,JMIN:midpoint),qiont(2,JMIN:midpoint) !units number/m3/s
     &        ,qiont(3,JMIN:midpoint) 
     &        ,sw_debug)
         if(sw_debug) then
!sms$ignore begin
           print *,'stoping ionize_ipe',mype
!sms$ignore end
         endif
      else 
         qiont=0.0
      endif                     !sw_aurora
      ret = gptlstop ('ionize_ipe')
 
!

      !.. Set up initial temperature and density profiles.
      !.. 0.1 < HPEQ < 1.0.
      ret = gptlstart ('CTIPINT PROFIN')
      IF(HPEQ.GT.0.001)
     >  CALL PROFIN(IHEPLS,INPLS,PCO,F107,N,TI,HPEQ,HEPRAT)
      ret = gptlstop  ('CTIPINT PROFIN')
      !.. This routine adjusts the H+ and He+ densities for depleted flux tubes      
      !..  if HPEQ is negative. 0.1 < -HPEQ < 1.0
      ret = gptlstart ('CTIPINT NEW_HP')
      IF(HPEQ.LT.-0.001) CALL NEW_HP(JMIN,JMAX,PCO,HPEQ,N,EFLAG)
      ret = gptlstop  ('CTIPINT NEW_HP')

      !.. Update solar EUV flux factors
      ret = gptlstart ('CTIPINT FACEUV')
      CALL FACEUV(F107,F107A,UVFAC,EUVFLUX)
      ret = gptlstop  ('CTIPINT FACEUV')

      !.. Update Schumann-Runge UV flux factors
      ret = gptlstart ('CTIPINT FACSR')
      CALL FACSR(UVFAC,F107,F107A)
      ret = gptlstop  ('CTIPINT FACSR')

      !.... evaluate primary EUV production
      ret = gptlstart ('CTIPINT PRIMPR')
      CALL PRIMPR(F107,F107A,UVFAC,COLUM,EUVFLUX)
      ret = gptlstop  ('CTIPINT PRIMPR')

      !.. electron density for photoelectron routine
      DO J=JMIN,JMAX
        EDEN(J)=XIONN(1,J)+XIONN(2,J)+XIONN(3,J)+XIONN(4,J)+XIONN(5,J)+
     >    XIONN(6,J)
      ENDDO

      !.. 2-stream photoelectron routine to get electron heating 
      !.. rate and secondary ion production
      ret = gptlstart ('CTIPINT PE2S')
      IF( sw_PE2S>0 )  !dbg20141210
     &   CALL PE2S(F107,F107A,N,TI,FPAS,-1.0E22,EDEN,UVFAC,COLUM,
     > IHEPLS,INPLS,INNO
     &,mp,lp,utime) !dbg20141209
!dbg      sw_PE2S=0  !dbg20141210.v4
      ret = gptlstop  ('CTIPINT PE2S')

      !-- Sum the EUV, photoelectron, and auroral production rate
      ret = gptlstart ('CTIPINT SUMPRD')
      CALL SUMPRD(JMIN,JMAX
     &,qiont(1:3,JMIN:JMAX) !units number/m3/s
     &)
      ret = gptlstop  ('CTIPINT SUMPRD')

      !.. Loop to calculate O+(4S) total ionization rate
      !.. PHION=total O+(4S) prod, including EUV, e*, dissoc of O2 and
      !.. minor ions. BCKPRD = small background production to eliminate 
      !,, instability below 200 km
      ret = gptlstart ('CTIPINT tot ion')
      DO J=JMIN,JMAX 
         PHION(J)=1.0E-22
         N(3,J)=1.0E-22
         DO I=1,9
            FD(I)=0.0
         ENDDO

         IF(Z(J).GE.80.AND.Z(J).LE.700) THEN
            !.. CALL cminor to get NO+, O+(2D), O2+, N2+ & O+(2P) densities
            CALL CMINOR(0,J,0,IHEPLS,INPLS,INNO,FD,7,N,TI,Z,EFLAG)
            BCKPRD=2.0E-10*N2N(J)*EXP(-5.0E-14*N2N(J)*TN(J))
            PHION(J)=SUMION(1,7,J)+SUMION(2,4,J)+SUMION(2,5,J)+FD(9)
            PHION(J)=PHION(J)+BCKPRD

         ELSE

           !.. Make sure minor ions and neutrals are zero above upper boundary
           DO I=5,ISPEC
             XIONN(I,J)=0
           ENDDO
           N2D(J)=0.0
           N2A(J)=0.0
           O1D(J)=0.0
           O1S(J)=0.0
           N2P(J)=0.0
           NNO(J)=0.0
         ENDIF
         !.. Sum minor ions N+, NO+, O2+, N2+ for electron density at low altitudes
         N(3,J)=XIONN(4,J)+XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J)
      ENDDO
      ret = gptlstop ('CTIPINT tot ion')

      !.. Debug write
      ret = gptlstart ('CTIPINT sw_output')
      IF (mype==peFort167.AND.sw_output_fort167.AND.mp==mpfort167.AND.  &
     &lp==lpfort167 ) THEN
!sms$ignore begin
        print*,'check unit#',LUN_FLIP1,LUN_FLIP3,LUN_FLIP2,LUN_FLIP4,   &
     &                       mype
!sms$ignore end
c      IF(JTI.EQ.1) THEN
        WRITE(UNIT=LUN_FLIP1,FMT=201)  
        WRITE(UNIT=LUN_FLIP3,FMT=201)
 201    FORMAT('   JMIN   JMAX   CTIPDIM  INNO  IHEPLS  INPLS')
        WRITE(UNIT=LUN_FLIP1,FMT='(22I7)')  JMIN,JMAX,CTIPDIM,INNO
     &,IHEPLS,INPLS

        WRITE(UNIT=LUN_FLIP1,FMT=202)  
        WRITE(UNIT=LUN_FLIP3,FMT=202)
 202    FORMAT(/8X,'PCO       DT          DTMIN       F107       F107A'
     >    ,10X,'FPAS        HPEQ      HEPRAT      COLFACX')
        WRITE(UNIT=LUN_FLIP1,FMT='(22F12.3)')  PCO,DT,DTMIN,F107,F107A
     &,FPAS,HPEQ,HEPRAT,COLFACX

        WRITE(UNIT=LUN_FLIP1,FMT=203)
        WRITE(UNIT=LUN_FLIP3,FMT=203)
 203    FORMAT(/5X,'Z      SL            GL      BM        GR '
     >   ,6X,'SZA        O         H       N2       O2       HE '
     >   ,6X,'N4S ')

        WRITE(UNIT=LUN_FLIP2,FMT=204)
        WRITE(UNIT=LUN_FLIP4,FMT=204)
 204    FORMAT(5X,'Z         TN       UN       NNO'
     >    ,6X,'EHT      TI       TE       O+       H+      Min+'
     >    ,5X,'He+      PHION    PRODO+     N+     EQN2D   NPLSPRD')

        !.. Northern Hemisphere
        DO J=JMIN,(JMAX/2)+1
        N(4,J)=XIONN(3,J)
        WRITE(UNIT=LUN_FLIP1,FMT='(F10.2,1P,E14.7,21E9.2)') Z(J),SL(J)
     &,GL(J),BM(J),GR(J),SZA(J),ON(J),HN(J),N2N(J),O2N(J),HE(J),N4S(J)

          WRITE(UNIT=LUN_FLIP2,FMT='(3F10.2,1P,22E9.2)') Z(J),TNX(J)
     &,UN(J),NNO(J)
     >    ,EHT(3,J),TI(1,J),TI(3,J),N(1,J),N(2,J),N(3,J),XIONN(3,J)
     >    ,PHION(J),SUMION(1,7,J),XIONN(4,J),EQN2D(J),NPLSPRD(J) 
!        WRITE(168,'(3F10.2,1P,9E9.2,E10.2,3E9.2)') Z(J),TNX(J),UN(J)
!     &,NNO(J)
!     >   ,EHT(3,J),TI(1,J),TI(3,J),N(1,J),N(2,J),N(3,J),N(4,J)
!     >   ,PHION(J)
!dbg20110404
!     &   ,XIONV(1,J)
!     &,SUMION(1,7,J),SUMION(2,4,J),SUMION(2,5,J)
        ENDDO

        !.. Southern Hemisphere
        DO J=JMAX,(JMAX/2)+1,-1
        N(4,J)=XIONN(3,J)
        WRITE(UNIT=LUN_FLIP3,FMT='(F10.2,1P,E14.7,21E9.2)') Z(J),SL(J)
     &,GL(J),BM(J),GR(J),SZA(J),ON(J),HN(J),N2N(J),O2N(J),HE(J),N4S(J)

          WRITE(UNIT=LUN_FLIP4,FMT='(3F10.2,1P,22E9.2)') Z(J),TNX(J)
     &,UN(J),NNO(J)
     >    ,EHT(3,J),TI(1,J),TI(3,J),N(1,J),N(2,J),N(3,J),XIONN(3,J)
     >    ,PHION(J),SUMION(1,7,J),XIONN(4,J),EQN2D(J),NPLSPRD(J) 
!        WRITE(171,'(3F10.2,1P,9E9.2,E10.2,3E9.2)') Z(J),TNX(J),UN(J)
!     &,NNO(J)
!     >   ,EHT(3,J),TI(1,J),TI(3,J),N(1,J),N(2,J),N(3,J),N(4,J)
!     >   ,PHION(J)
!!dbg20110404
!     &   ,XIONV(1,J)
!     &,SUMION(1,7,J),SUMION(2,4,J),SUMION(2,5,J)
        ENDDO
      END IF                    !( sw_output_fort167.AND...

      ret = gptlstop ('CTIPINT sw_output')
c      ENDIF

      !..  electron and ion temperature solution
      midpoint = (JMAX/2)+1
      ret = gptlstart ('CTIPINT TLOOPS')
      IF( sw_TEI>0) ! .AND. Z(midpoint)>100.00 )
     >  CALL TLOOPS(JMIN,JMAX,FLDIM,Z,N,TI,DT,DTMIN,EFLAG)   !$$$ 
      ret = gptlstop  ('CTIPINT TLOOPS')
      !.. O+, H+ solution
      ret = gptlstart ('CTIPINT DLOOPS')
      IF( sw_OHPLS>0) ! .AND. Z(midpoint)>120.00 )
     >  CALL DLOOPS(JMIN,JMAX,FLDIM,Z,N,TI,DT,DTMIN,EFLAG)   !$$$  
      ret = gptlstop  ('CTIPINT DLOOPS')
!----------------------
      if ( sw_optw_flip ) then
      !.. Recalculate the minor ion densities with the new O+ density
      !.. Added by PGR 2012-11-29
      DO J=JMIN,JMAX
       IF(Z(J).GE.80.AND.Z(J).LE.700) THEN
          CALL CMINOR(0,J,0,IHEPLS,INPLS,INNO,FD,7,N,TI,Z,EFLAG)
       ENDIF
      ENDDO
      endif !( sw_optw_flip ) then
!----------------------
      !.. He+ solution
      ret = gptlstart ('CTIPINT XION')
      IF(EFLAG(2,1).EQ.0.AND.IHEPLS.GT.0) ! .AND. Z(midpoint)>200.00 )
     & CALL XION(TI,DT,DTMIN,9,EFLAG)
      ret = gptlstop  ('CTIPINT XION')

      !.. N+ solution
!dbg20120301:
      IF ( sw_DEBUG_flip==1 ) then
!sms$ignore begin
         print *,'!dbg! apex ht=',z(midpoint),midpoint,lp,mp,mype
!sms$ignore end
      endif
      ret = gptlstart ('CTIPINT XION')
      IF(EFLAG(2,1).EQ.0.AND.INPLS.GT.0) CALL XION(TI,DT,DTMIN,11,EFLAG)
      ret = gptlstop  ('CTIPINT XION')

!nm20150322: output limiting flux: not so helpful for SED plume flux...
!
!      IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN
!      j802=87 !z_km=802km
!      if ( lp <=40 )   !mlat=38.85
!     & write(9009,"(2i3,9E9.2)") mp,lp, z(j802),on(j802),hn(j802)
!     &,tnx(j802) ,ti(1,j802),ti(3,j802),n(1,j802),n(2,j802)
!     &, (1.38E-16 * ( ti(1,j802)+ti(3,j802) )*0.5 / 
!     &   (1.662E-24*16*(-1.)*GR(j802) ) )*1.E-5 !scale height[cm-->km]
!      ENDIF

      ret = gptlstart ('CTIPINT transfer')
        !.. transfer densities from FLIP to CTIP variable
      DO J=JMIN,JMAX
        NNOX(J)=NNO(J)/M3_to_CM3
        DO I=1,ISPEC
          XIONNX(I,J)=XIONN(I,J)/M3_to_CM3
          XIONVX(I,J)=XIONV(I,J)/M_to_CM
        ENDDO
      ENDDO

      !.. Transfer Te and Ti and e heating to CTIPe variable
      DO J=JMIN,JMAX
      DO I=1,3
        TE_TIX(I,J)=TI(I,J)
        EHTX(I,J)=EHT(I,J)
      ENDDO
      ENDDO
      ret = gptlstop ('CTIPINT transfer')

!.. I=1 is used to turn on NHEAT writes
!nm033111:      IF(DEBUG.EQ.1.AND.JTI.EQ.115) I=1 
      IF ( sw_DEBUG_flip==1 ) THEN 
        I=1 !ON write
      ELSE
        I=0 !OFF
      END IF

        !.. Get neutral gas heating rate NHEAT
        ret = gptlstart ('CTIPINT get NHEAT')
        DO J=JMIN,JMAX
          IF(Z(J).GE.80.AND.Z(J).LE.700) THEN
            !.. electron density for photoelectron routine
            EDEN(J)=XIONN(1,J)+XIONN(2,J)+XIONN(3,J)+XIONN(4,J)+
     >      XIONN(5,J)+XIONN(6,J)
            CALL RATS(J,TI(3,J),TI(2,J),TN(J),RTS)  !.. Reaction rates
            !.. Neutral heating rate
!dbg      print *,'before call NEUT_HEATING: I=',I
            CALL NEUT_HEATING(I,JMIN,JMAX
     &       ,J,Z(J),RTS,TI(3,J),TI(2,J),TN(J),
     >        ON(J),O2N(J),N2N(J),HE(J),N4S(J),EDEN(J),N(1,J),XIONN(6,J)
     >       ,XIONN(5,J),N(2,J),XIONN(7,J),XIONN(4,J),NNO(J),N2D(J)
     >       ,N2P(J),N2A(J),XIONN(8,J),XIONN(9,J),O1D(J),O1S(J)
     >       ,EHT(3,J)   !.. Input: Electron heating rate
     &       ,NHEAT(J)   !.. OUTPUT: Total neutral heating rate
     &       ,O2DISF(J)  !.. OUTPUT: O2 dissociation frequency !PGR added index
     &,hrate_cgs(1:22,J) ) !.. OUTPUT: !nm20121020
!     &       ,hrate_cgs(2,J)  !..
!     &       ,hrate_cgs(3,J)  !..
!     &       ,hrate_cgs(4,J)  !..
!     &       ,hrate_cgs(5,J)  !..
!     &       ,hrate_cgs(6,J)  !..
!     &       ,hrate_cgs(7,J) )!..
          ENDIF
        ENDDO
        ret = gptlstop ('CTIPINT get NHEAT')
c      ENDIF

      !---------------------- DEBUG WRITE -----------------------------------
      !.. Debug write if DEBUG=1
      IF ( sw_DEBUG_flip==1 ) WRITE(172,88)    !nm20110923
 88   FORMAT('  J     ALT    TE     TI       O+        H+       He+'
     >  ,8X,'N+        NO+       O2+       NE      UVOX+    PEOX+'
     >  ,6X,'UVN2+      PEN2+     EHT       VO+       VH+      VHe+'
     >  ,7X,'VN+        NO       N2D       N4S      O1D       O1S'
     >  ,6X,'NHEAT')
      DO J=JMIN,JMAX
        IF(sw_DEBUG_flip.EQ.1)    !nm20110923
     >   WRITE(172,'(I5,I7,2F7.1,1P,29E10.2)') 
     >   J,NINT(Z(J)),TI(3,J),TI(2,J),XIONN(1,J),XIONN(2,J),XIONN(3,J),
     >   XIONN(4,J),XIONN(5,J),XIONN(6,J),EDEN(J),EUVION(1,1,J),
     >   PEPION(1,1,J),EUVION(1,1,J),PEPION(3,1,J),EHT(3,J),XIONV(1,J),
     >   XIONV(2,J),XIONV(3,J),XIONV(4,J),NNO(J),N2D(J),N4S(J),O1D(J),
     >   O1S(J),NHEAT(J)
      ENDDO

!nm20110715: for diagnostics only
!.. 2-stream photoelectron routine called to print fluxes
      ret = gptlstart ('CTIPINT PE2S')
      IF(sw_DEBUG_flip.EQ.1.AND.  !nm20110923
     &   sw_PE2S>0 )  !dbg20141210
     &  CALL PE2S(F107,F107A,N,TI,FPAS,300.0,EDEN,UVFAC,COLUM
     &    ,IHEPLS,INPLS,INNO)
      ret = gptlstop  ('CTIPINT PE2S')

      RETURN
      END
C:::::::::::::::::: WRITE_EFLAG ::::::::::::::::::::::::
C... This routine prints the information about error flags
C... Written by P. Richards September 2010
      SUBROUTINE WRITE_EFLAG(PRUNIT,   !.. Unit number to print results
     >                        EFLAG,   !.. Error flag array
     >                           mp, 
     >                           lp,utime,ltime)
      USE module_input_parameters,ONLY:sw_output_fort167,sw_ERSTOP_flip &
     &                                ,mype
      USE module_precision
      IMPLICIT NONE
      INTEGER PRUNIT,EFLAG(11,11)         !.. error flags
      INTEGER (KIND=int_prec),  INTENT(IN) :: mp,lp
      INTEGER (KIND=int_prec),  INTENT(IN) :: utime !universal time [sec]
      REAL    (KIND=real_prec), INTENT(IN) :: ltime !local time [hour]
!
      IF(EFLAG(1,1).NE.0) THEN
        WRITE(PRUNIT,11)mp,lp,mype
!(11)
!t        IF ( sw_ERSTOP_flip==1 )  STOP
      END IF
 11   FORMAT(/'  Convergence failure in Temperature solution (TLOOPS).'
     >  ,2X,'Time step less than minimum.mp=',i4,'lp=',i4,i7)
      IF(EFLAG(1,2).NE.0) then
!dbg110210:
        WRITE(PRUNIT,*)'EFLAG(1,2)',EFLAG(1,2),mype
        WRITE(PRUNIT,12)mp,lp
        IF ( sw_ERSTOP_flip==1 ) THEN
!sms$ignore begin
          print*,"(12)ERSTOP FLIP",mp,lp,mype
!sms$ignore end
          STOP
        END IF !( sw_ERSTOP_flip==1 )
      end IF !(EFLAG(1,2).NE.0)
 12   FORMAT(/'  Convergence failure in Temperature solution (TLOOPS).'
     >  ,2X,'Incorrect input to the band solver BDSLV.mp=',i3,'lp=',i4)     

      IF(EFLAG(2,1).NE.0) THEN
         WRITE(PRUNIT,21)lp,mp,ltime,UTIME,mype
!         print*,'JFM',prunit,mype,lp,mp,ltime,UTIME,mype
!(3)
!t         IF ( sw_ERSTOP_flip==1 )  STOP

      END IF
 21   FORMAT(/'  Convergence failure in O+ - H+ (DLOOPS).'
     &,'TimeStep less than minimum:lp=',i3,'mp=',i3,f7.2,2i7)
      IF(EFLAG(2,2).NE.0) THEN
         WRITE(PRUNIT,22)mp,lp,mype
         IF ( sw_ERSTOP_flip==1 ) THEN
!sms$ignore begin
           print *,"(22)ERSTOP FLIP",mp,lp,mype
!sms$ignore end
           STOP
        END IF


      END IF
 22   FORMAT(/'  Convergence failure in O+ - H+ solution (DLOOPS).'
     >  ,2X,'Incorrect input to the band solver BDSLV.mp',i3,'lp',i4,i6)

      IF(EFLAG(3,1).NE.0) THEN
         WRITE(PRUNIT,31)mp,lp,mype
!(5)
!t         IF ( sw_ERSTOP_flip==1 )  STOP
      END IF
 31   FORMAT(/'  Convergence failure in He+ solution (XION).'
     >  ,2X,'Time step less than minimum.',3i7)
      IF(EFLAG(3,2).NE.0) THEN
         WRITE(PRUNIT,32)mp,lp,mype
         IF ( sw_ERSTOP_flip==1 ) THEN
!sms$ignore begin
           print*,"(32)ERSTOP FLIP",mp,lp,mype
!sms$ignore end
           STOP
         END IF !( sw_ERSTOP_flip==1 ) THEN

      END IF !(EFLAG(3,2).NE.0) THEN
 32   FORMAT(/'  Convergence failure in He+ solution (XION).'
     >  ,2X,'Incorrect input to the band solver BDSLV.mp',i3,'lp',i4,i7)

      IF(EFLAG(4,1).NE.0) THEN
         WRITE(PRUNIT,41)mp,lp,mype
!t         IF ( sw_ERSTOP_flip==1 ) THEN
!t           print *,"(7)ERSTOP FLIP",mp,lp            
!t           STOP
!t        END IF
      END IF !      IF(EFLAG(4,1).NE.0) THEN
 41   FORMAT(/'  Convergence failure in N+ solution (XION).'
     >  ,2X,'Time step less than minimum.mp',i3,'lp',i4,i7)
      IF(EFLAG(4,2).NE.0) THEN
         WRITE(PRUNIT,42)mp,lp,mype
!(8)
         IF ( sw_ERSTOP_flip==1 )  THEN
!SMS$ignore begin
           print*,"(42)ERSTOP FLIP",mp,lp,mype
!SMS$ignore end
           STOP
         END IF !( sw_ERSTOP_flip==1 )  THEN

      END IF !(EFLAG(4,2).NE.0) THEN
 42   FORMAT(/'  Convergence failure in N+ solution (XION).'
     >  ,2X,'Incorrect input to the band solver BDSLV.mp',i3,'lp',i4,i7)

      IF(EFLAG(5,1).NE.0) WRITE(PRUNIT,51)mype
 51   FORMAT(/'  Convergence failure in CMINOR.',i10)

      IF(EFLAG(11,1).NE.0) WRITE(PRUNIT,111)mype
 111  FORMAT(/3X,'** CTIP dimension not equal to FLIP dimensions'
     >  /3X,'** Check dimensions in all FLIP modules',i7)
      RETURN
      END 

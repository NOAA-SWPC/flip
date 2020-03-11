! purpose: initialize flip module parameters
      SUBROUTINE initialize_module_parameters ( )
      USE module_precision
      USE module_physical_constants,ONLY: zero
! module parameters in CTIP-int
      USE FIELD_LINE_GRID, ONLY:Z,SZA,BM,SL,GR,GL
      USE ION_DEN_VEL, ONLY:XIONN,XIONV
      USE SOLVARR,ONLY:DELTA,RHS,WORK,S
      USE THERMOSPHERE, ONLY:ON,HN,N2N,O2N,HE,TN,UN,TINF,EHT,COLFAC
      USE AVE_PARAMS,ONLY:TEJ,TIJ,NUX,RBM,UNJ,DS,GRADTE
     &,GRADTI,GRAV,OLOSS,HLOSS,HPROD
      USE PRODUCTION,ONLY:EUVION,PEXCIT,PEPION,OTHPR1
     &,OTHPR2,SUMION,SUMEXC,PAUION,PAUEXC,NPLSPRD,PHION
      USE MINORNEUT,ONLY:N2D,N2P,N2A,O1D,O1S,EQN2D

!      REAL (KIND=real_prec), PARAMETER :: zero = 0.0_real_prec


!FIELD_LINE_GRID
        Z(:)=zero               !(FLDIM)
        SZA(:)=zero
        BM(:)=zero
        SL(:)=zero
        GR(:)=zero
        GL(:)=zero

!ION_DEN_VEL
        XIONN(:,:)=zero!(ISPEC,IDIM)
        XIONV(:,:)=zero

!(2) SOLVARR
        DELTA(:)=zero !(10*SDIM)
        RHS(:)=zero   !(10*SDIM)
        WORK(:)=zero   !(50*SDIM)
        S(:)=zero   !(50*SDIM)

!THERMOSPHERE
        ON(:)=zero              !(TDIM)
        HN(:)=zero
        N2N(:)=zero
        O2N(:)=zero
        HE(:)=zero
        TN(:)=zero
        UN(:)=zero
        TINF(:)=zero
        EHT(:,:)=zero           !(3,TDIM)
        COLFAC=zero 

!(3)AVE_PARAMS
      TEJ(:)=zero !(AVDIM)
      TIJ(:)=zero !(AVDIM)
      NUX(:,:)=zero !(2,AVDIM)
      RBM(:)=zero !(AVDIM)
      UNJ(:)=zero !(AVDIM)
      DS(:)=zero !(AVDIM)
      GRADTE(:)=zero !(AVDIM)
      GRADTI(:)=zero !(AVDIM)
      GRAV(:)=zero !(AVDIM)
      OLOSS(:)=zero !(AVDIM)
      HLOSS(:)=zero !(AVDIM)
      HPROD(:)=zero !(AVDIM)

!(1) production
        EUVION(:,:,:)=zero
	PEXCIT(:,:,:)=zero
	PEPION(:,:,:)=zero
	OTHPR1(:,:)=zero
	OTHPR2(:,:)=zero
	SUMION(:,:,:)=zero !(3,12,PDIM)
	SUMEXC(:,:,:)=zero !(3,12,PDIM)
	PAUION(:,:,:)=zero !(3,12,PDIM)
	PAUEXC(:,:,:)=zero !(3,12,PDIM)
	NPLSPRD(:)=zero !(PDIM)
	PHION(:)=zero !(PDIM)


!(4)MINORNEUT
        N2D(:)=zero !(NDIM)
	N2P(:)=zero !(NDIM),
	N2A(:)=zero !(NDIM)
	O1D(:)=zero !(NDIM)
	O1S(:)=zero !(NDIM)
	EQN2D(:)=zero !(NDIM)


      END       SUBROUTINE initialize_module_parameters

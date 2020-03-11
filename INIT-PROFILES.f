C........................ INIT-PROFILES.FOR.................................
C.... This file contains routines for setting up initial profiles and 
C.... Adjusting the H+ He+ depleted flux tube profiles
C::::::::::::::::::::::::::::: SET_HPEQ ::::::::::::::::::::::::::::::::
C... This function sets the default H+ density at the equator. It is used in both
C... setting up the initial profiles at the beginning of a run and also for
C... depleting flux tubes. 
C... Written by P. Richards Summer 2010
      SUBROUTINE SET_HPEQ(PCO,  !.. L-value 
     >                   HPEQ)  !.. Fraction of a full flux tube (see below)
      IMPLICIT NONE
      DOUBLE PRECISION PCO,HPEQ,HP_FULL
 
      !.. HPEQ enters as a fraction of a full flux tube but returns the
      !.. equatorial H+ density
      !If HPEQ out of range use a full flux tube
      IF(DABS(HPEQ).LT.0.1.OR.DABS(HPEQ).GT.1.0) HPEQ=1.0
      !... set default HPEQ for a full flux tube
      HP_FULL=5.0E4/PCO**3
      !... Take fraction HPEQ of HP_FULL
      HPEQ=DABS(HPEQ)*HP_FULL
      IF(HPEQ.LT.5) HPEQ=5   !.. don't let density get below 5 
      RETURN
      END
C::::::::::::::::::::::::::::::: PROFIN ::::::::::::::::::::::::::::::::::
C....... set up rough initial O+, H+, and temperature profiles
      SUBROUTINE PROFIN(IHEPLS,INPLS,PCO,F107,N,TI,HPEQ,HEPRAT)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P) ISPEC
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE module_input_parameters,ONLY: init_Te_max
      IMPLICIT NONE
      INTEGER J,JEQ
      INTEGER IHEPLS,INPLS  !.. switches He+ and N+ diffusive solutions on if > 0
      REAL F107
      DOUBLE PRECISION ALT,N(4,IDIM),TI(3,IDIM),HPEQ,XHOLD,HEPRAT,PCO

      CALL SET_HPEQ(PCO,HPEQ)     !.. Set the equatorial [H+]

      !..... Temperatures approximated by a log function
      DO J=JMIN,JMAX
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        TI(3,J)=904.0*DLOG(ALT)-3329.0
!nm20110630: The temperature needs to be limited on large flux tubes, otherwise they will stay too high. The large volume tubes will have too much energy and it will take too long to drain. 
!dbg        IF(TI(3,J).GT.5000) TI(3,J)=5000
!dbg        TI(1,J)=TI(3,J)
!dbg        TI(2,J)=TI(3,J)
        IF(TI(3,J) > init_Te_max) TI(3,J) = init_Te_max
        TI(1,J) = ( TN(J) + TI(3,J) ) * 0.50  !dbg20110810
        TI(2,J) = TI(1,J)  !dbg20110810
      ENDDO
      !.. Northern hemisphere O+
      JEQ=(JMAX+1)/2
      DO J=JMIN,JEQ
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        !...... Insert initial O+ profile: first low altitudes
        XHOLD=-5.5E-4*(ALT-250.0)**2
        IF(ALT.LE.500) N(1,J)=(F107/74.0)*5E5*EXP(XHOLD)
        !....... now high altitudes

!dbg20110810
!       if ( J > 190 .and. J < 210 )  
!     & print "('j',i6,' ALT=',F10.2,10E12.4)", j,alt,xhold,SL(J+1),SL(J)
!     &,GR(J),TI(3,J),TI(1,J)

        IF(ALT.GE.275) THEN 
          XHOLD=(SL(J)-SL(J-1))*1.9E-7*GR(J)/(TI(3,J)+TI(1,J))
          N(1,J)=N(1,J-1)*EXP(XHOLD)
        END IF 
!        if ( alt > 2500. .and. alt < 3500.) 
!       if ( J > 190 .and. J < 210 )  
!     & print "(10E12.4)",N(1,J),N(1,J-1),EXP(XHOLD)


      ENDDO

      !.. Southern hemisphere O+ 
      DO J=JMAX,JEQ+1,-1
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        !...... Insert initial O+ profile: first low altitudes
        XHOLD=-5.5E-4*(ALT-250.0)**2
        IF(ALT.LE.500) N(1,J)=(F107/74.0)*5E5*EXP(XHOLD)
        !....... now high altitudes
        IF(ALT.GE.275) THEN
          XHOLD=-(SL(J+1)-SL(J))*1.9E-7*GR(J)/(TI(3,J)+TI(1,J))
          N(1,J)=N(1,J+1)*EXP(XHOLD)
        END IF
      ENDDO

      !.. Now calculate H+, He+, and N+
      DO J=JMIN,JMAX
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        !...... Initial H+ profile 
        IF(ON(J).GT.0) N(2,J)=N(1,J)*HN(J)/ON(J)
        IF(N(2,J).GT.HPEQ.OR.Z(J).GE.2000) N(2,J)=HPEQ
        !-- He+ density is fraction of H+
        XIONN(3,J) = 0.0
        IF(IHEPLS.GT.0) XIONN(3,J)=HEPRAT*N(2,J)
        !-- N+ density is 10% of O+
        XIONN(4,J) = 0.0
        IF(INPLS.GT.0) XIONN(4,J) = 0.1 * N(1,J)
        XIONN(1,J)=N(1,J)
        XIONN(2,J)=N(2,J)
      ENDDO
      HPEQ=0.0   !.. switches off initial profiles
      RETURN
      END
C:::::::::::::::::: NEW_HP ::::::::::::::::::::::::
C... This routine adjusts the H+ and He+ density for storm depletions.
C... O+ and N+ are unchanged.
C... Written by P. Richards September 2010
      SUBROUTINE NEW_HP(JMIN,     !... First index on field line
     >                  JMAX,     !... Last index on field line
     >                   PCO,     !... L-shell
     >                  HPEQ,     !... New H+ value
     >                     N,     !... ion densities
     >                 EFLAG)     !.. Error flag array
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER J,JMIN,JMAX,JEQ
      INTEGER EFLAG(11,11)         !.. error flags
      DOUBLE PRECISION ALPHA,HP_MIN,HPEQ,PCO,N(4,IDIM)

      !.. equatorial density for a depleted flux tube. Not less than 30
      CALL SET_HPEQ(PCO,HPEQ)     !.. Set the equatorial [H+]

      !.. Equatorial density for a depleted flux tube (~10% full).
      HP_MIN=5000/PCO**3
      IF(HP_MIN.LT.30) HP_MIN=5    !.. don't let density get below 5 

      JEQ=(JMIN+JMAX)/2  !.. equatorial grid point

      !.. No need to deplete already depleted flux tube if equatorial density < min
      IF(N(2,JEQ).LT.HP_MIN) HPEQ=0.0
      IF(N(2,JEQ).LT.HP_MIN) RETURN
      IF(HPEQ.LT.HP_MIN) HPEQ=HP_MIN

      ALPHA=HPEQ/N(2,JEQ)   !.. reduction factor at equator

      !.. Don't delete an already depleted flux tube
      IF(ALPHA.GT.0.8) HPEQ=0.0
      IF(ALPHA.GT.0.8) RETURN

      !.. Reduce H+ and He+ densities. 
      DO J=JMIN,JMAX
        N(2,J)=N(2,J)*ALPHA
        XIONN(3,J)=XIONN(3,J)*ALPHA
      ENDDO
      !.. Debug write if EFLAG(11,11)=1
      IF(EFLAG(11,11).EQ.1) WRITE(6,661) N(2,(JMAX+1)/2),ALPHA,HPEQ
 661  FORMAT(' H+ and He+ reduced; new equatorial [H+]=',F10.2
     >   ,1X,' for large Kp, ALPHA=',F8.1,' HPEQ=',F8.1)
      HPEQ=0.0   !.. switches off reduced densities
      RETURN
      END 

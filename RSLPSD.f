C.................... RSLPSD.FOR ........................................... 
C::::::::::::::::::::::::::: CTIPFDEN :::::::::::::::::::::::::::
C..... This is the interface between CTIP and the FLIP model Te
      SUBROUTINE CTIPFDEN(FLDIM,JMIN,JMAX,Z,N,TI,DTIN,DTMIN,SL,
     >   EFLAG)
      IMPLICIT NONE
      INTEGER JMIN,JMAX,EFLAG(11,11),FLDIM
      REAL SEC    !.. UT (secs) may be used for debug purposes
      DOUBLE PRECISION Z(FLDIM),N(4,FLDIM),TI(3,FLDIM),SL(FLDIM)
      DOUBLE PRECISION DTIN,DTMIN  !.. Initial and minimum time steps

      CALL DLOOPS(JMIN,JMAX,FLDIM,Z,N,TI,DTIN,DTMIN,SL,EFLAG)
      RETURN
      END
C:::::::::::::::::::::::::::::: DLOOPS :::::::::::::::::::::::::::::::::;
C--- subroutine loops is the main sequencing control program. it calls sub-
C--- progs DFIJ @ TFIJ to obtain the error functions FIJ for the density
C--- and temperature d.e.'s. it sets up the Jacobian matrix in subroutine
C--- JMATRIX and solves for the increments using BDSLV. increments from
C--- the solver are tested to ensure non -ve densities (modified steepest
C--- descent).
      SUBROUTINE DLOOPS(JMIN,   !.. first point on the field line
     >                  JMAX,   !.. last point on the field line
     >                 FLDIM,   !.. Field line grid array dimension
     >                     Z,   !.. Altitude array
     >                     N,   !.. O+, H+, He+, minor ion densities array
     >                    TI,   !.. Ion and electron temperatures
     >                  DTIN,   !.. Time step from calling program
     >                 DTMIN,   !.. Minimum time step
     >                    SL,   !.. arc length from JMIN to grid point
     >                 EFLAG)   !.. OUTPUT: Error flag array
      USE SOLVARR       !... DELTA RHS WORK S, Variables for solver
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER FLDIM                      !.. Field line grid array dimension
      INTEGER NFLAG,EFLAG(11,11)         !.. error flags
      INTEGER I,J,JC,IPR,ITER            !.. Loop control variables
      INTEGER IDIV,KR,ION,IEQ,MIT,MAXIT  !.. solution variables
      INTEGER JBNN,JBNS,JBTN,JBTS        !.. boundary indices
      INTEGER JMIN,JEQ,JMAX,JK,JCON      !.. spatial grid indices
      INTEGER DCLON,DCLOS,DCUPP          !.. tests for convergence region
      DOUBLE PRECISION DT,DTIN,DTMIN,DTINC !.. Time step variables
      DOUBLE PRECISION ZLBDY               !.. lower boundary altitudes for N
      DOUBLE PRECISION DCRQ(2,FLDIM),DCR(2)      !.. solution variables
      DOUBLE PRECISION DCRP,DCRT,ADCR,DINC,F(20) !.. solution variables
      !.. see above for description of Z,N,TI,HE
      DOUBLE PRECISION Z(FLDIM),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION RATIN,RATIS,FACIRI,ALPHA,SL(FLDIM)
      DOUBLE PRECISION NSAVE(2,FLDIM),TISAV(3,FLDIM)  !.. N and T at previous time step

      DT=DTIN           !.. Set time step for dN/dt
      DTINC=0.0         !.. Used for reduced timestep
      JBNN=JMIN         !.. lower boundary point for O+ and H+
      JBNS=JMAX         !.. lower boundary point for O+ and H+
      JEQ=(JMIN+JMAX)/2 !.. Equatorial point
      EFLAG(2,1)=0      !.. Initialize error flag
      EFLAG(2,1)=0      !.. Initialize error flag
      ZLBDY=120.        !.. Lower boundary altitude

      !-- Save current densities for dn/dt, dT/dt. FACIRI is a factor for 
      !-- scaling FLIP to the IRI nmF2.
      DO J=JMIN,JMAX
        NSAVE(1,J)=N(1,J)
        NSAVE(2,J)=N(2,J)
        DCRQ(1,J)=1.0
        DCRQ(2,J)=1.0
      ENDDO

      !.. avden is called to obtain average midpoint densities and temps
      CALL AVDEN(TI,0,0.0D0,0.0D0)

      !.. OUTER LOOP: Return here on Non-Convergence with reduced time step
  10  CONTINUE
        !.. Update indices for the lower boundaries. 
        DO J=JMIN,JEQ-1
          IF(Z(J).LE.ZLBDY) JBNN=J      !.. North
        ENDDO
        DO J=JMAX,JEQ,-1
          IF(Z(J).LE.ZLBDY) JBNS=J      !.. South
        ENDDO

        !.. Set up boundary indices for the solution procedure
        MIT=JBNS-JBNN+1       !.. Number of points on field line
        IEQ=2*(MIT-2)         !.. Number of equations to set up      

        !.. Main loop: On each iteration the Jacobian is formed and solved for
        !.. the increments of to add to N. 
        MAXIT=10
        DO 220 ITER=1,MAXIT
          !.. boundary conditions on density. 
          DO J=JMIN,JBNN
            CALL HOEQ(FLDIM,1,J,JMAX,DT,N,TI)
          ENDDO
          DO J=JBNS,JMAX
            CALL HOEQ(FLDIM,1,J,JMAX,DT,N,TI)
          ENDDO

          !.. Compute FIJ values to use in calculating dN/dh
          DO J=2,MIT-1
            KR=2*(J-2)
            JC=J+JBNN-1
            CALL DFIJ(JC,0,IPR,DT,JEQ,N,TI,F,NSAVE)
            RHS(KR+1)=F(1)
            RHS(KR+2)=F(2)
          ENDDO

          !.. Create the Jacobian matrix
          CALL JMATRX(FLDIM,S,RHS,IEQ,IPR,DT,JEQ,Z,N,TI,F,JBNN,JBNS,MIT,
     >       NSAVE,TISAV)

          !.. Solve the linear system with the band solver
          !.. invert the jacobian matrix 'S' in the inversion routine BDSLV.
          !.. the increments are stored in array delta in this order
          !.. x(1...n,j),x(1...n,j+1),x(1...n,j+2),....x(1...n,jmax-1)
          CALL BDSLV(IEQ,3,S,0,RHS,DELTA,WORK,NFLAG)

          IF(NFLAG.NE.0) THEN
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(/A,I5,A,I5)')
     >       '  ITERATION =',ITER,'  RETURN FROM BDSLV =',NFLAG
            EFLAG(2,2)=-1   !.. Report problem to calling routine
            RETURN
          ENDIF

          IDIV=0
          DCR(1)=1
          DCR(2)=1
          DCRP=1.0

          !.. test n increments 'dinc' to ensure n>0 (mod steepest descent)
          DO  142 I=1,IPR
          DO 142 J=2,MIT-1
            DCRP=1.0
            JC=JBNN+J-1
            DCRQ(I,JC)=1.0
            ION=3
            IF(I.EQ.IPR) ION=2
            DINC=DELTA(2*J-ION)
            IF(DINC.LE.0) GO TO 142
            IF(ABS(DINC/N(I,JC)).GT.0.8) DCRP=0.5*ABS(N(I,JC)/DINC)
            IF((DCRP.LT.DCR(I)).AND.(Z(JC).GT.0.0)) DCR(I)=DCRP
            IF(ITER.GT.0) DCRQ(I,JC)=DCRP
 142      CONTINUE

          !. add iterative increment to the array 'N' and test for
          !. convergence when idiv=0. separate dcrs are switched off
          ADCR=DCR(1)
          IF(DCR(2).LT.DCR(1)) ADCR=DCR(2)
          DO I=1,IPR
            ION=3
            IF(I.EQ.IPR) ION=2
            DO J=2,MIT-1
              !..... first temperatures ....
              JC=JBNN+J-1
              DINC=DELTA(2*J-ION)
              N(I,JC)=N(I,JC)-DINC*ADCR
              IF(ABS(DINC/N(I,JC)).GT.1E-3)  IDIV=IDIV+1
            ENDDO
          ENDDO

          !. test to see if convergence has occured.
          IF(IDIV.EQ.0) GOTO 224
 220      CONTINUE

        !.. END OF SOLUTION LOOP     
 224    CONTINUE

        DCLON=1  !.. North Flag for non-convergence below 200 km
        DCLOS=1  !.. South Flag for non-convergence below 200 km
        DCUPP=1  !.. Flag for non-convergence above 200 km

        !.. Testing for convergence to see if need to reduce time step
        IF(IDIV.EQ.0) THEN
          !============== Convergence success ========
          EFLAG(2,1)=0     
          DTINC=DTINC+DT   !.. Used for reduced timestep
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,I5,9F14.2)') 
     >    '  O+,H+ ',ITER,DTINC,DTIN,DT,Z(JBNN)
          !.. Final convergence, store values and return
          IF(DTINC.GE.DTIN-1) THEN
             DO J=JMIN,JMAX
               XIONN(1,J)=N(1,J)
               XIONN(2,J)=N(2,J)
            ENDDO
            RETURN
          ENDIF
          !.. increase time step if convergence is easy
          IF(ITER.LT.5.AND.DTINC+2*DT.LE.DTIN) DT=2*DT
          IF(ITER.LT.5.AND.DTINC+2*DT.GT.DTIN) DT=DTIN-DTINC

          !-- Save current densities for dN/dt. 
          DO J=JMIN,JMAX
            NSAVE(1,J)=N(1,J)
            NSAVE(2,J)=N(2,J)
          ENDDO
        ELSE
          !============== Convergence failure ========
          !.. Test to see if difficulty only below 200 km
          DO J=JBNN,JBNS
            IF(Z(J).LT.200.0) THEN
              IF(J.LE.JEQ.AND.DCRQ(1,J).LT..99999999) DCLON=0
              IF(J.GT.JEQ.AND.DCRQ(1,J).LT..99999999) DCLOS=0
            ELSE
              IF(DCRQ(1,J).LT..99999999) DCUPP=0
            ENDIF
          ENDDO
          !-- Non-Convergence: Reduce time step and restore densities. 
          DT=DT/2  !.. reduce time step for non-convergence
          DO J=JMIN,JMAX
            N(1,J)=NSAVE(1,J)
            N(2,J)=NSAVE(2,J)
          ENDDO
          !.. Raise lower boundary if the problem is only below 200 km
          IF(DT.LE.DTIN/4.0.AND.DCUPP.EQ.1.AND.DCLON*DCLOS.EQ.0)
     >       ZLBDY=(ZLBDY+200)/2   

          !.. Check that DT is not too small
          IF(DT.LT.DTMIN) THEN
            EFLAG(2,1)=-1   !.. Report problem to calling routine
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,9I5)') 
     >        '  ERR FLAGS LPSD',DCLON,DCLOS,DCUPP,NINT(ZLBDY)
            RETURN
          ENDIF 
        ENDIF

      GOTO 10   !... END OF OUTER LOOP .........................

      RETURN
      END
C:::::::::::::::::::::::::: HOEQ :::::::::::::::::::::::::::::::::::::::::::::::::
C.....  finds the new chemical equilibrium densities
C.....  of o+ and h+ at point j for fraction implicitness
      SUBROUTINE HOEQ(FLDIM,ISW,J,JMAX,DT,N,TI)
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production, PHION
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      INTEGER FLDIM   !.. field line grid dimension
      DIMENSION RTS(99),N(4,FLDIM),TI(3,FLDIM)

      !.. Get reaction rates
      CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
      !....... evaluate loss rates
      LOPLS=(RTS(3)*N2N(J)+RTS(4)*O2N(J))
      LOPLS2=+RTS(2)*HN(J)
      LHPLS=RTS(1)*ON(J)
      N(1,J)=PHION(J)/LOPLS
      N(2,J)=N(1,J)*LOPLS2/LHPLS
       RETURN
       END
C::::::::::::::::::::::: JMATRX ::::::::::::::::::::::::::::::::::::
C......  JMATRX EVALUATES DELF/DELT USING STEPHANSONS METHOD.
C......  S(L,M) = DEL(FIJ) / DEL(N)
C......  WHERE N IS THE DENSITY OF ION IV AT ALTITUDE JV. L AND M
C......  ARE COMPUTED FROM JF...IV. THE ARRAY RHS CONTAINS VALUES OF
C...... FIJ SAVED FROM PREVIOUS COMPUTATION

      SUBROUTINE JMATRX(FLDIM,S,RHS,INEQ,IPR,DT,JEQ,Z,N,TI,F,JBNN,JBNS,
     >   MIT,NSAVE,TISAV)
      IMPLICIT DOUBLE PRECISION(A-H,N-Z)
      INTEGER FLDIM   !.. field line grid dimension
      DIMENSION RHS(INEQ),S(INEQ,8),N(4,FLDIM),TI(3,FLDIM),F(20)
     >  ,NSAVE(2,FLDIM),TISAV(3,FLDIM)

      DO 180 KZS=1,7
      DO 180 JZS=1,INEQ
 180  S(JZS,KZS)=0.0

      DO 80 JF=2,MIT-1
      J1=MAX0(2,JF-1)
      J2=MIN0(JF+1,MIT-1)
      DO 80 IV=1,2
      DO 80 JV=J1,J2

      L=2*(JF-2)
      IF(JF.LE.3)M=2*(JV-2)+IV
      IF(JF.GT.3)M=2*(JV-JF)+IV+4
      KRV=2*(JV-2)+IV

      JVC=JBNN+JV-1
      JFC=JBNN+JF-1

      !. H FOR TI CALCULATIONS ..........
      IF(IPR.EQ.2) GO TO 68
      H=1.E-4*TI(IV+1,JVC)
      TI(IV+1,JVC)=TI(IV+1,JVC)+H
      CALL TFIJ(JFC,1,IPR,DT,JEQ,N,TI,F,TISAV)
      TI(IV+1,JVC)=TI(IV+1,JVC)-H
       GO TO 70
 68    CONTINUE
      !..... DENSITIES ........
      H=1.E-4*N(IV,JVC)
      N(IV,JVC)=N(IV,JVC)+H
      CALL DFIJ(JFC,1,IPR,DT,JEQ,N,TI,F,NSAVE)
      N(IV,JVC)=N(IV,JVC)-H
 70     CONTINUE
      RH=1/H
      IF(JF.LE.3)  S(L+1,M)=(F(1)-RHS(L+1))*RH
      IF(JF.LE.3)  S(L+2,M)=(F(2)-RHS(L+2))*RH
      IF(JF.GT.3)  S(L+1,M-1)=(F(1)-RHS(L+1))*RH
      IF(JF.GT.3)  S(L+2,M-2)=(F(2)-RHS(L+2))*RH

   80 CONTINUE
      RETURN
      END
C:::::::::::::::::::::::::::: SMOOTH :::::::::::::::::::::::
C....... Smoothing out O+ profile when boundary above 120 km by 
C....... fitting an exponential from absolute lower boundary to the O+ 
C....... lower boundary. This profile is then blended with the
C....... local equilibrium values from subroutine HOEQ. The second
C....... call to HOEQ is to get the H+ density consistent with O+
      SUBROUTINE SMOOTH(FLDIM,JMAX,JBNN,JBS,Z,N,DT,TI)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      INTEGER FLDIM   !.. field line grid dimension
      DIMENSION N(4,FLDIM),TI(3,FLDIM),NSAVE(4,FLDIM),Z(FLDIM)
        JINC=1
        JBS=JMAX-JBNN+1
        JBSS=JMAX-JBS+1
        AL1=DLOG(N(1,JBNN+JINC)/N(1,JBS))/(Z(JBNN+JINC)-Z(JBS))
        AL2=DLOG(N(1,JBS-JINC)/N(1,JBSS))/(Z(JBS-JINC)-Z(JBSS))
      !.......... Northern hemisphere
        DO 511 J=JBS,JBNN+JINC
           N(1,J)=N(1,JBS)*EXP(AL1*(Z(J)-Z(JBS)))
           CALL HOEQ(FLDIM,1,J,JMAX,DT,NSAVE,TI)
           N(1,J)=(N(1,J)+NSAVE(1,J))*0.5
           !.. second call to get H+ with new O+ 
           CALL HOEQ(FLDIM,1,J,JMAX,DT,NSAVE,TI)
           N(2,J)=NSAVE(2,J)
 511    CONTINUE   

      !.......... Southern hemisphere
        DO 515 J=JBS-JINC,JBSS
           N(1,J)=N(1,JBSS)*EXP(AL2*(Z(J)-Z(JBSS)))
           CALL HOEQ(FLDIM,1,J,JMAX,DT,NSAVE,TI)
           N(1,J)=(N(1,J)+NSAVE(1,J))*0.5
           !.. second call to get H+ with new O+ 
           CALL HOEQ(FLDIM,1,J,JMAX,DT,NSAVE,TI)
           N(2,J)=NSAVE(2,J)
 515    CONTINUE   
      RETURN
      END

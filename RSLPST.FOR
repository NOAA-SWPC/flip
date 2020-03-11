C.................... RSLPST.FOR ........................................... 
      SUBROUTINE CTIPFTE(FLDIM,JMIN,JMAX,Z,N,TI,DTIN,DTMIN,EFLAG)
      IMPLICIT NONE
      INTEGER JMIN,JMAX,EFLAG(11,11),FLDIM
      DOUBLE PRECISION Z(FLDIM),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION DTIN,DTMIN  !.. Initial and minimum time steps
      CALL TLOOPS(JMIN,JMAX,FLDIM,Z,N,TI,DTIN,DTMIN,EFLAG)
      RETURN
      END
C:::::::::::::::::::::::::::::: TLOOPS :::::::::::::::::::::::::::::::::;
C--- subroutine loops is the main sequencing control program for the temperatures.
C--- It calls subroutine TFIJ to obtain the error functions FIJ for the density
C--- and temperature d.e.'s. It sets up the Jacobian matrix in subroutine
C--- JMATRIX and solves for the increments using BDSLV. increments from
C--- the solver are tested to ensure non -ve densities (modified steepest
C--- descent).
C---- This version carved out of the old RSLPSD by P. Richards May 2010
      SUBROUTINE TLOOPS(JMIN,   !.. first point on the field line
     >                  JMAX,   !.. last point on the field line
     >                 FLDIM,   !.. Field line grid array dimension
     >                     Z,   !.. Altitude array
     >                     N,   !.. O+, H+, He+, minor ion densities array
     >                    TI,   !.. Ion and electron temperatures
     >                  DTIN,   !.. Time step from calling program
     >                 DTMIN,   !.. Minimum time step
     >                 EFLAG)   !.. OUTPUT: Error flag array
      USE SOLVARR       !... DELTA RHS WORK S, Variables for solver
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      IMPLICIT NONE
      INTEGER FLDIM                      !.. Field line grid array dimension
      INTEGER NFLAG,EFLAG(11,11)         !.. error flags
      INTEGER I,J,JC,IPR,ITER            !.. Loop control variables
      INTEGER IDIV,KR,ION,IEQ,MIT,MAXIT  !.. solution variables
      INTEGER JBNN,JBNS,JBTN,JBTS        !.. boundary indices
      INTEGER JMIN,JEQ,JMAX              !.. spatial grid indices
      DOUBLE PRECISION DT,DTIN,DTMIN,DTINC !.. Time step variables
      DOUBLE PRECISION ZLBDY,ZTBDY       !.. lower boundary altitudes for N & T
      DOUBLE PRECISION DCRP,DCRT,DINC,F(20) !.. solution variables
      !.. see above for description of Z,N,TI,HE
      DOUBLE PRECISION Z(FLDIM),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION NSAVE(2,FLDIM),TISAV(3,FLDIM)
      DATA ZLBDY/120.0/,ZTBDY/100.0/

      DT=DTIN           !.. Set time step for dT/dt
      DTINC=0.0         !.. Used for reduced timestep
      JBTN=1            !.. lower boundary point for Te and Ti
      JBTS=JMAX         !.. lower boundary point for Te and Ti
      JBNN=1            !.. lower boundary point for O+ and H+
      JBNS=JMAX         !.. lower boundary point for O+ and H+
      JEQ=(JMIN+JMAX)/2 !.. Equatorial point
      EFLAG(1,1)=0      !.. Initialize error flag
      EFLAG(1,2)=0      !.. Initialize error flag

      !.. Update indices for the lower boundaries. 
      DO J=JMIN,JEQ-1
        IF(Z(J).LE.ZLBDY) JBNN=J
        IF(Z(J).LE.ZTBDY) JBTN=J
      ENDDO
      DO J=JMAX,JEQ,-1
        IF(Z(J).LE.ZLBDY) JBNS=J
        IF(Z(J).LE.ZTBDY) JBTS=J
      ENDDO
     
      !-- Save current temperatures for dT/dt. 
      DO J=JMIN,JMAX
         TISAV(1,J)=TI(1,J)
         TISAV(2,J)=TI(2,J)
         TISAV(3,J)=TI(3,J)
      ENDDO

      !.. Set up boundary indices for the solution procedure
      MIT=JBTS-JBTN+1       !.. Number of points on field line
      IEQ=2*(MIT-2)         !.. Number of equations to set up      

      !.. OUTER LOOP: Return here on Non-Convergence with reduced time step
  10  CONTINUE
        !.. Main loop: On each iteration the Jacobian is formed and solved for
        !.. the increments of to add to TI. 
        MAXIT=17
        DO 220 ITER=1,MAXIT
          !... boundary conditions on temperature in North
          DO J=JMIN,JBTN
             TI(1,J)=TN(J)
             TI(2,J)=TN(J)
             TI(3,J)=TN(J)
          ENDDO
          !... boundary conditions on temperature in South
          DO J=JBTS,JMAX
             TI(1,J)=TN(J)
             TI(2,J)=TN(J)
             TI(3,J)=TN(J)
          ENDDO
  
          !.. boundary conditions on density. Needed here for temperatures
          DO J=JMIN,JBNN
            CALL HOEQ(FLDIM,1,J,JMAX,DT,N,TI)
          ENDDO
          DO J=JBNS,JMAX
            CALL HOEQ(FLDIM,1,J,JMAX,DT,N,TI)
         ENDDO

          !.. Compute FIJ values to use in calculating dN/dh
          DO J=2,MIT-1
            KR=2*(J-2)
            JC=J+JBTN-1
            CALL TFIJ(JC,0,IPR,DT,JEQ,N,TI,F,TISAV)
            RHS(KR+1)=F(1)
            RHS(KR+2)=F(2)
          ENDDO

          !.. Create the Jacobian matrix
          CALL JMATRX(FLDIM,S,RHS,IEQ,IPR,DT,JEQ,Z,N,TI,F,JBTN,JBTS,MIT,
     >     NSAVE,TISAV)

          !.. Solve the linear system with the band solver
          !.. invert the jacobian matrix 'S' in the inversion routine BDSLV.
          !.. the increments are stored in array delta in this order
          !.. x(1...n,j),x(1...n,j+1),x(1...n,j+2),....x(1...n,jmax-1)
          CALL BDSLV(IEQ,3,S,0,RHS,DELTA,WORK,NFLAG)

          IF(NFLAG.NE.0) THEN
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(/A,I5,A,I5)')
     >        '  ITERATION =',ITER,'  RETURN FROM BDSLV =',NFLAG
            EFLAG(1,2)=-1   !.. Report problem to calling routine
            RETURN
          ENDIF

          IDIV=0
          DCRP=1.0
          DCRT=1.0

          !.. test TI increments 'DINC' to ensure TI>0 (mod steep descent)
          DO 137 I=1,IPR
          DO 137 J=2,MIT-1
            DCRP=1.0
            ION=3
            IF(I.EQ.IPR) ION=2
            DINC=DELTA(2*J-ION)
            !.. if DINC exceeds fraction of TI set the factor DCRT so that TI>0
            JC=JBTN+J-1
            IF(ABS(DINC/TI(I,JC)).GT.0.8) DCRP=0.5*ABS(TI(I,JC)/DINC)
            IF(DCRP.LT.DCRT) DCRT=DCRP
 137      CONTINUE

          !. add iterative increment to the array 'TI' and test for
          !. convergence when idiv=0. separate dcrs are switched off
          DO I=1,IPR
            ION=3
            IF(I.EQ.IPR) ION=2
            DO J=2,MIT-1
              JC=JBTN+J-1
              DINC=DELTA(2*J-ION)
              TI(I,JC)=TI(I,JC)-DINC*DCRT
              IF(ABS(DINC/TI(I,JC)).GT.1E-3)   IDIV=IDIV+1
            ENDDO
          ENDDO
          IF(IDIV.EQ.0) GO TO 224   !. test for convergence
 220    CONTINUE

        !.. END OF SOLUTION LOOP     
 224    CONTINUE

        !.. Testing for convergence to see if need to reduce time step
        IF(IDIV.EQ.0) THEN
          !============== Convergence success ========
          EFLAG(1,1)=0     
          DTINC=DTINC+DT   !.. Used for reduced timestep
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,I5,9F14.2)') 
     >       '  Te Ti ',ITER,DTINC,DTIN,DT,Z(JBTN)
          IF(DTINC.GE.DTIN-1) RETURN
          !.. increase time step if convergence is easy
          IF(ITER.LT.5.AND.DTINC+2*DT.LE.DTIN) DT=2*DT
          IF(ITER.LT.5.AND.DTINC+2*DT.GT.DTIN) DT=DTIN-DTINC

          !-- Save current temperatures for dt/dt. 
          DO J=JMIN,JMAX
            TISAV(1,J)=TI(1,J)
            TISAV(2,J)=TI(2,J)
            TISAV(3,J)=TI(3,J)
          ENDDO
        ELSE
        !============== Convergence failure ========
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,I5,9F14.2)') 
     >     ' Temp ',ITER,DTINC,DTIN,DT
          !-- Reduce time step and restore TI. 
          DT=DT/2  !.. reduce time step for non-convergence
          DO J=JMIN,JMAX
            TI(1,J)=TISAV(1,J)
            TI(2,J)=TISAV(2,J)
            TI(3,J)=TISAV(3,J)
          ENDDO
          !.. Check that DT is not too small
          IF(DT.LT.DTMIN) THEN
            EFLAG(1,1)=-1   !.. Report problem to calling routine
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,9I5)') '  ERROR IN LPST'
            RETURN
          ENDIF 
        ENDIF

      GOTO 10   !... END OF OUTER LOOP ....................

      RETURN
      END

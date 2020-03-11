C.................... RSDENA_EXB.FOR ..................................
C----------------------------------------------------------------
C   this program was written by eugene young(1978) but has been substantially
C   modified since then by p. richards. it is responsible for setting up the
C   continuity equations for minimization. it calls subr vel to get
C   velocities, chemo for chemical sources and sinks, and dave for
C   interpolated dn/dt. it integrates the continuity equations for the fluxes.
      SUBROUTINE DFIJ(J,JSJ,IPR,DT,JEQ,N,TI,F,NSAVE)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P) ISPEC
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      INTEGER ION
      DIMENSION NSAVE(2,FLDIM)
      DIMENSION VL(2),VU(2),FLU(2),FLL(2),ANM(2),PM(2),Q(2),L(2)
     > ,TINCR(2),N(4,FLDIM),TI(3,FLDIM),F(20)
      DATA FLU(1)/0.0/,FLU(2)/0.0/

      !.. ........ sts switches to steady state solution ........
      STS=1.0
      IPR=2
      ION=2  !.. number of ions to solve

      !.. . vel is called to obtain the velocities(fluxes) at the lower 1/2 pt
      FLL(1)=FLU(1)
      FLL(2)=FLU(2)
      IF(JSJ.NE.4) CALL VEL(J-1,ION,VL,FLL,N,TI,JSJ)

      !.. ..  call chemo to evaluate the source+sink term and call dave for dn/dt
      !.. ..  tincr=dn/dt, q=source, l=sink, pm(am)=future(past) cpt of tincr
      CALL CHEMO(JSJ,ION,J,Q,N,NSAVE,TI,L)
      CALL DAVE(JEQ,ION,J,ANM,PM,N,NSAVE)

      !.. .  calculate velocity at upper half point
  360 CALL VEL(J,ION,VU,FLU,N,TI,JSJ)
      BU=(BM(J)+BM(J+1))*0.5
      BL=(BM(J)+BM(J-1))*0.5

      !.. . this section sets up the continuity eqns f(i); tincr=dn/dt ; fgr=
      !.. . flux cpt ; q(i)=source, l(i)=sink.
      DO 380 I=1,ION
      TINCR(I)=(PM(I)-ANM(I))*STS/DT
      FGR=(FLU(I)/BU-FLL(I)/BL)

      F(I)=Q(I)-L(I)-TINCR(I)-FGR
      !.. ......... average velocity ........
      IF(JSJ.EQ.0) XIONV(I,J)=VU(I)           !..0.5*(VL(I)+VU(I))
      VL(I)=VU(I)

      !.. printing of factors in the continuity equations. Since 
      !.. all quantities have been divided by bm and integrated,
      !.. multiply by bd to get actual magnitudes
      IF(JSJ.EQ.4) THEN
        IF(IUN.NE.30) THEN
          IUN=6
          WRITE(IUN,665) 
        ENDIF
        BD=0.0     !..2.0*BM(J)/(SL(J+1)-SL(J-1))
        WRITE(IUN,666) J,Z(J),BD,FLU(I),FLL(I),FGR,Q(I),L(I),F(I)
     >    ,TINCR(I)
        IF(I.EQ.2) WRITE(IUN,*) '   '
      ENDIF
 380  CONTINUE
 665  FORMAT(6X,'J    ALT       BD         FLU         FLL       FGR'
     > ,9X,'Q         L          F          VEL       TINCR')
  666   FORMAT(1X,'FIJ',I4,F10.2,1P,22E11.3)
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VEL(J,ION,V,FLUX,N,TI,JSJ)
C-------------------------------------------------------------------
C   this program was written by eugene young (1978). it calculates the
C   velocities at the point z(j+1/2) from the ion momentum eqtns for
C   o+ and h+. parameters from the grid point z(j) are interpolated to
C   get their values  at the half point. the parameters from st maurice
C   and schunk pss 1975 are transferred from subr jp
C-------------------------------------------------------------------
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX RBM UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values  
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P) ISPEC
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DOUBLE PRECISION MASS(4),A(4)
      DIMENSION PP(4),QSIGN(2),V(2),FLUX(2),D(2),GAMMA(2),B(2),Q1(2)
      DIMENSION N(4,FLDIM),TI(3,FLDIM),F(20),QION(2)
      DATA QSIGN/-1.,1./
      DATA MASS/26.7616E-24,1.6726E-24,6.6904E-24,0.0/,A/16,1,4,0/

      PP(1)=DSQRT(N(1,J+1)*N(1,J))      !.. Midpoint O+ density
      PP(2)=DSQRT(N(2,J+1)*N(2,J))      !.. Midpoint H+ density

      NEU=N(1,J+1)+N(2,J+1)+N(3,J+1)+XIONN(3,J+1)  !.. upper e density
      NEL=N(1,J)+N(2,J)+N(3,J)+XIONN(3,J)          !.. lower e density
      !.. Midpoint electron densities
      NE=DSQRT(N(1,J+1)*N(1,J))+DSQRT(N(2,J+1)*N(2,J))+
     >   DSQRT(N(3,J+1)*N(3,J))+DSQRT(XIONN(3,J+1)*XIONN(3,J))

      !.. ... call jp to obtain d(i),gamma,al12,als12 @ b(i) used to calc v(i)
      CALL JP(J,JMAX,PP,D,GAMMA,AL12,ALS12,B,JSJ,MASS,A)
      !---------calculate altitude derivatives----------
      GRADNE=TEJ(J)*(NEU-NEL)/NE
      ALPP=(AL12-ALS12)*GRADTI(J)/(PP(1)+PP(2))

      !,, terms in the momentum eqn. consult st m @ s(1975) p910 eqn (19)

      DO 301 I=1,ION
      DLNI=(N(I,J+1)-N(I,J))/DS(J)/PP(I)
      !.. .      DLNI=DLOG(N(I,J+1)/N(I,J))/DS(J)
      GRA=-MASS(I)*GRAV(J)
      IK=3-I
      QION(I)=QSIGN(I)*(GAMMA(I)*GRADTE(J)-PP(IK)*ALPP)

      Q1(I)=DLNI+GRA+GRADTI(J)+GRADNE+GRADTE(J)+QION(I)

      !.. .... neutral wind term .....
      Q2=UNJ(J)*B(I)
      !.. .... form the rhs of the momtm eqn.
  729  F(I)=-D(I)*Q1(I)+Q2

      !....  print routine for momtm eqn terms ..........
      IF(JSJ.EQ.4) THEN
        IF(IUN.NE.31) THEN
          IUN=31
          WRITE(IUN,665) 
        ENDIF
        WRITE(IUN,605) J,NINT(Z(J)),DLNI,GRADTI(J),GRADNE
     >     ,GRADTE(J),GRA,QION(I),Q1(I),F(I),Q2,NE
     >     ,NINT(TI(3,J)),NINT(TI(2,J))
        IF(I.EQ.2) WRITE(IUN,*) '   '
      ENDIF
 301   CONTINUE

 665  FORMAT(6X,'J    ALT    DLNI      GRADTI    GRADNE'
     > ,3X,'GRADTE     GRAV      QION       Q1         F'
     > ,8X,'Q2        NE      Te    Ti')
 605   FORMAT(1X,'VEL',I4,I7,1P,10E10.2,2I6)

      !--------------now the velocities--------------
      DET=B(1)+B(2)+B(1)*B(2)
      V(1)=((F(1)+F(2))+B(2)*F(1))/DET
      V(2)=((F(1)+F(2))+B(1)*F(2))/DET
      FLUX(1)=V(1)*PP(1)
      FLUX(2)=V(2)*PP(2)
      !.. ....print velocities ...........
      IF(JSJ.EQ.-4) THEN
        IUN=32
        !...IF(J.GT.JMAX/2) IUN=9
        WRITE(IUN,606) J,Z(J),V(1),V(2),DET,F(1),F(2)
      ENDIF
 302  CONTINUE
 606  FORMAT(1X,'J=',I4,3X,'ALT=',F9.0,2X,'V1=',1P,E13.6,2X,'V2=',E13.6
     >  ,2X,'DET=',E13.6,2X,'F1=',E13.6,2X,'F2=',E13.6)
      RETURN
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C-------------------------------------------------------------------
C--- this routine provides the parameters on p910 of st maurice
C--- @ schunk pss 1975.  equations 20,21,22,23,24,25,26 and the
C--- ratio of neutral/ion coll freqs-b.
C-------------------------------------------------------------------
      SUBROUTINE JP(JI,JMAX,PP,D,GAMMA,AL12,ALS12,B,JSJ,MASS,A)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX RBM UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values  
      IMPLICIT  DOUBLE PRECISION(A-H,N-Z)
      DOUBLE PRECISION M12,MASS(4),A(4)
      DIMENSION AMR(2,2),NU(2,4)
     >   ,D1(2,2),D4(2,2),NUP(2,2),D(2),GAMMA(2),B(2),PP(4)
      DATA BK/1.3807E-16/,AMU/1.6726E-24/,M12/16.0/
      DATA AMR/8,0.94,0.94,0.5/, D1/0.725,-0.16,2.66,0.725/
     >  ,D4/-0.075,0.98,-0.079,-0.075/

      !.. ....... ion-ion collisions from st.m.@s p920 ............

      DO 113 I=1,2
      DO 113 J=1,2
      NU(I,J)=1.2726*PP(J)*DSQRT(AMR(I,J))/(A(I)*TIJ(JI)**1.5)
  113 CONTINUE

      !.. . nup=nu' (from st m @ s p920) ...........
      DO 20 I=1,2
      K=3-I
      DO 20 J=1,2
      IF(I.NE.J) NUP(I,J)=1.25*NU(I,J)*(D4(I,J)+1.5*AMR(I,J)/A(I))
      IF(I.EQ.J) NUP(I,J)=NU(I,J)+
     > 1.25*NU(I,K)*(D1(I,K)+1.5*AMR(I,K)/A(I))
20    CONTINUE

      N12=PP(1)/PP(2)
      X=NU(1,2)/NUP(1,1)*(1.-NUP(2,1)/NUP(2,2))
      X=X/(1.-NUP(1,2)*NUP(2,1)/(NUP(1,1)*NUP(2,2)))
      AL12=1.875*X*(N12+1.)/(M12+1.)
      ALS12=AL12*M12*M12*(NUP(1,1)-NUP(1,2))/(NUP(2,2)-NUP(2,1))
      DELTA=3./(5.*(N12+1.)*(M12+1.))*(AL12+N12*M12*ALS12)
      RDELTA=1.0/(1.0-DELTA)


      DO 40 I=1,2
      GAMMA(I)=0.0
      J=3-I
      !--- diffusion coeff d(i) on page 910 and ratio of coll freqs b(i) -------
      D(I)=BK*TIJ(JI)*RDELTA/MASS(I)/NU(I,J)
      B(I)=NUX(I,JI)*RDELTA/NU(I,J)
 40   CONTINUE
     
      !..... output diagnostics
      IF(JSJ.EQ.4) THEN
        IF(IUN.NE.33) THEN
          IUN=33
          WRITE(IUN,665) 
        ENDIF
        WRITE(IUN,99) JI,PP(1),PP(2),NU(1,2),NU(2,1),D(1)
     >  ,D(2),B(1),B(2),DELTA,ALS12,AL12
        IF(I.EQ.2) WRITE(IUN,*) '   '
      ENDIF
 665  FORMAT(6X,'J   PP(1)      PP(2)   NU(1,2)    NU(2,1)    D(1)'
     > ,6X,'D(2)     B(1)     B(2)      DELTA     ALS12    AL12')
 99   FORMAT(1X,'JP',I4,1P,22E10.3)
      RETURN
      END
C::::::::::::::::::::::::::::: AVDEN :::::::::::::::::::::::::::::::::::
C.. ..................................................................
C   this program evaluates the interpolated densities at the midpoints
C   it also evaluates the ion-neutral collision frequencies nux(i,j)
C   and o+ and h+ production and loss rates
C.. .....................................................................
      SUBROUTINE AVDEN(TI,JSJ,ZLO,ZHI)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !.. TEJ TIJ NUX RBM UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values  
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION TI(3,FLDIM),RTS(99)
      DATA BK/1.3807E-16/
      DO 10 J=JMIN,JMAX-1
      RBM(J)=1.0/BM(J)
      DS(J)=SL(J+1)-SL(J)
      TIJ(J)=0.5*(TI(1,J)+TI(1,J+1))
      TEJ(J)=0.5*(TI(3,J)+TI(3,J+1))/DS(J)/TIJ(J)
      GRAV(J)=0.5*(GR(J)+GR(J+1))/BK/TIJ(J)
      GRADTI(J)=(TI(1,J+1)-TI(1,J))/DS(J)/TIJ(J)
      GRADTE(J)=(TI(3,J+1)-TI(3,J))/DS(J)/TIJ(J)
      UNJ(J)=0.5*(UN(J+1)+UN(J))

      !..  ion-neutral coll freqs taken from schunk and nagy rev. geophys. v18
      !..  p813, 1980. first,  resonant charge exchange from table 5 p823
      OA=SQRT(ON(J)*ON(J+1))
      HA=SQRT(HN(J)*HN(J+1))
      N2A=SQRT(N2N(J)*N2N(J+1))
      O2A=SQRT(O2N(J)*O2N(J+1))
      TNJ=0.5*(TN(J)+TN(J+1))
      TR=(TIJ(J)+TNJ)*0.5
      SQT=SQRT(TR)
      RHPH=2.65E-10*HA*SQT*(1-0.083*DLOG10(TR))**2
      !.. .... COLFAC is the scaling factor for O+ - O collision frequency  
      !.. .... should be 1.7 according to Burnside 1987 (in press)
      ROPO= COLFAC * 3.67E-11*OA*SQT*(1-0.064*DLOG10(TR))**2 
      RHPO=6.61E-11*OA*SQRT(TIJ(J))*(1-0.047*DLOG10(TIJ(J)))**2
      ROPH=4.63E-12*HA*SQRT(TNJ+TIJ(J)/16.0)
      !.. ..... non-resonant ion neutral interactions table 6
      CHPN=(33.6*N2A+32*O2A)*1.0E-10
      COPN=(6.28*N2A+6.64*O2A)*1.0E-10
      NUX(1,J)=ROPO+ROPH+COPN
      NUX(2,J)=RHPH+RHPO+CHPN

      !.. .... o+ and h+ production and loss rates
      CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
      OLOSS(J)=RTS(2)*HN(J)+RTS(3)*N2N(J)+RTS(4)*O2N(J)
      HLOSS(J)=RTS(1)*ON(J)
      HPROD(J)=RTS(2)*HN(J)
      IF(JSJ.NE.4.OR.Z(J).LT.ZLO.OR.Z(J).GT.ZHI) GO TO 10
        IUN=34
        !... IF(J.GT.JMAX/2) IUN=9
        WRITE(IUN,99) J,Z(J),OA,HA,N2A,O2A,TIJ(J),TNJ,TR,RHPH,ROPO,RHPO
        WRITE(IUN,99) J,Z(J),ROPH,CHPN,COPN,NUX(1,J),NUX(2,J),OLOSS(J)
     >   ,HLOSS(J),HPROD(J)
 10   CONTINUE
 99   FORMAT(1X,'AVDEN',I4,F7.0,1P,22E10.3)
      RETURN
      END
C:::::::::::::::::::::::::::: TERLIN :::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE TERLIN(J,QJM1,QJ,QJP1,SL,SU,QL,QM,QU)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      !.. .. ql and qu are interpolation points
      !.. .. of q in the (zb,zj) and (zj,za)intervals respectively
      !.. .. qm is the average value of q in the (zj-yl,zj+yl) interval
      QL=.5*(QJM1+QJ)
      QU=.5*(QJ+QJP1)
      !**** The integration is to be done only from the lower midpt to the
      !.. upper midpt. lower midpt. is the midpt of points j and j-1 and
      !.. upper one is the midpt of points j and j+1. the dist between the
      !.. j and j-1 = sl and that between j and j+1 = su. the integration 
      !.. is done only over sl/2 and su/2. hence the extra division by 2.
      QM=.25*((QL+QJ)*SL+(QU+QJ)*SU)
      RETURN
      END
C::::::::::::::::::::::::::::::::: TEREXP ::::::::::::::::::::::::::::::::::::::::::
C.. .... This subroutine calculates the integral of Q represented by
C.. .... the 3 values QJM1, QJ, QJP1 from SL=(S(J)+S(J-1))/2 to
C.. .... SU=(S(J+1)+S(J))/2. The integral is done in 2 parts,
C.. .... from SL to SJ and from SJ to SU.
C.. .... For exponential interpolation, the mean value
C.. .... of A and B is SQRT(A*B)
      SUBROUTINE TEREXP(J,QJM1,QJ,QJP1,SL,SU,QL,QM,QU)
      IMPLICIT DOUBLE PRECISION (A-H,L,N-Z)
      IF((QJM1*QJ.GT.0).AND.(QJ*QJP1.GT.0))GO TO 100

 100  CONTINUE
      !.. ...... find values of Q at midpoints
      QL=DSQRT(QJM1*QJ)
      QU=DSQRT(QJ*QJP1)
      !.. .... integral from SL to SU
      QM=0.5*(SU*(QJ-QU)/DLOG(QJ/QU)+SL*(QL-QJ)/DLOG(QL/QJ))      
      RETURN
      END
C:::::::::::::::::::::::::::::: DAVE :::::::::::::::::::::::::::::::::::::::::::::
C.. . calculates ante and post values of half interval, and average
C.. . of (ion density)/(magnetic field)
      SUBROUTINE DAVE(JEQ,ION,J,ANM,PM,N,NSAVE)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX RBM UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values  
      IMPLICIT DOUBLE PRECISION(A-H,N-Z)
      DIMENSION ANM(2),PM(2),N(4,FLDIM),NSAVE(2,FLDIM)

      DO 100 I=1,ION
      !.. .     ante values
      BB=NSAVE(I,J-1)*RBM(J-1)
      C=NSAVE(I,J)*RBM(J)
      A=NSAVE(I,J+1)*RBM(J+1)
      IF(J.LT.JEQ) CALL TEREXP(J,BB,C,A,DS(J-1),DS(J),QL,QM,QU)
      IF(J.GE.JEQ) CALL TEREXP(J,A,C,BB,DS(J),DS(J-1),QL,QM,QU)
      ANM(I)=QM
      !.. ...... post values
      BB=N(I,J-1)*RBM(J-1)
      C=N(I,J)*RBM(J)
      A=N(I,J+1)*RBM(J+1) 
      IF(J.LT.JEQ) CALL TEREXP(J,BB,C,A,DS(J-1),DS(J),QL,QM,QU)
      IF(J.GE.JEQ) CALL TEREXP(J,A,C,BB,DS(J),DS(J-1),QL,QM,QU)
      PM(I)=QM
 100  CONTINUE
      RETURN
      END
C:::::::::::::::::::::::::::::::: CHEMO :::::::::::::::::::::::::::::::::::::::::::
C.. . this program determines the interpolated production and loss
C.. . processes. it calls rates to get the rate constants and terd
C.. . to do the interpolation
      SUBROUTINE CHEMO(JSJ,ION,JI,SOURCE,N,NSAVE,TI,SINK)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !.. TEJ TIJ NUX RBM UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production, PHION
      USE AVE_PARAMS    !.. midpoint values  
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION N(4,FLDIM),TI(3,FLDIM),NSAVE(2,FLDIM),Q(2,3),L(2,3)
     > ,SOURCE(2),SINK(2)

      !-- Q and L are chemical source and sink
      DO 300 K=1,3
      J=K+JI-2
      OPROD=HLOSS(J)*N(2,J)+PHION(J)
      HPRODN=HPROD(J)*N(1,J)
      Q(1,K)=OPROD*RBM(J)
      Q(2,K)=HPRODN*RBM(J)
      L(1,K)=OLOSS(J)*N(1,J)*RBM(J)
      L(2,K)=HLOSS(J)*N(2,J)*RBM(J)

 300  CONTINUE

      !.. ..... terd is called to interpolate
      DO 400 I=1,ION
      CALL TERLIN(JI,Q(I,1),Q(I,2),Q(I,3),DS(JI-1),DS(JI),QL,QM,QU)
      CALL TERLIN(JI,L(I,1),L(I,2),L(I,3),DS(JI-1),DS(JI),LL,LM,LU)
      SOURCE(I)=QM
      SINK(I)=LM
 400  CONTINUE
      IF(JSJ.EQ.4) THEN
        IUN=35
        !... IF(JI.GT.JMAX/2) IUN=9
        WRITE(IUN,99) JI,Z(JI),Q(2,1),Q(2,2),Q(2,3),SOURCE(1)
     >   ,SOURCE(2),SINK(1),SINK(2),SL(JI),DS(JI)
      ENDIF
      RETURN
 99   FORMAT(1X,'CHEMO',I4,F7.0,1P,22E10.3)
      END
C.................... RSTEMC_EXB.FOR .............................
C.... this subroutine sets up the error functions of the time-dependent
C.... ion and electron energy equations for subsequent solution by the 
C.... Newton iterative procedure. It calls a function for Ti and one for
C.... Te. This version was created by Silvy John to implement the
C.... flux preserving method (1996). The flux preserving approach
C.... is explained in the paper by Torr et al., JGU, December 1990
      SUBROUTINE TFIJ(J,ILJ,IPR,DT,JEQ,N,TI,F,TISAV)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      IMPLICIT DOUBLE PRECISION(A-H,K-L,N-Z)
      DIMENSION N(4,FLDIM),TI(3,FLDIM),F(20),TISAV(3,FLDIM)
      IPR=3             !... signifies Temp sol
      ID=3                 !... unit for diagnostics in North
      IF(J.GT.JEQ) ID=9    !... unit for diagnostics in South

      !.... set up the ion and electron energy equations
      CALL SOLVE_IONHEAT(J,ILJ,ID,DT,JEQ,N,F,TI,TISAV)
      CALL SOLVE_ELECTRON_HEAT(J,ILJ,ID,DT,JEQ,N,F,TI,TISAV)
      RETURN
      END

C::::::::::::::::::::::::::: SOLVE_ELECTRON_HEAT :::::::::::::::::::::
C.. This subroutine sets up the integrated electron energy equation for
C.. the Newton solver F(2)
C.. It calculates the electron cooling rate, the derivative of 
C.. temperature,T, wrt to time , the heating rate, heat_flux and the
C.. loss rate for three adjacent grid points j,j-1 and j+1. These are 
C.. stored and later used to determine the values at the upper and lower
C.. midpoints by interpolation. The upper midpoint corresponds to the
C.. point midway between j and j+1 whereas lower midpoint corresponds to
C.. that between j-1 and j. The other terms are Q_UPPR and Q_LOWR which
C.. are the integrated values of the heat flow between j+1 and j && j and
C.. j-1 respectively.The basic equation is 
C.. F(2)= DT/DT(DERIVT) + DIV(HEAT_FLUX) + ELEC_ION_COOLING_RATE(E_CLRT)
C..       + ELEC_NEUT_COOLING(HLOSS_E) - HEATING_RATE(HT_RAT) = 0
C.. TERLIN interpolates and the integration is done between the midpoints.
      SUBROUTINE SOLVE_ELECTRON_HEAT(JI,ILJ,ID,DT,JEQ,N,F,TI,TISAV)
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,K-L,N-Z)
      DIMENSION N(4,FLDIM),TI(3,FLDIM),TISAV(3,FLDIM),F(20),
     > DER_T(3),CLRT(3),HLOSE(3),DEL_S(2),
     > HT_ERG(3),KE(3),kNT(3)
      DATA  BOLTZ /1.3807E-16/

      !... evaluate terms at J-1, J, J+1 for integration, JI is the
      !... actual grid point. All terms divided by magnetic field (BM)
      DO 300 I=1,3
        J=I+JI-2
        NE = N(1,J)+N(2,J)+N(3,J)+XIONN(3,J)   !... electron density
        !... temperature derivative dT/dt 
        DER_T(I) = 1.5*BOLTZ*NE*(TI(3,J)-TISAV(3,J))/BM(J)/DT
        !... electron heating rate in ergs
        HT_ERG(I) = EHT(3,J) * 1.6E-12/BM(J)
        !.. determine the thermal conductivity KE 
        CALL GET_THERMCONDUCTIVITY(TI,NE,KE,I,J)
        !.. ion and neutral cooling rates
        CALL GET_LOSS(I,J,TI,N,HLOSE,CLRT)
        !.. Convection term d(vel)/ds coefficient
        kNT(I)=BOLTZ*NE*TI(3,J)
 300  CONTINUE

      !.. distances along field line between grid points for the integration 
      DEL_S(2) = SL(JI+1) - SL(JI)
      DEL_S(1) = SL(JI) - SL(JI-1)

      !... Now perform the integration. DERIVT, E_CLRT, LOSS_E, and 
      !... HT_RAT are the four values after integration 
      CALL TERLIN(JI,DER_T(1),DER_T(2),DER_T(3),DEL_S(1),DEL_S(2),
     >  DTL,DERIVT,DTU)
      CALL TERLIN(JI,CLRT(1),CLRT(2),CLRT(3),DEL_S(1),DEL_S(2),
     >  CLRT_L,E_CLRT,CLRT_U)
      CALL TERLIN(JI,HLOSE(1),HLOSE(2),HLOSE(3),DEL_S(1),DEL_S(2),
     >  HLOSE_L,LOSS_E,HLOSE_U)
      CALL TERLIN(JI,HT_ERG(1),HT_ERG(2),HT_ERG(3),DEL_S(1),DEL_S(2),
     >  EH_L,HT_RAT,EH_U)

      !.... Calculate the divergence of the heat flux term
      KU_AVG= (KE(3) + KE(2))/2.0
      KL_AVG = (KE(2) + KE(1))/2.0
      BU_AVG = (BM(JI) + BM(JI+1)) * 0.5
      BL_AVG = (BM(JI) + BM(JI-1)) * 0.5
      Q_UPPR = KU_AVG *(TI(3,JI+1)-TI(3,JI))/DEL_S(2)/BU_AVG
      Q_LOWR = KL_AVG *(TI(3,JI)-TI(3,JI-1))/DEL_S(1)/BL_AVG
      !.... Calculate the divergence of the velocity term (SMALL)
      QV_UPPR = 0.5*(kNT(3)+kNT(2))*
     >   (XIONV(1,JI+1)-XIONV(1,JI))/DEL_S(2)/BU_AVG
      QV_LOWR = 0.5*(kNT(2)+kNT(1))*
     >   (XIONV(1,JI)-XIONV(1,JI-1))/DEL_S(1)/BL_AVG

      !... now form the error function for Te
      F(2) = DERIVT + E_CLRT + LOSS_E - HT_RAT -  Q_UPPR + Q_LOWR
 
      !... This section for diagnostic print
      IF(ILJ.EQ.4) THEN
        WRITE(ID,501)
        WRITE(ID,5555) Z(JI),F(2),DERIVT,E_CLRT,LOSS_E,Q_UPPR,Q_LOWR
     >    ,HT_RAT,QV_UPPR,QV_LOWR
      ENDIF
 501  FORMAT(1X,'Te: ALT     F        dn/dt      C_ei      C_en'
     >   ,6X,'Q_UPPR   Q_LOWR    H_tot')
 5555 FORMAT(F9.1,1P,22E10.2)
      RETURN
      END


C:::::::::::::::::::::::::::::: GET_LOSS ::::::::::::::::::::::::::::::::::::::::::::
C   PURPOSE:	calculates the cooling rates, hlose for electrons
C   ARGUMENT LIST: 
C   REFERENCES: Rees and Roble, Rev. Geophys.(1975) P220 
C               Schunk and Nagy Rev Geophys (1978) p366
C   NAME		USE
C   ----		---		
C   I,J			LOOP VARIABLES
C   TI			TEMPERATURE
C   NE			ELECTRON DENSITY
C   HLOSE		loss OF ELECTRONS
C   BM			MAGNETIC FIELD  			
      SUBROUTINE GET_LOSS(I,J,TI,N,HLOSE,CLRT)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,K-L, N-Z)
      DIMENSION HLOSE(3),TI(3,FLDIM),CLRT(3),N(4,FLDIM)

      !.... electron-ion cooling rate
      NE = N(1,J)+N(2,J)+N(3,J)+XIONN(3,J)
      NI = N(1,J)/16 + N(2,J) + N(3,J)/30.0

      !.. Electron- ion cooling rate KEI. The Coulomb logarithm from Itikawa
      !..  JATP 1975, page 1001. Previously assumed to be 15.0
      DEBYE = 1.602E-12 * SQRT(4.0*3.1415926*NE/(1.381E-16*TI(3,J)))
      COULOG = 0.43086 + DLOG(TI(3,J)) - DLOG(DEBYE)
      CLRT(I) = 1.232e-17*NE*NI*(TI(3,J)-TI(2,J))/BM(J)/(TI(3,J)**1.5)
      CLRT(I)= CLRT(I)*COULOG/15.0               !.. Itikawa correction
      CLRT(I)= CLRT(I)                      !.. fudge factor

      !... Neglect losses to neutrals above 1000 km
      IF(O2N(J).LT.1.0E-10) THEN
         HLOSE(I)=0.0
         RETURN
      ENDIF

      !.. These terms are evaluated here for efficiency because
      !.. they are used in many places
      TE2=TI(3,J)
      RTE2=1.0/TE2
      TNJ=TN(J)
      RTNJ=1.0/TNJ
      SQTE=SQRT(TE2)
      TDIF=TE2-TNJ
      TDIFRT=TDIF*RTE2*RTNJ

       !..  elec-neut elas coll loss rates s&N(1978) rev geophys p366
       LEN2=1.77E-19*NE*N2N(J)*(1.0-1.2E-4*TE2)*TE2*TDIF
       LEO2=1.21E-18*NE*O2N(J)*(1.0+3.6E-2*SQTE)*SQTE*TDIF
       LEO=7.9E-19*NE*ON(J)*(1.0+5.7E-4*TE2)*SQTE*TDIF
       LEH=9.63E-16*NE*HN(J)*(1.0-1.35E-4*TE2)*SQTE*TDIF

       !...  rotational loss rates b&k 268 ####
       LRN2=2.0E-14*NE*N2N(J)*TDIF/SQTE
       LRO2=7.0E-14*NE*O2N(J)*TDIF/SQTE

       !...  n2 vib loss rates b&k p268 and s&n p364  +++++
       EF=1.06E+4+7.51E+3*TANH(1.1E-3*(TE2-1800.0))
       GE=3300.0+1.233*(TE2-1000.0)-2.056E-4*(TE2-1000.0)*
     >   (TE2-4000.0)
       LVN2=-2.99E-12*NE*N2N(J)*EXP(5.0E-4*RTE2*EF*(TE2-2000.0))
     > *(EXP(-GE*TDIFRT)-1.0)

       !... o2 vib loss rate s&n p364  $$$$
       HS=3300.0-839.0*SIN(1.91E-4*(TE2-2700.0))
       LVO2=-5.196E-13*NE*O2N(J)*EXP(1.4286E-3*RTE2*HS*(TE2-700.0))
     >        *(EXP(-2770.0*TDIFRT)-1.0)

       !... fine structure cooling by atomic oxygen
       D1=EXP(-228.0*RTNJ)
       D2=EXP(-326.0*RTNJ)
       E1=EXP(-228.0*RTE2)
       E2=EXP(-326.0*RTE2)
       E3=EXP(-98.0*RTE2)*D1
       LF1=8.49E-6*TE2**0.519*(0.02*(D1-E1)-5.91E-9*
     >     TDIF*(2.019*D1+(228.0*RTE2+2.019)*E1))
       LF2=7.7E-6*TE2**0.3998*(0.028*(D2-E2)-5.91E-9*
     >     TDIF*(1.8998*D2+(326.0*RTE2+1.8998)*E2))
       LF3=2.22E-7*TE2**0.768*(0.008*(D2-E3)-5.91E-9*
     >     TDIF*(2.268*D2+(98.0*RTE2+2.268)*E3))
       ZFO=5.0+3.0*D1+D2
       LFO=-8.629E-6*NE*ON(J)*(LF1+LF2+LF3)/ZFO
 
       !.. excitation of O(1D)  S&N P365. 
       DE=2.4E4+(TE2-1500.0)*(0.3-1.947E-5*(TE2-4000.0))
       EXH=3.3333E-4*DE*(TE2-3000.0)*RTE2
       IF(EXH.GT.70.0)EXH=70.0
       EXH2=22713.0*TDIFRT
       IF(EXH2.GT.70.0)EXH2=70.0
       LF1D=-1.57E-12*NE*ON(J)*EXP(EXH)*(EXP(-EXH2)-1.0)

 14    KEN=LEN2+LEO2+LEO+LEH+LRN2+LRO2
       HLOSE(I)=(KEN+LVN2+LVO2+LFO+LF1D)*1.6E-12/BM(J)
       RETURN
       END

C:::::::::::::::::::::::::: GET_THERMCONDUCTIVITY ::::::::::::::::::::::::::::::::::::::::::::::::
C.....calculates ke, the thermal conductivity co-efficient
      SUBROUTINE GET_THERMCONDUCTIVITY(TI,NET,KE,I,J)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      IMPLICIT DOUBLE PRECISION(A-H,K-L,N-Z)
      DIMENSION KE(3),TE(3), TI(3,FLDIM)
      TE(I)=TI(3,J)    !.. electron temperature
      SQTE=SQRT(TE(I))
      KN2=(2.82E-17*SQTE-3.41E-21*SQTE*TE(I))*N2N(J)
      KO2=(2.2E-16+7.92E-18*SQTE)*O2N(J)
      KO=3.4E-16*ON(J)
      KHE=5.6E-16*HE(J)
      KH=(5.47E-15-7.45E-19*TE(I))*HN(J)
      !.. KN is Bank's factor for neutral effects on conductivity.
      !.. 1.0E+4 is added to avoid numerical problems for low e densities
      KN=3.22E+4*TE(I)**2*(KN2+KO2+KO+KHE+KH)/(NET + 1.0E+4)
      KE(I)=1.232E-6*(TE(I)**2.5)/(1.0+KN)   !-- Thermal conductivity 
      RETURN
      END
C:::::::::::::::::::::::::: SOLVE_IONHEAT ::::::::::::::::::::::::::
C.. This subroutine is similar to SOLVE_IONHEAT. Look at the Te function
C.. for more explanation. It determines the interpolated production and
C.. loss processes for ions. It calls RATES to get the rate constants
C.. and TERLIN to do the interpolation
      SUBROUTINE SOLVE_IONHEAT(JI,ILJ,ID,DT,JEQ,N,F,TI,TISAV)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !.. TEJ TIJ NUX RBM UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values  
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,K-L,N-Z)
      DOUBLE PRECISION N(4,FLDIM),TI(3,FLDIM),TISAV(3,FLDIM),F(20)
     > ,DER_T(3),CLRT(3),HLOSI(3),KI(3), DEL_S(2)
     > ,HT_ERG(3),COOL(22),kNT(3)

      DATA BOLTZ/1.3807E-16 /
  
      !... evaluate terms at j-1, j, J+1 for integration, JI is the
      !... actual grid point. All terms divided by magnetic field (BM)
      DO 300 I=1,3
        J=I+JI-2
        TNJ = TN(J)
        NE = N(1,J) + N(2,J) + N(3,J)+XIONN(3,J)
        NI = N(1,J)/16 + N(2,J) + N(3,J)/30.0
        RBM(J) = 1/BM(J)
        DER_T(I) = 1.5*BOLTZ*NE*(TI(2,J) - TISAV(2,J))*RBM(J)/DT
        !.. direct heating rate for ions (ring current?)
        HT_ERG(I) = EHT(1,J) * 1.6E-12/BM(J)

        !.. electron-ion heating/cooling rate
        !.. Electron- ion cooling rate KEI. The Coulomb logarithm from Itikawa
        !..  JATP 1975, page 1001. Previously assumed to be 15.0
        DEBYE = 1.602E-12 * SQRT(4.0*3.1415926*NE/(1.381E-16*TI(3,J)))
        COULOG = 0.43086 + DLOG(TI(3,J)) - DLOG(DEBYE)
        CLRT(I)=1.232E-17*NE*NI*(TI(3,J)-TI(2,J))*RBM(J)/(TI(3,J)**1.5)
        CLRT(I)= CLRT(I)*COULOG/15.0               !.. Itikawa correction
        CLRT(I)= CLRT(I)                      !.. fudge factor

        !.. loss rate coeff to neutrals (kin) rees & roble(1975) P220 ,,,
        !.. RCE = resonant charge exchange, pol = polarization interaction
        RCE=(0.21*N(1,J)*ON(J)+1.4*N(2,J)*HN(J))
        POL1=N(1,J)*(6.6*N2N(J)+2.8*HE(J)+5.8*O2N(J)+5.6*ON(J))
     >       +N(2,J)*(5.5*HE(J)+1.9*HN(J))
        POL2=(0.36*N(1,J)*HN(J)+0.4*N(2,J)*ON(J))*SQRT(TNJ)
        KIN=(RCE*SQRT(TNJ+TI(2,J))+POL1+POL2)*1.6E-26
        CALL ION_NEUTRAL(0,N(1,J),N(2,J),ON(J),N2N(J),O2N(J)
     >   ,HN(J),HE(J),TN(J),TI(2,J),Z(J),COOL)

        HLOSI(I)=+KIN*(TI(2,J)-TNJ)*RBM(J)
        HLOSI(I)=+COOL(1)*(TI(2,J)-TNJ)*RBM(J)

        !.. determine the thermal conductivity KI 
        KI(I)=1.84E-8*(N(1,J)+4*N(2,J))*(TI(2,J)**2.5)/NE
        !.. Convection term d(vel)/ds coefficient
        kNT(I)=BOLTZ*NE*TI(2,J)
 300  CONTINUE

       !.. distances between grid points
       DEL_S(2) = SL(JI+1) - SL(JI)
       DEL_S(1) = SL(JI) - SL(JI-1)

      !.. Do interpolations. derivT, E_CLRT, loss , drft_trm are the
      !.. values after integration 
      CALL TERLIN(JI,DER_T(1),DER_T(2),DER_T(3),DEL_S(1),DEL_S(2),
     >	DTL,DERIV_T,DTU)
      CALL TERLIN(JI,CLRT(1),CLRT(2),CLRT(3),DEL_S(1),DEL_S(2),
     >  CLRT_L,E_CLRT,CLRT_U)
      CALL TERLIN(JI,HLOSI(1),HLOSI(2),HLOSI(3),DEL_S(1),DEL_S(2),
     >  HLOSI_L,LOSS,HLOSI_U)
      CALL TERLIN(JI,HT_ERG(1),HT_ERG(2),HT_ERG(3),DEL_S(1),DEL_S(2),
     >  EH_L,HT_RAT,EH_U)

      !.. ion thermal conductivity
      KU_AVG= (KI(3) + KI(2))/2.0
      KL_AVG = (KI(2) + KI(1))/2.0
      TU_AVG = (TI(2,JI+1)+TI(2,JI))/2.0
      TL_AVG = (TI(2,JI)+TI(2,JI-1))/2.0
      BU_AVG = (BM(JI) + BM(JI+1)) * 0.5
      BL_AVG = (BM(JI) + BM(JI-1)) * 0.5
      Q_UPPR = KU_AVG *(TI(2,JI+1)-TI(2,JI))/DEL_S(2)/BU_AVG
      Q_LOWR = KL_AVG *(TI(2,JI)-TI(2,JI-1))/DEL_S(1)/BL_AVG
      !.... Calculate the divergence of the velocity term (SMALL)
      QV_UPPR = 0.5*(kNT(3)+kNT(2))*
     >   (XIONV(1,JI+1)-XIONV(1,JI))/DEL_S(2)/BU_AVG
      QV_LOWR = 0.5*(kNT(2)+kNT(1))*
     >   (XIONV(1,JI)-XIONV(1,JI-1))/DEL_S(1)/BL_AVG

      !.. integrated energy equation for ions
      F(1) = DERIV_T-E_CLRT+LOSS-Q_UPPR+Q_LOWR-HT_RAT
      
      !.. diagnostic print
      IF(ILJ.EQ.4) THEN
        WRITE(ID,501)
        WRITE(ID,5555) Z(JI),F(1),DERIV_T,E_CLRT,LOSS,Q_UPPR,Q_LOWR
     >      ,QV_UPPR,QV_LOWR
      ENDIF
 501  FORMAT(' Ti: ALT     F        dn/dt      C_ie      C_in'
     >   ,6X,'Q_UPPR   Q_LOWR')
 5555   FORMAT(F9.1,1P,22E10.2)
      RETURN
      END
C:::::::::::::::: ION_NEUTRAL ::::::::::::::::::::::
C.... loss rate coeff to neutrals (kin) rees & roble(1975) P220 ,,,
      SUBROUTINE ION_NEUTRAL(IWR,OPLUS,HPLUS,ON,N2N,O2N,HN,HEN
     >   ,TN,TI,ALT,COOL)
      IMPLICIT NONE
      INTEGER IWR,IK
      DOUBLE PRECISION OPLUS,HPLUS,ON,N2N,O2N,HN,HEN,TN,TI,KIN
     > ,RCE,POL1,POL2,ALT2,COOL(22),ALT
        !.. Resonant charge exchange
        COOL(2)=1.6E-26*(0.21*OPLUS*ON)*SQRT(TI+TN)
        COOL(3)=1.6E-26*(1.4*HPLUS*HN)*SQRT(TI+TN)
        !... Polarization interactions
        COOL(4)=1.6E-26*OPLUS*(6.6*N2N)
        COOL(5)=1.6E-26*OPLUS*(2.8*HEN)
        COOL(6)=1.6E-26*OPLUS*(5.8*O2N)
        COOL(7)=1.6E-26*OPLUS*(5.6*ON)
        COOL(8)=1.6E-26*HPLUS*(5.5*HEN)
        COOL(9)=1.6E-26*HPLUS*(1.9*HN)
        COOL(10)=1.6E-26*(0.36*OPLUS*HN)*SQRT(TN)
        COOL(11)=1.6E-26*(0.40*HPLUS*ON)*SQRT(TN)
        !... total cooling rate
        COOL(1)=COOL(2)+COOL(3)+COOL(4)+COOL(5)+COOL(6)+COOL(7)+
     >    COOL(8)+COOL(9)+COOL(10)+COOL(11)
      RETURN
      END
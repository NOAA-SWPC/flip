C............................... NEUT_HEATING.FOR ...................................
      SUBROUTINE NEUT_HEATING(JPR,  !.. Input: Turns printing on and off
     >                         IJ,  !.. Input: Altitude index
     >                        ALT,  !.. Input: Altitude (km)
     >                        RTS,  !.. Input: Reaction rates
     >                   TE,TI,TN,  !.. Input: Electron and ion temperatures
     >            OXN,O2N,N2N,HEN,  !.. Input: O, O2, N2, and He densities (cm-3)
     >                        N4S,  !.. Input: N4S should be 0.5*MSIS N density (cm-3)
     >                         NE,  !.. Input: electron density (cm-3)
     >       OXPLUS,O2PLUS,NOPLUS,  !.. Input: O+, O2+, NO+ densities (cm-3)
     >         HPLUS,N2PLUS,NPLUS,  !.. Input: N2+ and N+ densities (cm-3)
     >                NNO,N2D,N2P,  !.. Input: NO, N(2D), and N(2P) density (cm-3)
     >              N2A,OP2D,OP2P,  !.. Input: N2(A3), O+(2D), O+(2D)
     >                    O1D,O1S,  !.. Input: O(1D), O(1S)
     >                      NHEAT)  !.. OUTPUT: Total neutral heating rate
      USE PRODUCTION    !.. EUV, photoelectron, and auroral production, PHION
      IMPLICIT NONE
      INTEGER IJ,K          !.. loop control variables
      INTEGER JPR,JPT       !.. Turns on printing of production and loss
      DOUBLE PRECISION ALT              !.. Altitude (km)
      DOUBLE PRECISION RTS(99)          !.. Reaction rates
      DOUBLE PRECISION TE,TN,TI         !.. Electron and ion temperatures
	!..  H+, He+, O+, N2+, NO+, O2+, N+, 
      DOUBLE PRECISION HPLUS,OXPLUS,N2PLUS,NOPLUS,O2PLUS,NPLUS
      !.. O2,O,N2,NO,He,N4S
      DOUBLE PRECISION O2N,OXN,N2N,NNO,HEN,N4S
      !.. Ne, N(2P),N(2D),O+(2P),O+(2D), total minor ion densities
      DOUBLE PRECISION NE,N2P,N2D,OP2D,OP2P,NMINOR
      DOUBLE PRECISION O2DEL,O2B1,O2ADEL,O2ASIG  !.. O2 3-body states
      DOUBLE PRECISION DISNP            !.. Photodissociatin of N2 to give N+
      DOUBLE PRECISION UVDISN           !.. Photodissociatin of N2 to give N
	DOUBLE PRECISION N2A,O1D,O1S      !.. N2(A), O(1), O(1S) densities    
      DOUBLE PRECISION HR(99)           !.. Individual heating rates
      DOUBLE PRECISION BOD3,O2SS        !.. O2 3-body states
      DOUBLE PRECISION HRATE(22)        !.. Total heating rates
      DOUBLE PRECISION PO1DSR           !.. Schumann-Runge production of O(1D)
      DOUBLE PRECISION PEO1D,PO1SEL     !.. Photoelectron production of O(1D), O(1S)
      DOUBLE PRECISION PO1D,LO1D        !.. O(1D) producion and loss rates
      DOUBLE PRECISION COOL(22)         !.. Electron cooling rates
      DOUBLE PRECISION FOL              !.. Electron OX fine structure cooling rate
      DOUBLE PRECISION TLSS             !.. Electron cooling rate to neutrals
      DOUBLE PRECISION EHT              !.. Electron heating
      DOUBLE PRECISION NHEAT            !.. Neutral heating rate

      DATA O2DEL,O2B1,O2ADEL,O2ASIG/4*0.0/  !.. O2 3-body states
      DATA JPT/0/
      JPT=JPT+1

      PO1DSR=OTHPR1(3,IJ)
      PEO1D=PEXCIT(1,1,IJ)
      PO1SEL=PEXCIT(1,2,IJ)
      UVDISN=OTHPR1(1,IJ)

      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(121,150)
 150  FORMAT(/3X,'FOR HEATING EFFICIENCY: O2 dissociative energy'
     > ,1X,'losses '
     > ,/3X,'ALT   S-R    N2D+O2  N4S+O2  N++O2   E+O2+   N+O2+  TOTAL')
      HR(1)=5.08*PO1DSR
      HR(2)=5.08*N2D*O2N*RTS(16)
      HR(3)=5.08*RTS(7)*O2N*N4S
      HR(4)=5.08*(RTS(22)+RTS(25)+RTS(30)+RTS(65)+RTS(66))*O2N*NPLUS
      HR(5)=5.08*RTS(6)*NE*O2PLUS
      HR(6)=5.08*RTS(21)*O2PLUS*N4S
      HRATE(3)=HR(1)+HR(2)+HR(3)+HR(4)+HR(5)+HR(6)
      IF(JPR.GT.0) WRITE(121,88) ALT,(HR(K),K=1,6),HRATE(3)

      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(122,154)
 154  FORMAT(/2X,' Metastable kinetic heating'
     > ,/2X,'ALT   N2D+O   N2D+O2  N2D+e   O1D+N2  O1D+O2 O+2D+N2'
     > ,2X,'O+N2+  e+O+2D  O+O+2P  N2+O+2P e+O+2P   O+N2A   O+N2P'
     > ,3X,'total')
      HR(1)=2.38*N2D*OXN*RTS(15)
      HR(2)=1.84*N2D*O2N*RTS(16)
      HR(3)=2.38*N2D*NE*RTS(8)
      HR(4)=1.96*O1D*N2N*RTS(33)
      HR(5)=1.96*O1D*O2N*RTS(34)
      HR(6)=1.33*OP2D*RTS(19)*N2N
      HR(7)=2.3*N2PLUS*OXN*RTS(99)
      HR(8)=3.31*OP2D*NE*RTS(12)
      HR(9)=5.00*RTS(26)*OXN*OP2P
      HR(10)=3.02*RTS(20)*N2N*OP2P
      HR(11)=1.69*RTS(13)*NE*OP2P
      HR(12)=6.14*RTS(36)*OXN*N2A
      HR(13)=3.0*RTS(37)*N2P*OXN
      !.. n2 vibrational heating via o and o2 quenching not included
      HRATE(4)=HR(1)+HR(2)+HR(3)+HR(4)+HR(5)+HR(6)+HR(7)+HR(8)
     > +HR(9)+HR(10)+HR(11)+HR(12)+HR(13)
      IF(JPR.GT.0) WRITE(122,88) ALT,(HR(K),K=1,13),HRATE(4)
     >  ,O1D,O1S,N2D,N2P,N2A,NNO,OP2D,OP2P

      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(123,155)
 155   FORMAT(/2X,'Ground state kinetic heating'
     > ,/1X,'ALT  O+N2+   E+N2+   O2+N2+  E+NO+   E+O2+   N+O2+   N+O2+'
     > ,1X,'  NO+O2+  N2+O+   O2+O+   N+NO    N+O2    O2+N+   O2+N+'
     > ,1X,'   O+N+    total')
      HR(1)=0.70*RTS(10)*OXN*N2PLUS
      HR(2)=3.44*RTS(11)*NE*N2PLUS
      HR(3)=3.53*RTS(17)*O2N*N2PLUS
      HR(4)=0.90*NOPLUS*NE*RTS(5)
      HR(5)=4.42*RTS(6)*NE*O2PLUS
      HR(6)=4.19*RTS(21)*O2PLUS*N4S
      HR(7)=0.05*RTS(21)*O2PLUS*N4S
      HR(8)=2.81*RTS(23)*O2PLUS*NNO
      HR(9)=1.10*OXPLUS*N2N*RTS(3)
      HR(10)=1.55*OXPLUS*O2N*RTS(4)
      HR(11)=1.85*RTS(9)*N4S*NNO
      HR(12)=1.42*RTS(7)*O2N*N4S
      HR(13)=6.67*RTS(30)*O2N*NPLUS
      HR(14)=0.11*RTS(25)*O2N*NPLUS
      HR(15)=0.93*RTS(31)*OXN*NPLUS
      HRATE(5)=HR(1)+HR(2)+HR(3)+HR(4)+HR(5)+HR(6)+HR(7)+HR(8)
     > +HR(9)+HR(10)+HR(11)+HR(12)+HR(13)+HR(14)+HR(15)
      IF(JPR.GT.0) WRITE(123,88) ALT,(HR(K),K=1,15),HRATE(5)

      !.. Electron cooling
      NMINOR=NOPLUS+O2PLUS+N2PLUS
      CALL EHCRATS(125,JPR,ALT,TE,TI,TN,OXPLUS,HPLUS
     > ,NMINOR,OXN,N2N,O2N,HPLUS,EHT,TLSS,FOL,COOL)

      !.. Heating from 3 body collisions O+O+M
      HRATE(7)=(300/TN)**2*OXN**2*(O2N+N2N)*4.7E-33

      IF(COOL(5).LT.0) COOL(5)=0.0
      HRATE(8)=COOL(5)

      !..o2 dissoc heating. Energy left after producing O(1D)
      !.. HRATE(9)= S-R deposition rate. It is calculated in RSPRIM.FOR
      HRATE(9)=OTHPR1(5,IJ)-7.04*PO1DSR
	IF(HRATE(9).LT.0.0) HRATE(9)=0.0

      !.. EUV and photoelectron N2 dissociation rate to form N and N+
      DISNP= EUVION(3,4,IJ)+EUVION(3,5,IJ)+EUVION(3,6,IJ)
     >   +0.1*(PEPION(3,1,IJ)+PEPION(3,2,IJ)+PEPION(3,3,IJ))  !.. Rydberg diss       
      !.. kinetic heating from dissociation of N2 to N + N and N + N+
      !.. by EUV and photoelectrons. EuV dissociation assumed to release
	!.. 1 eV of KE and photoelectron dissociation assumed to release 2 eV
      HRATE(12)=1.0*UVDISN+2.0*DISNP

      !.. Total neutral heating rate
      NHEAT=HRATE(3)+HRATE(4)+HRATE(5)+HRATE(7)+HRATE(8)+HRATE(12)

      IF(JPT.EQ.1.AND.JPR.GT.0) WRITE(126,159)
 159   FORMAT(/2X,'Total neutral heating rates'
     > ,/3X,'ALT   O2_diss     Meta     Ground    3body    Electron'
     > ,2X,'SRO2dis   UVN2dis  Total')
      IF(JPR.GT.0) WRITE(126,'(F7.2,1P,22E10.2)') ALT,HRATE(3),
     >  HRATE(4),HRATE(5),HRATE(7),HRATE(8),HRATE(9),HRATE(12),NHEAT

 88   FORMAT(F6.1,1P,22E8.1)
 89   FORMAT(I4,1P,8E8.1,0PF6.0,F4.1)

      RETURN
       END
C::::::::::::::::::::::: EHCRATS :::::::::::::::::::::::::::::::::::
C.... Subroutine for printing thermal electron cooling rates
C.... This function needs to be eliminated eventually and replaced
C.... by GET_loss. It is used to print cooling rates in the Hyyyyddd
C.... file of the FLIP run and also in the FLIPRINT phase 
      SUBROUTINE EHCRATS(ID,JPR,ALT,TE,TI,TN,OPLS,HPLUS,NMIN,OX,N2,O2,HN
     > ,EHT,TLSS,LFO,COOL)
      IMPLICIT DOUBLE PRECISION(A-H,K,L,N,O-Z)
      DIMENSION FD(9),COOL(22) !.. heating and cooling rates
C
      
      IF(ID.NE.IDS.AND.ID.NE.0.AND.JPR.GT.0) WRITE(ID,757)
      IDS=ID
 757  FORMAT(/,8X,'Thermal electron cooling(C_) and heating(H_) rates'
     > /,6X,'ALT  Te   Ti   Tn     C_ei      C_rN2     C_rO2    C_vN2',
     > '    C_vO2     C_fsO      C_O1D     C_tot      H_tot   H_N2D'
     > '  HOp2D')
        TLSS=0.0
        LFO=0.0

        !.. Don't calculate cooling rates if neutral densities too low
        IF(O2.LT.1.0E-10) RETURN     
        SQTE=SQRT(TE)
        NE=OPLS+HPLUS+NMIN
        NI=OPLS/16.0+HPLUS+NMIN/30.0

        !.. Electron- ion cooling rate KEI. The Coulomb logarithm from Itikawa
        !..  JATP 1975, page 1001. Previously assumed to be 15.0
        DEBYE = 1.602E-12 * SQRT(4.0*3.1415926*NE/(1.381E-16*TE))
        COULOG = 0.43086 + LOG(TE) - LOG(DEBYE)
        KEI=1.232E-17*NE*NI*(TE-TI)/(TE**1.5)/1.6E-12    !.. e - i cooling
        KEI= KEI*COULOG/15.0                     !.. Itikawa correction

        RTE2=1.0/TE
        RTNJ=1.0/TN
 13     TDIF=TE-TN
        TDIFRT=TDIF*RTE2*RTNJ
        LEN2=1.77E-19*NE*N2*(1.0-1.2E-4*TE)*TE*TDIF
        LEO2=1.21E-18*NE*O2*(1.0+3.6E-2*SQTE)*SQTE*TDIF
        LEO=7.9E-19*NE*OX*(1.0+5.7E-4*TE)*SQTE*TDIF
        LEH=9.63E-16*NE*HN*(1.0-1.35E-4*TE)*SQTE*TDIF
C
C###  ROTATIONAL loss RATES B&K 268 ####
        LRN2=2.0E-14*NE*N2*TDIF/SQTE
        LRO2=7.0E-14*NE*O2*TDIF/SQTE
C
C+++  N2 VIB loss RATES B&K P268 AND S&N P364  +++++
        EF=1.06E+4+7.51E+3*TANH(1.1E-3*(TE-1800.0))
        GE=3300.0+1.233*(TE-1000.0)-2.056E-4*(TE-1000.0)*(TE-4000.0)
        EXH=5.0E-4*RTE2*EF*(TE-2000.0)
        IF(EXH.GT.70.0)EXH=70.0
        EXH2=GE*TDIFRT
        IF(EXH2.GT.70.0)EXH2=70.0
        LVN2=-2.99E-12*NE*N2*EXP(EXH)
     > *(EXP(-EXH2)-1.0)

C$$$$ O2 VIB loss RATE S&N P364  $$$$
        HS=3300.0-839.0*SIN(1.91E-4*(TE-2700.0))
        EXH=1.4286E-3*RTE2*HS*(TE-700.0)
        IF(EXH.GT.70.0)EXH=70.0
        EXH2=2770.0*TDIFRT
        IF(EXH2.GT.70.0)EXH2=70.0
        LVO2=-5.196E-13*NE*O2*EXP(EXH2)
     >        *(EXP(-EXH2)-1.0)
C
C;;;;;  FINE STRUCTURE EXCITATIONS S&N P365 ;;;;;;
C;;;;; NOTE THAT D1,D2,E1 ETC. MAY NEED TO BE CHANGED
C;;;; NOTE THAT A TERM IS ADDED TO NE AT LOW ALTS FOR O2+ , NO+
        HD1=228.0*RTNJ
        IF(HD1.GT.70.0)HD1=70.0
        D1=EXP(-HD1)
        HD2=326.0*RTNJ
        IF(HD2.GT.70.0)HD2=70.0
        D2=EXP(-HD2)
        HE1=228.0*RTE2
        IF(HE1.GT.70.0)HE1=70.0
        E1=EXP(-HE1)
        HE2=326.0*RTE2
        IF(HE2.GT.70.0)HE2=70.0
        E2=EXP(-HE2)
        HE3=98.0*RTE2
        IF(HE3.GT.70.0)HE3=70.0
        E3=EXP(-HE3)*D1
        LF1=8.49E-6*TE**0.519*(0.02*(D1-E1)-5.91E-9*
     >     TDIF*(2.019*D1+(228.0*RTE2+2.019)*E1))
        LF2=7.7E-6*TE**0.3998*(0.028*(D2-E2)-5.91E-9*
     >     TDIF*(1.8998*D2+(326.0*RTE2+1.8998)*E2))
        LF3=2.22E-7*TE**0.768*(0.008*(D2-E3)-5.91E-9*
     >     TDIF*(2.268*D2+(98.0*RTE2+2.268)*E3))
        ZFO=5.0+3.0*D1+D2
        LFO=-8.629E-6*NE*OX*(LF1+LF2+LF3)/ZFO
C....... LFO=0.333*LFO
C++++  EXCITATION OF O(1D) BY PE'S S&N P365. - SMALL ONLY CALC FOR PRINT
        DE=2.4E4+(TE-1500.0)*(0.3-1.947E-5*(TE-4000.0))
        EXH=3.3333E-4*DE*(TE-3000.0)*RTE2
        IF(EXH.GT.70.0)EXH=70.0
        EXH2=22713.0*TDIFRT
        IF(EXH2.GT.70.0)EXH2=70.0
        LF1D=-1.57E-12*NE*OX*EXP(EXH)*(EXP(-EXH2)-1.0)

C......... add other losses to electron-ion
      TLSS=LEN2+LEO2+LEO+LEH+LRO2+LRN2+LVO2+LVN2+LF1D
	IF(TLSS.LT.0) TLSS=0.0

C......... electron heating from N(2D) .....
      DO III=1,8
        FD(III)=0.0
      ENDDO
      IF(JPR.GT.0) WRITE(ID,54) ALT,INT(TE),INT(TI),INT(TN),KEI
     > ,LRN2,LRO2,LVN2,LVO2,LFO,LF1D,TLSS+LFO+KEI,EHT,FD(8)*2.4,
     > FD(6)*3.31
 54    FORMAT(F9.1,3I5,1P,22E10.2)

      COOL(1)=EHT
      COOL(2)=KEI
      COOL(3)=LVN2
      COOL(4)=LFO+LF1D
      COOL(5)=TLSS+LFO+KEI
      RETURN
      END

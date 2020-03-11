C.................... RSJACA.FOR;3 ..........12-MAR-1993 09:02:35.48 
      SUBROUTINE BDSLV(N,M,S,KS,B,X,WORK,NFLAG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION MUL(50000)
      DIMENSION S(N,1),WORK(N,1),SCALE(N),INDEX(N),B(N),
     >  X(N)
      NFLAG=0
      IW=2*M+1
      M1=M+1
C
      CALL BNDX(N,M,S,KS,B,X,WORK,SCALE,INDEX,MUL,NFLAG,IW,M1)
      RETURN
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE BNDX(N,M,A,KS,B,X,S,SCALE,INDEX,MUL,NFLAG,IW,M1)
      USE module_IO,ONLY: PRUNIT
      USE module_input_parameters,ONLY:sw_ERSTOP_flip,mype
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION MUL(M1,N)
      DIMENSION A(N,IW),B(N),X(N),SCALE(N),INDEX(N),S(N,IW)
      KDIM=0
103   IF(KS.NE.0) GO TO 126
      IBW=2*M+1
C ...... modification here by PR Aug 91
      IF(M.GT.N-1) THEN
!dbg20120306         WRITE(6,918)
!sms$ignore begin
         WRITE(PRUNIT,918) mype
!sms$ignore end
918      FORMAT('    IN BDSLV &&&&&&& BANDWIDTH IS TOO LARGE',i10)
         NFLAG=3
         RETURN
      ENDIF
C
      DO 105 J=1,IBW
      DO 104 I=1,N
      S(I,J)=A(I,J)
104   CONTINUE
105   CONTINUE
C
C
      NP1=N+1
      MP1=M+1
C
      DO 108 I=1,M
      MP1PI=MP1+I
      IF(MP1PI-IBW) 106,106,109
106   NP1MI=NP1-I
      DO 107 J=MP1PI,IBW
      S(I,J)=0.
      S(NP1MI,J)=0.
107   CONTINUE
108   CONTINUE
C
109   DO 115 I=1,N
      BIG=0.
      DO 111 J=1,IBW
      IF(BIG-ABS(S(I,J))) 110,111,111
110   BIG=ABS(S(I,J))
111   CONTINUE
      IF(BIG) 114,112,114
!dbg20120306  112 WRITE(6,919)I
!sms$ignore begin
  112 WRITE(PRUNIT,919)I,mype
!sms$ignore end
919   FORMAT('    IN BDSLV, ROW',I6,' IS ZERO IN INPUT MATRIX',i10)
      NFLAG=2
!dbg20140610: code must stop if this error happens!
!nm20170111: commented out to get WAM-IPE keep going at the initial trial stage by loosening error stop.
!      sw_ERSTOP_flip=1
      RETURN
114   SCALE(I)=1./BIG
115   CONTINUE
C
      LOW=M
      NM1=N-1
      DO 125 K=1,NM1
      LOW=MIN0(LOW+1,N)
      BIG=0.
      DO 117 I=K,LOW
      SIZE=ABS(S(I,1))*SCALE(I)
      IF(SIZE-BIG) 117,117,116
116   BIG=SIZE
      IPIV=I
117   CONTINUE
      IF(BIG) 119,118,119
119   INDEX(K)=IPIV
      IF(IPIV-K) 120,122,120
120   SCALE(IPIV)=SCALE(K)
      DO 121 J=1,IBW
      TEMP=S(K,J)
      S(K,J)=S(IPIV,J)
      S(IPIV,J)=TEMP
121   CONTINUE
122   PIVOT=S(K,1)
      KOUNT=0
      KP1=K+1
      DO 124 I=KP1,LOW
      KOUNT=KOUNT+1
      FACT=-S(I,1)/PIVOT
      MUL(KOUNT,K)=FACT
      IDIM=KOUNT*K
      IF(KDIM.LT.IDIM) KDIM=IDIM
      DO 123 J=2,IBW
      S(I,J-1)=S(I,J)+FACT*S(K,J)
123   CONTINUE
      S(I,IBW)=0.
124   CONTINUE
125   CONTINUE
C
      IF(S(N,1)) 126,118,126
C  ..... PR mod in Aug 91
!dbg20120306: 118   WRITE(6,917)
118   WRITE(PRUNIT,917) mype
917   FORMAT('  IN BDSLV &&&&&&&&   ZERO PIVOT ELEMENT',i10)
      NFLAG=1
      RETURN
C
126   DO 127 I=1,N
      X(I)=B(I)
127   CONTINUE
      LOW=M
      DO 130 K=1,NM1
      IPIV=INDEX(K)
      IF(IPIV.EQ.K) GO TO 128
      TEMP=X(K)
      X(K)=X(IPIV)
      X(IPIV)=TEMP
128   KOUNT=0
      KP1=K+1
      LOW=MIN0(LOW+1,N)
      DO 129 I=KP1,LOW
      KOUNT=KOUNT+1
      X(I)=X(I)+MUL(KOUNT,K)*X(K)
129   CONTINUE
130   CONTINUE
      X(N)=X(N)/S(N,1)
      LOW=1
      DO 132 IBACK=1,NM1
      I=N-IBACK
      SUM=0.
      LOW=MIN0(LOW+1,IBW)
      DO 131 JBACK=2,LOW
      J=I+JBACK-1
      SUM=SUM+S(I,JBACK)*X(J)
131   CONTINUE
      X(I)=(X(I)-SUM)/S(I,1)
132   CONTINUE
      IF(KDIM.GT.50000) WRITE(6,90) KDIM
 90   FORMAT(' WARNING!! DIMENSION OF MUL IN BDSLV TOO SMALL. KDIM=',I9)
      RETURN
      END

 
C*********************************************************************
 
C...PYX4JT
C...Selects the kinematical variables of four-jet events.
 
      SUBROUTINE PYX4JT(NJET,CUT,KFL,ECM,KFLN,X1,X2,X4,X12,X14)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
C...Local arrays.
      DIMENSION WTA(4),WTB(4),WTC(4),WTD(4),WTE(4)
 
C...Common constants. Colour factors for QCD and Abelian gluon theory.
      PMQ=PYMASS(KFL)
      QME=(2D0*PMQ/ECM)**2
      CT=LOG(1D0/CUT-5D0)
      IF(MSTJ(109).EQ.0) THEN
        CF=4D0/3D0
        CN=3D0
        TR=2.5D0
      ELSE
        CF=1D0
        CN=0D0
        TR=15D0
      ENDIF
 
C...Choice of process (qqbargg or qqbarqqbar).
  100 NJET=4
      IT=1
      IF(PARJ(155).GT.PYR(0)) IT=2
      IF(MSTJ(101).LE.-3) IT=-MSTJ(101)-2
      IF(IT.EQ.1) WTMX=0.7D0/CUT**2
      IF(IT.EQ.1.AND.MSTJ(109).EQ.2) WTMX=0.6D0/CUT**2
      IF(IT.EQ.2) WTMX=0.1125D0*CF*TR/CUT**2
      ID=1
 
C...Sample the five kinematical variables (for qqgg preweighted in y34).
  110 Y134=3D0*CUT+(1D0-6D0*CUT)*PYR(0)
      Y234=3D0*CUT+(1D0-6D0*CUT)*PYR(0)
      IF(IT.EQ.1) Y34=(1D0-5D0*CUT)*EXP(-CT*PYR(0))
      IF(IT.EQ.2) Y34=CUT+(1D0-6D0*CUT)*PYR(0)
      IF(Y34.LE.Y134+Y234-1D0.OR.Y34.GE.Y134*Y234) GOTO 110
      VT=PYR(0)
      CP=COS(PARU(1)*PYR(0))
      Y14=(Y134-Y34)*VT
      Y13=Y134-Y14-Y34
      VB=Y34*(1D0-Y134-Y234+Y34)/((Y134-Y34)*(Y234-Y34))
      Y24=0.5D0*(Y234-Y34)*(1D0-4D0*SQRT(MAX(0D0,VT*(1D0-VT)*
     &VB*(1D0-VB)))*CP-(1D0-2D0*VT)*(1D0-2D0*VB))
      Y23=Y234-Y34-Y24
      Y12=1D0-Y134-Y23-Y24
      IF(MIN(Y12,Y13,Y14,Y23,Y24).LE.CUT) GOTO 110
      Y123=Y12+Y13+Y23
      Y124=Y12+Y14+Y24
 
C...Calculate matrix elements for qqgg or qqqq process.
      IC=0
      WTTOT=0D0
  120 IC=IC+1
      IF(IT.EQ.1) THEN
        WTA(IC)=(Y12*Y34**2-Y13*Y24*Y34+Y14*Y23*Y34+3D0*Y12*Y23*Y34+
     &  3D0*Y12*Y14*Y34+4D0*Y12**2*Y34-Y13*Y23*Y24+2D0*Y12*Y23*Y24-
     &  Y13*Y14*Y24-2D0*Y12*Y13*Y24+2D0*Y12**2*Y24+Y14*Y23**2+2D0*Y12*
     &  Y23**2+Y14**2*Y23+4D0*Y12*Y14*Y23+4D0*Y12**2*Y23+2D0*Y12*Y14**2+
     &  2D0*Y12*Y13*Y14+4D0*Y12**2*Y14+2D0*Y12**2*Y13+2D0*Y12**3)/
     &  (2D0*Y13*Y134*Y234*Y24)+(Y24*Y34+Y12*Y34+Y13*Y24-
     &  Y14*Y23+Y12*Y13)/(Y13*Y134**2)+2D0*Y23*(1D0-Y13)/
     &  (Y13*Y134*Y24)+Y34/(2D0*Y13*Y24)
        WTB(IC)=(Y12*Y24*Y34+Y12*Y14*Y34-Y13*Y24**2+Y13*Y14*Y24+2D0*Y12*
     &  Y14*Y24)/(Y13*Y134*Y23*Y14)+Y12*(1D0+Y34)*Y124/(Y134*Y234*Y14*
     &  Y24)-(2D0*Y13*Y24+Y14**2+Y13*Y23+2D0*Y12*Y13)/(Y13*Y134*Y14)+
     &  Y12*Y123*Y124/(2D0*Y13*Y14*Y23*Y24)
        WTC(IC)=-(5D0*Y12*Y34**2+2D0*Y12*Y24*Y34+2D0*Y12*Y23*Y34+
     &  2D0*Y12*Y14*Y34+2D0*Y12*Y13*Y34+4D0*Y12**2*Y34-Y13*Y24**2+
     &  Y14*Y23*Y24+Y13*Y23*Y24+Y13*Y14*Y24-Y12*Y14*Y24-Y13**2*Y24-
     &  3D0*Y12*Y13*Y24-Y14*Y23**2-Y14**2*Y23+Y13*Y14*Y23-
     &  3D0*Y12*Y14*Y23-Y12*Y13*Y23)/(4D0*Y134*Y234*Y34**2)+
     &  (3D0*Y12*Y34**2-3D0*Y13*Y24*Y34+3D0*Y12*Y24*Y34+
     &  3D0*Y14*Y23*Y34-Y13*Y24**2-Y12*Y23*Y34+6D0*Y12*Y14*Y34+
     &  2D0*Y12*Y13*Y34-2D0*Y12**2*Y34+Y14*Y23*Y24-3D0*Y13*Y23*Y24-
     &  2D0*Y13*Y14*Y24+4D0*Y12*Y14*Y24+2D0*Y12*Y13*Y24+
     &  3D0*Y14*Y23**2+2D0*Y14**2*Y23+2D0*Y14**2*Y12+
     &  2D0*Y12**2*Y14+6D0*Y12*Y14*Y23-2D0*Y12*Y13**2-
     &  2D0*Y12**2*Y13)/(4D0*Y13*Y134*Y234*Y34)
        WTC(IC)=WTC(IC)+(2D0*Y12*Y34**2-2D0*Y13*Y24*Y34+Y12*Y24*Y34+
     &  4D0*Y13*Y23*Y34+4D0*Y12*Y14*Y34+2D0*Y12*Y13*Y34+2D0*Y12**2*Y34-
     &  Y13*Y24**2+3D0*Y14*Y23*Y24+4D0*Y13*Y23*Y24-2D0*Y13*Y14*Y24+
     &  4D0*Y12*Y14*Y24+2D0*Y12*Y13*Y24+2D0*Y14*Y23**2+4D0*Y13*Y23**2+
     &  2D0*Y13*Y14*Y23+2D0*Y12*Y14*Y23+4D0*Y12*Y13*Y23+2D0*Y12*Y14**2+
     &  4D0*Y12**2*Y13+4D0*Y12*Y13*Y14+2D0*Y12**2*Y14)/
     &  (4D0*Y13*Y134*Y24*Y34)-(Y12*Y34**2-2D0*Y14*Y24*Y34-
     &  2D0*Y13*Y24*Y34-Y14*Y23*Y34+Y13*Y23*Y34+Y12*Y14*Y34+
     &  2D0*Y12*Y13*Y34-2D0*Y14**2*Y24-4D0*Y13*Y14*Y24-
     &  4D0*Y13**2*Y24-Y14**2*Y23-Y13**2*Y23+Y12*Y13*Y14-
     &  Y12*Y13**2)/(2D0*Y13*Y34*Y134**2)+(Y12*Y34**2-
     &  4D0*Y14*Y24*Y34-2D0*Y13*Y24*Y34-2D0*Y14*Y23*Y34-
     &  4D0*Y13*Y23*Y34-4D0*Y12*Y14*Y34-4D0*Y12*Y13*Y34-
     &  2D0*Y13*Y14*Y24+2D0*Y13**2*Y24+2D0*Y14**2*Y23-
     &  2D0*Y13*Y14*Y23-Y12*Y14**2-6D0*Y12*Y13*Y14-
     &  Y12*Y13**2)/(4D0*Y34**2*Y134**2)
        WTTOT=WTTOT+Y34*CF*(CF*WTA(IC)+(CF-0.5D0*CN)*WTB(IC)+
     &  CN*WTC(IC))/8D0
      ELSE
        WTD(IC)=(Y13*Y23*Y34+Y12*Y23*Y34-Y12**2*Y34+Y13*Y23*Y24+2D0*Y12*
     &  Y23*Y24-Y14*Y23**2+Y12*Y13*Y24+Y12*Y14*Y23+Y12*Y13*Y14)/(Y13**2*
     &  Y123**2)-(Y12*Y34**2-Y13*Y24*Y34+Y12*Y24*Y34-Y14*Y23*Y34-Y12*
     &  Y23*Y34-Y13*Y24**2+Y14*Y23*Y24-Y13*Y23*Y24-Y13**2*Y24+Y14*
     &  Y23**2)/(Y13**2*Y123*Y134)+(Y13*Y14*Y12+Y34*Y14*Y12-Y34**2*Y12+
     &  Y13*Y14*Y24+2D0*Y34*Y14*Y24-Y23*Y14**2+Y34*Y13*Y24+Y34*Y23*Y14+
     &  Y34*Y13*Y23)/(Y13**2*Y134**2)-(Y34*Y12**2-Y13*Y24*Y12+Y34*Y24*
     &  Y12-Y23*Y14*Y12-Y34*Y14*Y12-Y13*Y24**2+Y23*Y14*Y24-Y13*Y14*Y24-
     &  Y13**2*Y24+Y23*Y14**2)/(Y13**2*Y134*Y123)
        WTE(IC)=(Y12*Y34*(Y23-Y24+Y14+Y13)+Y13*Y24**2-Y14*Y23*Y24+Y13*
     &  Y23*Y24+Y13*Y14*Y24+Y13**2*Y24-Y14*Y23*(Y14+Y23+Y13))/(Y13*Y23*
     &  Y123*Y134)-Y12*(Y12*Y34-Y23*Y24-Y13*Y24-Y14*Y23-Y14*Y13)/(Y13*
     &  Y23*Y123**2)-(Y14+Y13)*(Y24+Y23)*Y34/(Y13*Y23*Y134*Y234)+
     &  (Y12*Y34*(Y14-Y24+Y23+Y13)+Y13*Y24**2-Y23*Y14*Y24+Y13*Y14*Y24+
     &  Y13*Y23*Y24+Y13**2*Y24-Y23*Y14*(Y14+Y23+Y13))/(Y13*Y14*Y134*
     &  Y123)-Y34*(Y34*Y12-Y14*Y24-Y13*Y24-Y23*Y14-Y23*Y13)/(Y13*Y14*
     &  Y134**2)-(Y23+Y13)*(Y24+Y14)*Y12/(Y13*Y14*Y123*Y124)
        WTTOT=WTTOT+CF*(TR*WTD(IC)+(CF-0.5D0*CN)*WTE(IC))/16D0
      ENDIF
 
C...Permutations of momenta in matrix element. Weighting.
  130 IF(IC.EQ.1.OR.IC.EQ.3.OR.ID.EQ.2.OR.ID.EQ.3) THEN
        YSAV=Y13
        Y13=Y14
        Y14=YSAV
        YSAV=Y23
        Y23=Y24
        Y24=YSAV
        YSAV=Y123
        Y123=Y124
        Y124=YSAV
      ENDIF
      IF(IC.EQ.2.OR.IC.EQ.4.OR.ID.EQ.3.OR.ID.EQ.4) THEN
        YSAV=Y13
        Y13=Y23
        Y23=YSAV
        YSAV=Y14
        Y14=Y24
        Y24=YSAV
        YSAV=Y134
        Y134=Y234
        Y234=YSAV
      ENDIF
      IF(IC.LE.3) GOTO 120
      IF(ID.EQ.1.AND.WTTOT.LT.PYR(0)*WTMX) GOTO 110
      IC=5
 
C...qqgg events: string configuration and event type.
      IF(IT.EQ.1) THEN
        IF(MSTJ(109).EQ.0.AND.ID.EQ.1) THEN
          PARJ(156)=Y34*(2D0*(WTA(1)+WTA(2)+WTA(3)+WTA(4))+4D0*(WTC(1)+
     &    WTC(2)+WTC(3)+WTC(4)))/(9D0*WTTOT)
          IF(WTA(2)+WTA(4)+2D0*(WTC(2)+WTC(4)).GT.PYR(0)*(WTA(1)+WTA(2)+
     &    WTA(3)+WTA(4)+2D0*(WTC(1)+WTC(2)+WTC(3)+WTC(4)))) ID=2
          IF(ID.EQ.2) GOTO 130
        ELSEIF(MSTJ(109).EQ.2.AND.ID.EQ.1) THEN
          PARJ(156)=Y34*(WTA(1)+WTA(2)+WTA(3)+WTA(4))/(8D0*WTTOT)
          IF(WTA(2)+WTA(4).GT.PYR(0)*(WTA(1)+WTA(2)+WTA(3)+WTA(4))) ID=2
          IF(ID.EQ.2) GOTO 130
        ENDIF
        MSTJ(120)=3
        IF(MSTJ(109).EQ.0.AND.0.5D0*Y34*(WTC(1)+WTC(2)+WTC(3)+
     &  WTC(4)).GT.PYR(0)*WTTOT) MSTJ(120)=4
        KFLN=21
 
C...Mass cuts. Kinematical variables out.
        IF(Y12.LE.CUT+QME) NJET=2
        IF(NJET.EQ.2) GOTO 150
        Q12=0.5D0*(1D0-SQRT(1D0-QME/Y12))
        X1=1D0-(1D0-Q12)*Y234-Q12*Y134
        X4=1D0-(1D0-Q12)*Y134-Q12*Y234
        X2=1D0-Y124
        X12=(1D0-Q12)*Y13+Q12*Y23
        X14=Y12-0.5D0*QME
        IF(Y134*Y234/((1D0-X1)*(1D0-X4)).LE.PYR(0)) NJET=2
 
C...qqbarqqbar events: string configuration, choose new flavour.
      ELSE
        IF(ID.EQ.1) THEN
          WTR=PYR(0)*(WTD(1)+WTD(2)+WTD(3)+WTD(4))
          IF(WTR.LT.WTD(2)+WTD(3)+WTD(4)) ID=2
          IF(WTR.LT.WTD(3)+WTD(4)) ID=3
          IF(WTR.LT.WTD(4)) ID=4
          IF(ID.GE.2) GOTO 130
        ENDIF
        MSTJ(120)=5
        PARJ(156)=CF*TR*(WTD(1)+WTD(2)+WTD(3)+WTD(4))/(16D0*WTTOT)
  140   KFLN=1+INT(5D0*PYR(0))
        IF(KFLN.NE.KFL.AND.0.2D0*PARJ(156).LE.PYR(0)) GOTO 140
        IF(KFLN.EQ.KFL.AND.1D0-0.8D0*PARJ(156).LE.PYR(0)) GOTO 140
        IF(KFLN.GT.MSTJ(104)) NJET=2
        PMQN=PYMASS(KFLN)
        QMEN=(2D0*PMQN/ECM)**2
 
C...Mass cuts. Kinematical variables out.
        IF(Y24.LE.CUT+QME.OR.Y13.LE.1.1D0*QMEN) NJET=2
        IF(NJET.EQ.2) GOTO 150
        Q24=0.5D0*(1D0-SQRT(1D0-QME/Y24))
        Q13=0.5D0*(1D0-SQRT(1D0-QMEN/Y13))
        X1=1D0-(1D0-Q24)*Y123-Q24*Y134
        X4=1D0-(1D0-Q24)*Y134-Q24*Y123
        X2=1D0-(1D0-Q13)*Y234-Q13*Y124
        X12=(1D0-Q24)*((1D0-Q13)*Y14+Q13*Y34)+Q24*((1D0-Q13)*Y12+
     &  Q13*Y23)
        X14=Y24-0.5D0*QME
        X34=(1D0-Q24)*((1D0-Q13)*Y23+Q13*Y12)+Q24*((1D0-Q13)*Y34+
     &  Q13*Y14)
        IF(PMQ**2+PMQN**2+MIN(X12,X34)*ECM**2.LE.
     &  (PARJ(127)+PMQ+PMQN)**2) NJET=2
        IF(Y123*Y134/((1D0-X1)*(1D0-X4)).LE.PYR(0)) NJET=2
      ENDIF
  150 IF(MSTJ(101).LE.-2.AND.NJET.EQ.2) GOTO 100
 
      RETURN
      END

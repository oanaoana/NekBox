c-----------------------------------------------------------------------
      subroutine cbcmesh
C
C     Generate boundary conditions (CBC arrays) for mesh solver
C
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
C
      CHARACTER CBM*1,CBF*3,CBT*3,CB*3
C
      IFLD   = 0
      NFACE  = 2*NDIM
      IFMELT = .FALSE.
      IF (IFTMSH(IFLD)) IFMELT=.TRUE.
C
      DO 100 IEL=1,NELT
      DO 100 IFC=1,NFACE
         CBM = CBC(IFC,IEL,0)
         CBF = CBC(IFC,IEL,1)
         CBT = CBC(IFC,IEL,2)
C
         IF (CBT(1:1).EQ.'M') THEN
             IFLD = 2
             CB   = CBT
             GOTO 200
         ENDIF
         IF (CBF(1:1).EQ.'M' .OR. CBF(1:1).EQ.'m') THEN
             IFLD = 1
             CB   = CBF
C             IF (CBF.EQ.'mv ' .AND. CBT.EQ.'E  ') THEN
C                 IFTMSH(0)=.TRUE.
C             ENDIF
             IF (CBF.EQ.'mv ' .AND. CBM.EQ.'+'  ) THEN
                CB = 'mvn'
                CBC(IFC,IEL,1) = CB
             ENDIF
             GOTO 200
         ENDIF
         IF (CBF(1:1).EQ.'V' .OR. CBF(1:1).EQ.'v' .OR.
     $       CBF(1:1).EQ.'W' ) THEN
             IFLD = 1
             CB   = 'FIX'
             IF (IFMELT .OR. CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IF (CBT.EQ.'T  ' .OR. CBT.EQ.'t  ') THEN
             IFLD = 2
             CB   = 'FIX'
             IF (CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IF (CBF.EQ.'P  ' .OR. CBF.EQ.'E  ') THEN
             IFLD = 1
             CB   = CBF
             IF (CBM.EQ.'-') CB='FIX'
             IF (CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IF (CBT.EQ.'P  ' .OR. CBT.EQ.'E  ') THEN
             IFLD = 2
             CB   = CBT
             IF (CBM.EQ.'-') CB='FIX'
             IF (CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IFLD = 1
         IF (CBF.EQ.'   ') IFLD = 2
         CB   = 'SYM'
         IF (CBM.EQ.'-') CB = 'FIX'
C
  200    CBC(IFC,IEL,0) = CB
         DO 250 I=1,5
  250    BC(I,IFC,IEL,0)=BC(I,IFC,IEL,IFLD)
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine admeshv
C
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
C
      COMMON /SCRUZ/ FM1(LX1,LY1,LZ1,LELT)
     $             , FM2(LX1,LY1,LZ1,LELT)
     $             , FM3(LX1,LY1,LZ1,LELT)
     $             , PHI(LX1,LY1,LZ1,LELT)
C
      NTOT1=NX1*NY1*NZ1*NELV
C
      CALL RZERO (FM1,NTOT1)
      CALL RZERO (FM2,NTOT1)
      CALL RZERO (FM3,NTOT1)
C
      CALL DIVWS (FM1,VX,PHI,NELV,1)
      CALL DIVWS (FM2,VY,PHI,NELV,2)
      CALL ADD2  (BFX,FM1,NTOT1)
      CALL ADD2  (BFY,FM2,NTOT1)
      IF (NDIM.EQ.3) THEN
         CALL DIVWS (FM3,VZ,PHI,NELV,3)
         CALL ADD2  (BFZ,FM3,NTOT1)
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine admesht
C
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
C
      COMMON /SCRUZ/ FMT(LX1,LY1,LZ1,LELT)
     $             , PHI(LX1,LY1,LZ1,LELT)
C
      IFLD = 0
      NEL  = NELFLD(IFLD)
      NTOT1= NX1*NY1*NZ1*NEL
C
      CALL RZERO   (FMT,NTOT1)
      CALL DIVWS   (FMT,T(1,1,1,1,IFIELD-1),PHI,NEL,1)
      CALL ADDCOL3 (BQ(1,1,1,1,IFIELD-1),FMT,VTRANS(1,1,1,1,IFIELD),
     &              NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine divws (fms,sfv,phi,nel,idir)
C
      include 'SIZE'
      include 'GEOM'
      include 'MASS'
      include 'MVGEOM'
      include 'WZ'
      include 'INPUT'
C
      COMMON /SCRSF/ PHR(LX1,LY1,LZ1,LELT)
     $             , PHS(LX1,LY1,LZ1,LELT)
     $             , PHT(LX1,LY1,LZ1,LELT)

      DIMENSION FMS(LX1,LY1,LZ1,1)
     $        , SFV(LX1,LY1,LZ1,1)
     $        , PHI(LX1,LY1,LZ1,1)
C
      NXYZ1 = NX1*NY1*NZ1
      NTOT1 = NXYZ1*NEL
C
      CALL COL3    (PHI,SFV,WX,NTOT1)
      CALL URST    (PHI,PHR,PHS,PHT,NEL)
      CALL ADDCOL3 (FMS,RXM1,PHR,NTOT1)
      CALL ADDCOL3 (FMS,SXM1,PHS,NTOT1)
      IF (NDIM.EQ.3) CALL ADDCOL3 (FMS,TXM1,PHT,NTOT1)
C
      CALL COL3    (PHI,SFV,WY,NTOT1)
      CALL URST    (PHI,PHR,PHS,PHT,NEL)
      CALL ADDCOL3 (FMS,RYM1,PHR,NTOT1)
      CALL ADDCOL3 (FMS,SYM1,PHS,NTOT1)
      IF (NDIM.EQ.3) CALL ADDCOL3 (FMS,TYM1,PHT,NTOT1)
C      
      IF (NDIM.EQ.3) THEN
         CALL COL3    (PHI,SFV,WZ,NTOT1)
         CALL URST    (PHI,PHR,PHS,PHT,NEL)
         CALL ADDCOL3 (FMS,RZM1,PHR,NTOT1)
         CALL ADDCOL3 (FMS,SZM1,PHS,NTOT1)
         CALL ADDCOL3 (FMS,TZM1,PHT,NTOT1)
      ENDIF
C
      CALL COL2    (FMS,BM1,NTOT1)
      CALL INVCOL2 (FMS,JACM1,NTOT1)
C
      IF (IFAXIS) CALL AXIFMS (FMS,SFV,PHI,NEL,IDIR)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine axifms (fms,sfv,phi,nel,idir)
C
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
      include 'MVGEOM'
      include 'WZ'
      COMMON /SCRSF/ PHR(LX1,LY1,LZ1,LELT)
     $             , PHS(LX1,LY1,LZ1,LELT)
     $             , PHT(LX1,LY1,LZ1,LELT)
C
      DIMENSION FMS(LX1,LY1,LZ1,1)
     $        , PHI(LX1,LY1,LZ1,1)
     $        , SFV(LX1,LY1,LZ1,1)
     $        , WYS(LX1)
      EQUIVALENCE (WYS(1),PHT(1,1,1,1))
C
      NXYZ1 = NX1*NY1*NZ1
      NTOT1 = NXYZ1*NEL
      CALL COL3 (PHI,SFV,WY,NTOT1)
C
      DO 100 IEL=1,NEL
      IF ( IFRZER(IEL) ) THEN
         IF (IDIR.EQ.1) THEN
            CALL MXM (WY(1,1,1,IEL),NX1,DATM1,NY1,WYS,1)
            DO 220 IX=1,NX1
            FMS(IX,1,1,IEL)= FMS(IX,1,1,IEL) + WXM1(IX)*WAM1(1)*
     $                       WYS(IX)*SFV(IX,1,1,IEL)*JACM1(IX,1,1,IEL)
  220       CONTINUE
         ENDIF
         DO 320 IX=1,NX1
         DO 320 IY=2,NY1
            FMS(IX,IY,1,IEL)=FMS(IX,IY,1,IEL) + PHI(IX,IY,1,IEL) *
     $                       BM1(IX,IY,1,IEL) / YM1(IX,IY,1,IEL)
  320    CONTINUE
      ELSE
         CALL ADDCOL4 (FMS(1,1,1,IEL),PHI(1,1,1,IEL),JACM1(1,1,1,IEL),
     $                 W2CM1,NXYZ1)
      ENDIF
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine updcoor
C-----------------------------------------------------------------------
C
C     Subroutine to update geometry for moving boundary problems
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TSTEP'
C
      IFIELD = 0
      NEL    = NELFLD(IFIELD)
C
C     Update collocation points coordinates
C
      CALL UPDXYZ (NEL)
C
C     Shift lagged mesh velocity
C
      CALL LAGMSHV (NEL)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine mvbdry (nel)
C
C     Routine to evaluate mesh velocities at all moving boundaries
C
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MVGEOM'
      include 'SOLN'
      include 'TSTEP'
      COMMON /SCRSF/ WVX(LX1,LY1,LZ1,LELT)
     $             , WVY(LX1,LY1,LZ1,LELT)
     $             , WVZ(LX1,LY1,LZ1,LELT)
      COMMON /SCRCH/ WTX(LX1,LY1,LZ1,LELT)
     $             , WTY(LX1,LY1,LZ1,LELT)
      COMMON /SCRMG/ WTZ(LX1,LY1,LZ1,LELT)
     $             , RNX(LX1,LY1,LZ1,LELT)
     $             , RNY(LX1,LY1,LZ1,LELT)
     $             , RNZ(LX1,LY1,LZ1,LELT)
      COMMON /SCRUZ/ DSA(LX1,LY1,LZ1,LELT)
     $             , QNI(LX1,LY1,LZ1,LELT)
     $             , SMT(LX1,LY1,LZ1,LELT)
     $             , TA (LX1,LY1,LZ1,LELT)
C
      LOGICAL IFALGN,IFNORX,IFNORY,IFNORZ,IFDSMV,IFREGW
      CHARACTER CB*3
C
      IFIELD = 0
      NXYZ1  = NX1*NY1*NZ1
      NTOT1  = NX1*NY1*NZ1*NEL
      NFACE  = 2*NDIM
      CALL RZERO3  (RNX,RNY,RNZ,NTOT1)
C
      DO 100 IEL=1,NEL
      DO 100 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MS ' .OR. CB.EQ.'ms ' .OR. 
     $       CB.EQ.'MSI' .OR. CB.EQ.'msi' .OR. 
     $       CB.EQ.'MM ' .OR. CB.EQ.'mm ' .OR. 
     $       CB.EQ.'mv ' .OR. CB.EQ.'mvn' .OR.
     $       CB.EQ.'MLI') THEN
             CALL FACEXV (UNX(1,1,IFC,IEL),UNY(1,1,IFC,IEL),
     $                    UNZ(1,1,IFC,IEL),RNX(1,1,1,IEL),
     $                    RNY(1,1,1,IEL),RNZ(1,1,1,IEL),IFC,1)
         ENDIF
  100 CONTINUE
C
      CALL DSSUM (RNX,NX1,NY1,NZ1)
      CALL DSSUM (RNY,NX1,NY1,NZ1)
      IF (NDIM.EQ.3) CALL DSSUM (RNZ,NX1,NY1,NZ1)
      CALL UNITVEC (RNX,RNY,RNZ,NTOT1)
C
      CALL RZERO3 (WVX,WVY,WVZ,NTOT1)
      CALL RZERO3 (WTX,WTY,WTZ,NTOT1)
      DO 1000 ISWEEP=1,2
C
      IFREGW = .FALSE.
      IFDSMV = .FALSE.
      CALL RZERO  (DSA,NTOT1)
      CALL RZERO  (TA,NTOT1)
C
      IF (IFFLOW) THEN
C
      IFIELD = 1
      CALL RZERO  (SMT,NTOT1)
      DO 210 IEL=1,NELV
      DO 210 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'mv ' .OR. CB.EQ.'mvn'  .OR.
     $       CB.EQ.'MM ' .OR. CB.EQ.'mm '  .OR.
     $       CB.EQ.'MS ' .OR. CB.EQ.'ms ') THEN
            IFREGW = .TRUE.
            CALL FACEC3 (WVX(1,1,1,IEL),WVY(1,1,1,IEL),WVZ(1,1,1,IEL),
     $                   VX(1,1,1,IEL),VY(1,1,1,IEL),VZ(1,1,1,IEL),IFC)
            IF (CB.NE.'mv ')
     $      CALL NORCMP (WVX(1,1,1,IEL),WVY(1,1,1,IEL),WVZ(1,1,1,IEL),
     $                   RNX(1,1,1,IEL),RNY(1,1,1,IEL),RNZ(1,1,1,IEL),
     $                 IFC) 
         ENDIF
         IF (CB.EQ.'MSI' .OR. CB.EQ.'msi') THEN
            IFDSMV = .TRUE.
            CALL FACSMT (SMT(1,1,1,IEL),IFC)
         ENDIF
  210 CONTINUE

      if (istep.eq.0) call opcopy(wx,wy,wz,wvx,wvy,wvz)
c     call outpost(wvx,wvy,wvz,wtx,wtx,'   ')
c     call outpost(wx,wy,wz,wtx,wtx,'   ')
c     write(6,*) 'quit1'
c     call exitt

      iregw = 0                   ! Global handshake on logicals ifregw, ifdsmv
      if (ifregw) iregw = 1       ! pff  9/11/07
      iregw = iglmax(iregw,1)
      if (iregw.eq.1) ifregw = .true.

      idsmv = 0
      if (ifdsmv) idsmv = 1
      idsmv = iglmax(idsmv,1)
      if (idsmv.eq.1) ifdsmv = .true.

C
      IF (IFDSMV) THEN
      CALL DSSUM (SMT,NX1,NY1,NZ1)
      DO 215 IEL=1,NELV
      DO 215 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MSI' .OR. CB.EQ.'msi') THEN
         CALL FACEC3 (WTX(1,1,1,IEL),WTY(1,1,1,IEL),WTZ(1,1,1,IEL),
     $                VX(1,1,1,IEL),VY(1,1,1,IEL),VZ(1,1,1,IEL),IFC)
         CALL FACEMV (WTX(1,1,1,IEL),WTY(1,1,1,IEL),WTZ(1,1,1,IEL),
     $                RNX(1,1,1,IEL),RNY(1,1,1,IEL),RNZ(1,1,1,IEL),
     $                SMT(1,1,1,IEL),IFC) 
         ENDIF
  215 CONTINUE
         CALL DSSUM (WTX,NX1,NY1,NZ1)
         CALL DSSUM (WTY,NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSSUM (WTZ,NX1,NY1,NZ1)
      ENDIF
C
      ENDIF
C
      IF (IFMELT .AND. ISTEP.GT.0) THEN
         IFIELD = 2
         CALL RZERO (SMT,NTOT1)
         CALL CQNET (QNI,TA,NEL)
      DO 220 IEL=1,NELT
      DO 220 IFC=1,NFACE
         CB   = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MLI') THEN
            CALL FACSMT (SMT(1,1,1,IEL),IFC)
         ENDIF
         IF (CB.EQ.'MLI' .OR. CB.EQ.'MCI') THEN
            CALL FACEXS (AREA(1,1,IFC,IEL),TA,IFC,1)
            CALL ADD2   (DSA(1,1,1,IEL),TA,NXYZ1)
         ENDIF
  220 CONTINUE
         CALL DSSUM (SMT,NX1,NY1,NZ1)
         CALL DSSUM (DSA,NX1,NY1,NZ1)
      DO 280 IEL=1,NELT
      DO 280 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MLI') THEN
            RHOLA = -0.5 * BC(5,IFC,IEL,IFIELD)
            CALL FACEMT (WTX(1,1,1,IEL),WTY(1,1,1,IEL),WTZ(1,1,1,IEL),
     $                   RNX(1,1,1,IEL),RNY(1,1,1,IEL),RNZ(1,1,1,IEL),
     $                   QNI(1,1,1,IEL),DSA(1,1,1,IEL),SMT(1,1,1,IEL),
     $                   RHOLA,IFC)
         ENDIF
  280 CONTINUE
         CALL DSSUM (WTX,NX1,NY1,NZ1)
         CALL DSSUM (WTY,NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSSUM (WTZ,NX1,NY1,NZ1)
      ENDIF
C
      IFIELD = 0
      DO 330 IEL=1,NEL
      DO 330 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'SYM') THEN
            CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
            IF (IFREGW) THEN
               IF (IFNORX) CALL FACEV (WVX,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORY) CALL FACEV (WVY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORZ) CALL FACEV (WVZ,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (.NOT.IFALGN) 
     $         CALL FACZQN (WVX(1,1,1,IEL),WVY(1,1,1,IEL),
     $                      WVZ(1,1,1,IEL),IFC,IEL)
            ENDIF
            IF (IFDSMV .OR. IFMELT) THEN
               IF (IFNORX) CALL FACEV (WTX,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORY) CALL FACEV (WTY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORZ) CALL FACEV (WTZ,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (.NOT.IFALGN) 
     $         CALL FACZQN (WTX(1,1,1,IEL),WTY(1,1,1,IEL),
     $                      WTZ(1,1,1,IEL),IFC,IEL)
            ENDIF
         ENDIF
  330 CONTINUE
C
      DO 350 IEL=1,NEL
      DO 350 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'FIX') THEN
            IF (IFREGW) THEN
               CALL FACEV (WVX,IEL,IFC,0.0,NX1,NY1,NZ1)
               CALL FACEV (WVY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (NDIM.EQ.3) CALL FACEV (WVZ,IEL,IFC,0.0,NX1,NY1,NZ1)
            ENDIF
            IF (IFDSMV .OR. IFMELT) THEN
               CALL FACEV (WTX,IEL,IFC,0.0,NX1,NY1,NZ1)
               CALL FACEV (WTY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (NDIM.EQ.3) CALL FACEV (WTZ,IEL,IFC,0.0,NX1,NY1,NZ1)
            ENDIF
         ENDIF
  350 CONTINUE 
C
      IF (ISWEEP.EQ.1) THEN
         IF (IFREGW) THEN
            CALL DSOP (WVX,'MXA',NX1,NY1,NZ1)
            CALL DSOP (WVY,'MXA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WVZ,'MXA',NX1,NY1,NZ1)
         ENDIF
         IF (IFDSMV .OR. IFMELT) THEN
            CALL DSOP (WTX,'MXA',NX1,NY1,NZ1)
            CALL DSOP (WTY,'MXA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WTZ,'MXA',NX1,NY1,NZ1)
         ENDIF
      ELSE
         IF (IFREGW) THEN
            CALL DSOP (WVX,'MNA',NX1,NY1,NZ1)
            CALL DSOP (WVY,'MNA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WVZ,'MNA',NX1,NY1,NZ1)
         ENDIF
         IF (IFDSMV .OR. IFMELT) THEN
            CALL DSOP (WTX,'MNA',NX1,NY1,NZ1)
            CALL DSOP (WTY,'MNA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WTZ,'MNA',NX1,NY1,NZ1)
         ENDIF
      ENDIF
C
 1000 CONTINUE
C
      CALL RMASK (WX,WY,WZ,NEL)
c     call outpost(wx,wy,wz,wtx,wtx,'   ')
c     write(6,*) 'quit2'
c     call exitt
C
      IF (IFREGW) THEN
         CALL ADD2  (WX,WVX,NTOT1)
         CALL ADD2  (WY,WVY,NTOT1)
         IF (NDIM.EQ.3) CALL ADD2  (WZ,WVZ,NTOT1)
      ENDIF
      IF (IFDSMV .OR. IFMELT) THEN
         CALL ADD2  (WX,WTX,NTOT1)
         CALL ADD2  (WY,WTY,NTOT1)
         IF (NDIM.EQ.3) CALL ADD2  (WZ,WTZ,NTOT1)
      ENDIF
 
c     call outpost(wx,wy,wz,wtx,wtx,'   ')
c     call exitt
 
      return
      end
c-----------------------------------------------------------------------
      subroutine norcmp (wt1,wt2,wt3,rnx,rny,rnz,ifc)
C
      include 'SIZE'
      COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)
C
      DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1)
     $        , RNX(LX1,LY1,LZ1),RNY(LX1,LY1,LZ1),RNZ(LX1,LY1,LZ1)
C
      NXYZ1 = NX1*NY1*NZ1
C
      CALL COPY (R1,WT1,NXYZ1)
      CALL COPY (R2,WT2,NXYZ1)
      IF (NDIM.EQ.3) CALL COPY (R3,WT3,NXYZ1)
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C
      IF (NDIM.EQ.2) THEN
         DO 200 J2=JS2,JF2,JSKIP2
         DO 200 J1=JS1,JF1,JSKIP1
            WN          = R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                    R2(J1,J2,1)*RNY(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
  200    CONTINUE
      ELSE
         DO 300 J2=JS2,JF2,JSKIP2
         DO 300 J1=JS1,JF1,JSKIP1
            WN          = R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                    R2(J1,J2,1)*RNY(J1,J2,1) +
     $                    R3(J1,J2,1)*RNZ(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
            WT3(J1,J2,1) = WN *RNZ(J1,J2,1)
  300    CONTINUE
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine facemv (wt1,wt2,wt3,rnx,rny,rnz,smt,ifc)
C
      include 'SIZE'
      COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)
C
      DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1)
     $        , RNX(LX1,LY1,LZ1),RNY(LX1,LY1,LZ1),RNZ(LX1,LY1,LZ1)
     $        , SMT(LX1,LY1,LZ1)
C
      NXYZ1 = NX1*NY1*NZ1
C
      CALL COPY (R1,WT1,NXYZ1)
      CALL COPY (R2,WT2,NXYZ1)
      IF (NDIM.EQ.3) CALL COPY (R3,WT3,NXYZ1)
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C
      IF (NDIM.EQ.2) THEN
         DO 200 J2=JS2,JF2,JSKIP2
         DO 200 J1=JS1,JF1,JSKIP1
            WN          = ( R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                      R2(J1,J2,1)*RNY(J1,J2,1) ) / SMT(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
  200    CONTINUE
      ELSE
         DO 300 J2=JS2,JF2,JSKIP2
         DO 300 J1=JS1,JF1,JSKIP1
            WN          = ( R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                      R2(J1,J2,1)*RNY(J1,J2,1) +
     $                      R3(J1,J2,1)*RNZ(J1,J2,1) ) / SMT(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
            WT3(J1,J2,1) = WN *RNZ(J1,J2,1)
  300    CONTINUE
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine faczqn (wt1,wt2,wt3,ifc,iel)
C
      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'
      COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)
C
      DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1)
C
      NXYZ1 = NX1*NY1*NZ1
      CALL COPY (R1,WT1,NXYZ1)
      CALL COPY (R2,WT2,NXYZ1)
      IF (NDIM.EQ.3) CALL COPY (R3,WT3,NXYZ1)
C
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
      I = 0
C
      IF (NDIM.EQ.2) THEN
      DO 200 J2=JS2,JF2,JSKIP2
      DO 200 J1=JS1,JF1,JSKIP1
         I = I+1
         W1           = R1(J1,J2,1)*T1X(I,1,IFC,IEL) + 
     $                  R2(J1,J2,1)*T1Y(I,1,IFC,IEL)
         WT1(J1,J2,1) = W1 *T1X(I,1,IFC,IEL)
         WT2(J1,J2,1) = W1 *T1Y(I,1,IFC,IEL)
  200 CONTINUE
      ELSE
      DO 300 J2=JS2,JF2,JSKIP2
      DO 300 J1=JS1,JF1,JSKIP1
         I = I+1
         W1           = R1(J1,J2,1)*T1X(I,1,IFC,IEL) + 
     $                  R2(J1,J2,1)*T1Y(I,1,IFC,IEL) +
     $                  R3(J1,J2,1)*T1Z(I,1,IFC,IEL)
         WT1(J1,J2,1) = W1 *T1X(I,1,IFC,IEL)
         WT2(J1,J2,1) = W1 *T1Y(I,1,IFC,IEL)
         WT3(J1,J2,1) = W1 *T1Z(I,1,IFC,IEL)
  300 CONTINUE
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine facsmt (smt,ifc)
C
      include 'SIZE'
      DIMENSION SMT(LX1,LY1,LZ1)
C
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         SMT(J1,J2,1)=SMT(J1,J2,1) + 1.0
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine cqnet (qni,ta,nel)
C
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT)
     $             , H2(LX1,LY1,LZ1,LELT)
C
      DIMENSION QNI(LX1,LY1,LZ1,1)
     $        , TA (LX1,LY1,LZ1,1)
C
      INTLOC = -1
      IMSHL  =  2
      NTOT1  = NX1*NY1*NZ1*NEL
C
      CALL SETHLM (H1,H2,INTLOC)
      CALL AXHELM (TA,T(1,1,1,1,IFIELD-1),H1,H2,IMSHL,1)
      CALL SUB3   (QNI,TA,BQ(1,1,1,1,IFIELD-1),NTOT1)
      CALL DSSUM  (QNI,NX1,NY1,NZ1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine facemt (w1,w2,w3,rnx,rny,rnz,qni,dsa,smt,rhola,ifc)
C
      include 'SIZE'
      include 'GEOM'
C
      DIMENSION  W1 (LX1,LY1,LZ1)
     $        ,  W2 (LX1,LY1,LZ1)
     $        ,  W3 (LX1,LY1,LZ1)
     $        ,  RNX(LX1,LY1,LZ1)
     $        ,  RNY(LX1,LY1,LZ1)
     $        ,  RNZ(LX1,LY1,LZ1)
     $        ,  QNI(LX1,LY1,LZ1)
     $        ,  DSA(LX1,LY1,LZ1)
     $        ,  SMT(LX1,LY1,LZ1)
C
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C
      IF (NDIM.EQ.2) THEN
      DO 200 J2=JS2,JF2,JSKIP2
      DO 200 J1=JS1,JF1,JSKIP1
         AA = QNI(J1,J2,1) / ( DSA(J1,J2,1)*SMT(J1,J2,1)*RHOLA )
         W1(J1,J2,1) = RNX(J1,J2,1) * AA
         W2(J1,J2,1) = RNY(J1,J2,1) * AA
  200 CONTINUE
      ELSE
      DO 300 J2=JS2,JF2,JSKIP2
      DO 300 J1=JS1,JF1,JSKIP1
         AA = QNI(J1,J2,1) / ( DSA(J1,J2,1)*SMT(J1,J2,1)*RHOLA )
         W1(J1,J2,1) = RNX(J1,J2,1) * AA
         W2(J1,J2,1) = RNY(J1,J2,1) * AA
         W3(J1,J2,1) = RNZ(J1,J2,1) * AA
  300 CONTINUE
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine elasolv (nel)
C
C     Elastostatic solver for mesh deformation
C
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MVGEOM'
      include 'SOLN'
      include 'TSTEP'
C
      COMMON /SCRNS/ DW1  (LX1,LY1,LZ1,LELT)
     $             , DW2  (LX1,LY1,LZ1,LELT)
     $             , DW3  (LX1,LY1,LZ1,LELT)
     $             , AW1  (LX1,LY1,LZ1,LELT)
     $             , AW2  (LX1,LY1,LZ1,LELT)
     $             , AW3  (LX1,LY1,LZ1,LELT)
      COMMON /SCRVH/ H1   (LX1,LY1,LZ1,LELT)
     $             , H2   (LX1,LY1,LZ1,LELT)
      common /scruz/ prt  (lx1,ly1,lz1,lelt)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
C     Set up parameters for mesh solver
C
c     return
c     if (istep.gt.1) then
c        ifldx = ifield
c        ifield = 1
c        call incomprn (wx,wy,wz,prt) ! project U onto div-free space
c        ifield = ifldx
c        return
c     endif
c
c     call quickmv

      if (ifusermv) return  ! Compute wx,wy,wz in userchk.
c
      NTOT1  = NX1*NY1*NZ1*NEL
      MAXIT  = 1000
      MATMOD = -1
      IFH2   = .FALSE.
      IFSOLV = .TRUE.
      IMSOLV = 0
      VNU    = 0.0
      VNU    = 0.4
      VNU    = param(47)
      vnu    = max(0.00,vnu)
      vnu    = min(0.49,vnu)
C
C     Set up elastic material constants
C
      CE = 1./(1. + VNU)
      C2 = VNU * CE / (1. - 2.*VNU)
      C3 = 0.5 * CE
      CALL CFILL (H1,C2,NTOT1)
      CALL CFILL (H2,C3,NTOT1)
C
C     Solve for interior mesh velocities
C
      CALL MESHTOL (AW1,TOLMSH,NEL,IMSOLV)
      IF (IMSOLV.EQ.1) return
C
      CALL AXHMSF  (AW1,AW2,AW3,WX,WY,WZ,H1,H2,MATMOD)
c     call outpost(wx,wy,wz,h1,h2,'   ')
c     call exitt
      CALL CHSIGN  (AW1,NTOT1)
      CALL CHSIGN  (AW2,NTOT1)
      IF (NDIM.EQ.3) CALL CHSIGN (AW3,NTOT1)
      CALL HMHZSF  ('NOMG',DW1,DW2,DW3,AW1,AW2,AW3,H1,H2,
     $               W1MASK,W2MASK,W3MASK,WMULT,TOLMSH,
     $               MAXIT,MATMOD)
C
C     Update mesh velocities
C
      CALL ADD2 (WX,DW1,NTOT1)
      CALL ADD2 (WY,DW2,NTOT1)
      IF (NDIM.EQ.3) CALL ADD2 (WZ,DW3,NTOT1)

c     if (istep.gt.1) then
c        ifldx = ifield
c        ifield = 1
c        call incomprn (wx,wy,wz,prt) ! project U onto div-free space
c        ifield = ifldx
c        return
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine meshtol (ta,tolmsh,nel,imsolv)
C
      include 'SIZE'
      include 'EIGEN'
      include 'MVGEOM'
      include 'TSTEP'
      DIMENSION TA(LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NEL
      TOLAB = TOLREL
C
      DELTA  = 1.0E-9
      X      = 1.0 + DELTA
      Y      = 1.0
      DIFF   = ABS(X - Y)
      IF (DIFF .EQ. 0.0) EPS = 1.0E-05
      IF (DIFF .GT. 0.0) EPS = 1.0E-12
C
      CALL OPDOT  (TA,WX,WY,WZ,WX,WY,WZ,NTOT1)
C
      WDOT = GLMAX(TA,NTOT1)
      WMAX = SQRT(WDOT)
      IF (WMAX .LT. EPS) THEN
          IMSOLV = 1
          return
      ENDIF
C
      TOLMSH = TOLAB * WMAX * SQRT(EIGAA)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine updxyz (nel)
C
      include 'SIZE'
      include 'TSTEP'
      include 'MVGEOM'
      include 'GEOM'
      COMMON /SCRSF/ UX(LX1,LY1,LZ1,LELT)
     $             , UY(LX1,LY1,LZ1,LELT)
     $             , UZ(LX1,LY1,LZ1,LELT)
      DIMENSION ABM(3)
C
      NTOT1 = NX1*NY1*NZ1*NEL
C
      DO 10 I=1,NBD
   10 ABM(I) = DT*ABMSH(I)
C
      IF (ISTEP.EQ.0) THEN
         CALL COPY (UX,WX,NTOT1)
         CALL COPY (UY,WY,NTOT1)
         IF (NDIM.EQ.3) CALL COPY (UZ,WZ,NTOT1)
      ELSE
         CALL CMULT2 (UX,WX,ABM(1),NTOT1)
         CALL CMULT2 (UY,WY,ABM(1),NTOT1)
         IF (NDIM.EQ.3) CALL CMULT2 (UZ,WZ,ABM(1),NTOT1)
         DO 100 ILAG=2,NBD
            CALL ADD2S2 (UX,WXLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
            CALL ADD2S2 (UY,WYLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
            IF (NDIM.EQ.3) 
     $      CALL ADD2S2 (UZ,WZLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
  100    CONTINUE
      ENDIF
C
      CALL ADD2 (XM1,UX,NTOT1)
      CALL ADD2 (YM1,UY,NTOT1)
      IF (NDIM.EQ.3) CALL ADD2 (ZM1,UZ,NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine lagmshv (nel)
C-----------------------------------------------------------------------
C
C     Keep old mesh velocity
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MVGEOM'
      include 'TSTEP'
C
      NTOT1 = NX1*NY1*NZ1*NEL
C
      DO 100 ILAG=NBDINP-1,2,-1
         CALL COPY (WXLAG(1,1,1,1,ILAG),WXLAG(1,1,1,1,ILAG-1),NTOT1)
         CALL COPY (WYLAG(1,1,1,1,ILAG),WYLAG(1,1,1,1,ILAG-1),NTOT1)
         IF (NDIM.EQ.3)
     $   CALL COPY (WZLAG(1,1,1,1,ILAG),WZLAG(1,1,1,1,ILAG-1),NTOT1)
 100  CONTINUE
C
         CALL COPY (WXLAG(1,1,1,1,1),WX,NTOT1)
         CALL COPY (WYLAG(1,1,1,1,1),WY,NTOT1)
         IF (NDIM.EQ.3)
     $   CALL COPY (WZLAG(1,1,1,1,1),WZ,NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine facec3 (a1,a2,a3,b1,b2,b3,ifc)
C
C     Copy the face (IFC) of B1,B2,B3 to A1,A2,A3.
C     IFACE is the input in the pre-processor ordering scheme.
C
      include 'SIZE'
      DIMENSION A1(LX1,LY1,LZ1)
     $        , A2(LX1,LY1,LZ1)
     $        , A3(LX1,LY1,LZ1)
     $        , B1(LX1,LY1,LZ1)
     $        , B2(LX1,LY1,LZ1)
     $        , B3(LX1,LY1,LZ1)
C
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         A1(J1,J2,1)=B1(J1,J2,1)
         A2(J1,J2,1)=B2(J1,J2,1)
         A3(J1,J2,1)=B3(J1,J2,1)
  100 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine ptbgeom
C-----------------------------------------------------------------------
C
C     Subroutine to impose perturbation to geometry before solution
C     for moving boundary problems
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'TSTEP'
      include 'INPUT'
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C
      IF (ISTEP .EQ. 0) return
      IFIELD = 0
      NEL    = NELFLD(IFIELD)
      NTOT1  = NX1*NY1*NZ1*NEL
      IMESH  = 1
      IF ( IFTMSH(IFIELD) ) IMESH = 2
C
      CALL IBDGEOM (NEL)
      CALL ELASOLV (NEL)
      CALL UPDXYZ  (NEL)
      CALL GEOM1 (XM3,YM3,ZM3)
      CALL GEOM2
      CALL UPDMSYS (0)
      CALL VOLUME
      CALL SETINVM
      CALL LAGMASS
C
      CALL RZERO (WX,NTOT1)
      CALL RZERO (WY,NTOT1)
      IF (NDIM.EQ.3) CALL RZERO (WZ,NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine ibdgeom (nel)
C
C     Routine to evaluate mesh velocities at all moving boundaries
C
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'MVGEOM'
      include 'TSTEP'
C
      CHARACTER CB*1
C
      NFACE  = 2*NDIM
      NTOT1  = NX1*NY1*NZ1*NEL
C
      CALL RZERO (WX,NTOT1)
      CALL RZERO (WY,NTOT1)
      CALL RZERO (WZ,NTOT1)
C
      DO 1000 ISWEEP=1,2
C
      IFLD = 0
      DO 110 IEL=1,NEL
      DO 110 IFC=1,NFACE
         ieg = lglel(iel)
         CB  = CBC(IFC,IEL,IFLD)
      IF (CB.EQ.'M' .OR. CB.EQ.'m') THEN
         CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFC)
         DO 140 IZ=KZ1,KZ2
         DO 140 IY=KY1,KY2
         DO 140 IX=KX1,KX2
         CALL INIGEOM (WX(IX,IY,IZ,IEL),WY(IX,IY,IZ,IEL),
     $                 WZ(IX,IY,IZ,IEL),XM1(IX,IY,IZ,IEL),
     $                 YM1(IX,IY,IZ,IEL),ZM1(IX,IY,IZ,IEL),
     $                 IFC,IEG)
  140    CONTINUE
      ENDIF
  110 CONTINUE
C
      IF (ISWEEP.EQ.1) THEN
         CALL DSOP (WX,'MXA',NX1,NY1,NZ1)
         CALL DSOP (WY,'MXA',NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSOP (WZ,'MXA',NX1,NY1,NZ1)
      ELSE
         CALL DSOP (WX,'MNA',NX1,NY1,NZ1)
         CALL DSOP (WY,'MNA',NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSOP (WZ,'MNA',NX1,NY1,NZ1)
      ENDIF
C
 1000 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine inigeom (ux,uy,uz,x,y,z,iside,iel)
C
      include 'SIZE'
      include 'TSTEP'
C
      UX  =  0.0
      UY  =  0.0
      UZ  =  0.0
C
      return
      end
c-----------------------------------------------------------------------
      subroutine quickmv
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
c
      if (if3d) then
         call quickmv3d
      else
         call quickmv2d
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine quickmv2d
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
c
      integer e,ex,ey,ez,eg
      common /surfa/ zsurf(lx1,lz1,lelx,lely)
     $             , wsurf(lx1,lz1,lelx,lely)
      real nxs,nys,nzs
c
      icount = 0
      do ex=1,nelx
      do ix=1,nx1
         zsurf(ix,1,ex,1) = -1.e20
         wsurf(ix,1,ex,1) = -1.e20
         ey=nely
         eg  = ex + nelx*(ey-1) 
         mid = gllnid(eg)
         e   = gllel (eg)
         if (mid.eq.nid) then
            zsurf(ix,1,ex,1) = ym1(ix,ny1,1,e)
            vxs              = vx (ix,ny1,1,e)
            vys              = vy (ix,ny1,1,e)
            nxs              = unx(ix,1,3,e)          ! Face 3 is on top in 2D
            nys              = uny(ix,1,3,e)
            gamma_s          = (vxs*nxs + vys*nys)/(nys) 
            wsurf(ix,1,ex,1) = gamma_s           ! vertical component of wsurf
         endif
         zsurf(ix,1,ex,1) = glmax(zsurf(ix,1,ex,1),1)
         wsurf(ix,1,ex,1) = glmax(wsurf(ix,1,ex,1),1)
         icount = icount+1

c        write(6,6) ex,e,ix,xm1(ix,ny1,1,e),ym1(ix,ny1,1,e)
c    $   ,vxs,vys,nxs,nys,gamma_s,wsurf(ix,1,ex,1),zsurf(ix,1,ex,1)
c   6        format(3i3,1p9e12.4,' srf')

      enddo
      enddo
      zmin = glmin(ym1,nx1*ny1*nz1*nelv)
c
      do ex=1,nelx
      do ix=1,nx1
         do ey=1,nely
            eg  = ex + nelx*(ey-1) 
            mid = gllnid(eg)
            e   = gllel (eg)
            if (mid.eq.nid) then
               do iy=1,ny1
                  wy (ix,iy,1,e) = wsurf(ix,1,ex,1)
     $               * (ym1(ix,iy,1,e)-zmin)/(zsurf(ix,1,ex,1)-zmin)
               enddo
            endif
         enddo
      enddo
      enddo
c
      n = nx1*ny1*nz1*nelv
      call rzero(wx,n)
      call dsavg(wy)

c     call opcopy(vx,vy,vz,wx,wy,wz)
c     call outpost(vx,vy,vz,pr,t,'   ')
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine quickmv3d
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
c
      integer e,ex,ey,ez,eg
      common /surfa/ zsurf(lx1,lz1,lelx,lely)
     $             , wsurf(lx1,lz1,lelx,lely)
      real nxs,nys,nzs
c
      icount = 0
      do ey=1,nely
      do ex=1,nelx
      do iy=1,ny1
      do ix=1,nx1
         zsurf(ix,iy,ex,ey) = -1.e20
         wsurf(ix,iy,ex,ey) = -1.e20
         ez  = nelz
         eg  = ex + nelx*(ey-1) + nelx*nely*(ez-1)
         mid = gllnid(eg)
         e   = gllel (eg)
         if (mid.eq.nid) then
            zsurf(ix,iy,ex,ey) = zm1(ix,iy,nz1,e)
            vxs                = vx (ix,iy,nz1,e)
            vys                = vy (ix,iy,nz1,e)
            vzs                = vz (ix,iy,nz1,e)
            nxs                = unx(ix,iy,6,e)     ! Face 6 is on top in 3D
            nys                = uny(ix,iy,6,e)
            nzs                = unz(ix,iy,6,e) 
            gamma_s            = (vxs*nxs+vys*nys+vzs*nzs)/(nzs)
            wsurf(ix,iy,ex,ey) = gamma_s  ! vertical component of wsurf
         endif
         zsurf(ix,iy,ex,ey) = glmax(zsurf(ix,iy,ex,ey),1)
         wsurf(ix,iy,ex,ey) = glmax(wsurf(ix,iy,ex,ey),1)
         icount = icount+1
      enddo
      enddo
      enddo
      enddo

      n = nx1*ny1*nz1*nelv
      zmin = glmin(zm1,n)
c
      do ey=1,nely
      do ex=1,nelx
      do iy=1,ny1
      do ix=1,nx1
         do ez=1,nelz
            eg  = ex + nelx*(ey-1) + nelx*nely*(ez-1)
            mid = gllnid(eg)
            e   = gllel (eg)
            if (mid.eq.nid) then
               do iz=1,nz1
                  wz (ix,iy,iz,e) = wsurf(ix,iy,ex,ey)
     $               * (zm1(ix,iy,iz,e)-zmin)/(zsurf(ix,iy,ex,ey)-zmin)
               enddo
            endif
         enddo
      enddo
      enddo
      enddo
      enddo
c
      n = nx1*ny1*nz1*nelv
      call rzero(wx,n)
      call rzero(wy,n)
      call dsavg(wz)
c
      return
      end
c-----------------------------------------------------------------------
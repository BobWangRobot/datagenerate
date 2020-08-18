      SUBROUTINE RANMAR(RVEC,LEN)
c modlib.f (ccp4 cvs 2010-05-25, but has been stable for a long time)
C     ===========================
C
C     Universal random number generator proposed by Marsaglia and Zaman
C     in report FSU-SCRI-87-50
C     slightly modified by F. James, 1988 to generate a vector
C     of pseudorandom numbers RVEC of length LEN
C     and making the COMMON block include everything needed to
C     specify completely the state of the generator.
C     Transcribed from CERN report DD/88/22.
C     Rather inelegant messing about added by D. Love, Jan. 1989 to
C     make sure initialisation always occurs.
C     *** James says that this is the preferred generator.
C     Gives bit-identical results on all machines with at least
C     24-bit mantissas in the flotaing point representation (i.e.
C     all common 32-bit computers. Fairly fast, satisfies very
C     stringent tests, has very long period and makes it very
C     simple to generate independly disjoint sequences.
C     See also RANECU.
C     The state of the generator may be saved/restored using the
C     whole contents of /RASET1/.
C     Call RANMAR to get a vector, RMARIN to initialise.
C
C  Argument list
C  -------------
C
C     VREC (O)                 (REAL)   Random Vector
C
C     LEN  (I)              (INTEGER)   Length of random vector
C
C
C  For ENTRY point RMARIN
C  ----------------------
C
C     Initialisation for RANMAR.  The input values should
C     be in the ranges: 0<=ij<=31328, 0<=kl<=30081
C     This shows the correspondence between the simplified input seeds
C     IJ, KL and the original Marsaglia-Zaman seeds i,j,k,l
C     To get standard values in Marsaglia-Zaman paper,
C     (I=12, J=34, K=56, L=78) put IJ=1802, KL=9373
C
C     IJ   (I)              (INTEGER)   Seed for random number generator
C
C     KL   (I)              (INTEGER)   Seed for randon number generator
C
C_END_RANMAR
C
C     ..
C     .. Agruments ..
      REAL RVEC(*)
      INTEGER LEN,IJ,KL
C     ..
C     .. Common Variables ..
      REAL C,CD,CM,U
      INTEGER I97,J97
C     ..
C     .. Local Scalars ..
      REAL S,T,UNI
      INTEGER I,II,IVEC,J,JJ,K,L,M
      LOGICAL INITED
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common Blocks ..
      COMMON /RASET1/ U(97),C,CD,CM,I97,J97
C     ..
C     .. Save Statement ..
      SAVE INITED, /RASET1/
C     ..
C     .. Data Statement ..
      DATA INITED /.FALSE./
C
C---- If initialised, fill RVEC and RETURN. If not, do initialisation
C     and return here later.
C
 1    IF (INITED) THEN
        DO 100 IVEC=1,LEN
          UNI=U(I97)-U(J97)
          IF (UNI.LT.0.) UNI=UNI+1.
          U(I97)=UNI
          I97=I97-1
          IF (I97.EQ.0) I97=97
          J97=J97-1
          IF (J97.EQ.0) J97=97
          C=C-CD
          IF (C.LT.0.) C=C+CM
          UNI=UNI-C
          IF (UNI.LT.0.) UNI=UNI+1.
          RVEC(IVEC)=UNI
 100    CONTINUE
        RETURN
      ENDIF
      I=MOD(1802/177,177)+2
      J=MOD(1802,177)+2
      K=MOD(9373/169,178)+1
      L=MOD(9373,169)
C
      GOTO 10
C
C---- Initialise and return without filling RVEC
C
      ENTRY RMARIN(IJ,KL)
      I=MOD(IJ/177,177)+2
      J=MOD(IJ,177)+2
      K=MOD(KL/169,178)+1
      L=MOD(KL,169)
      INITED=.TRUE.

 10   CONTINUE
      DO 2 II=1,97
        S=0.
        T=.5
        DO 3 JJ=1,24
          M=MOD(MOD(I*J,179)*K,179)
          I=J
          J=K
          K=M
          L=MOD(53*L+1,169)
          IF (MOD(L*M,64).GE.32) S=S+T
          T=0.5*T
 3      CONTINUE
        U(II)=S
 2    CONTINUE
      C=362436./16777216.
      CD=7654321./16777216.
      CM=16777213./16777216.
      I97=97
      J97=33
      IF (.NOT. INITED) THEN
        INITED=.TRUE.
        GOTO 1
      ENDIF

      END

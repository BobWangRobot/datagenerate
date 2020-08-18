c ccplib.f
C     ==================================
      REAL FUNCTION STHLSQ(IH, IK, IL)
C     ==================================
C     
      IMPLICIT NONE
C     
      EXTERNAL STHLSQ1
C
C     PARAMETERS
      INTEGER IH, IK, IL
C
C     reso RETURN value
      REAL RESO
C
      CALL STHLSQ1(RESO, IH, IK, IL)
      STHLSQ=RESO
      RETURN
      END

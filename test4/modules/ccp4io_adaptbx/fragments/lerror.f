c parser.f
C     =======================================
      SUBROUTINE LERROR(ERRFLG,IFAIL,ERRMSG)
C     =======================================
C---- General error reporting subroutine, for the MTZ routines, etc
C
C---- Arguments:
C
C     ERRFLG  (I)  INTEGER         =1 output meesage as warning
C                                  =2 output message as fatal
C
C     IFAIL   (I)  INTEGER         =0 return after fatal error
C                                  =-1 STOP after reporting fatal error
C
C     ERRMSG  (I)  CHARACTER*(*)   character string containing error
C                                  message to output
C
C_END_LERROR
C
C     .. Scalar Arguments ..
      INTEGER ERRFLG,IFAIL
      CHARACTER ERRMSG* (*)
      STOP 'LERROR not implemented'
CC     ..
CC     ..
CC     .. External Subroutines ..
C      EXTERNAL BLANK,PUTLIN
CC     ..
CC
C      IF (ERRFLG.EQ.1) THEN
CC
CC---- Output a warning message and return
CC
C        CALL BLANK('ERRWIN',1)
C        CALL PUTLIN('***  Warning','ERRWIN')
C        CALL PUTLIN(ERRMSG,'ERRWIN')
C        CALL BLANK('ERRWIN',1)
CC
C      ELSE IF (ERRFLG.EQ.2) THEN
CC
CC---- Output a fatal message, and quit or return depending on IFAIL
CC
C        CALL BLANK('ERRWIN',1)
C        CALL PUTLIN('***  Error','ERRWIN')
C        CALL PUTLIN(ERRMSG,'ERRWIN')
C        IF (IFAIL.LT.0) THEN
C          call ccperr(1,'*** Program Terminated ')
C        ELSE
C          CALL BLANK('ERRWIN',1)
C        END IF
C        RETURN
C      ELSE
CC
CC---- Bad errflg, output message and continue
CC
C        CALL BLANK('ERRWIN',1)
C        CALL PUTLIN('*** Unrecognised  error','ERRWIN')
C        CALL PUTLIN(ERRMSG,'ERRWIN')
C        CALL PUTLIN('Program continuing ...','ERRWIN')
C        CALL BLANK('ERRWIN',1)
CC
C      END IF
      END

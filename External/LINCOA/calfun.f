C **********************************************************************
C This file was sent to me upon my request to Professor M. Powell in
C 9/2014. This file did not come with any copyright
C message, but we want to acknowledge that this is a product of
C Professor Powell. We made changes to this file to make it 
C compatible with our software interface.
C Reference:
C M.J.D. Powell (2014), "On fast trust region methods for quadratic
C     models with linear constraints", Report No. DAMTP 2014/NA02,
C     University of Cambridge.
C **********************************************************************
C Changes made:
C   - change name from CALFUN to LCALFUN
C   - change body
C **********************************************************************
      SUBROUTINE LCALFUN (N,X,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON FMODE, OFMODE, KFMODE, NFMODE, LFMODE
C     COMMON FMAX
      IF (LFMODE .EQ. 9999) THEN
         CALL LincoaEvalFunc(N,X,F)
      ELSE IF (LFMODE .EQ. 8888) THEN
         CALL OUUEvalFunc(N, X, F);
      ENDIF
C     ZERO=0.0D0
C     F=FMAX
C     V12=X(1)*X(5)-X(4)*X(2)
C     V13=X(1)*X(8)-X(7)*X(2)
C     V14=X(1)*X(11)-X(10)*X(2)
C     V23=X(4)*X(8)-X(7)*X(5)
C     V24=X(4)*X(11)-X(10)*X(5)
C     V34=X(7)*X(11)-X(10)*X(8)
C     DEL1=V23*X(12)-V24*X(9)+V34*X(6)
C     IF (DEL1 .LE. ZERO) GOTO 10
C     DEL2=-V34*X(3)-V13*X(12)+V14*X(9)
C     IF (DEL2 .LE. ZERO) GOTO 10
C     DEL3=-V14*X(6)+V24*X(3)+V12*X(12)
C     IF (DEL3 .LE. ZERO) GOTO 10
C     DEL4=-V12*X(9)+V13*X(6)-V23*X(3)
C     IF (DEL4 .LE. ZERO) GOTO 10
C     TEMP=(DEL1+DEL2+DEL3+DEL4)**3/(DEL1*DEL2*DEL3*DEL4)
C     F=DMIN1(TEMP/6.0D0,FMAX)
C  10 CONTINUE
      RETURN
      END


C-----------------------------------------------------------------------
C Plot the samples input pair 0, 1 during run
C number represents which processor is running the sample
C circle means the sample has been completed
C-----------------------------------------------------------------------
      SUBROUTINE PLOTSAMPLES2D(N,X2,Y2,NSTAT,DXMIN,DXMAX,DYMIN,DYMAX,
     +                         FLAG) 
      INTEGER I, J, IER, N, NMAX, NSTAT(*)
      PARAMETER (NMAX=20000)
      REAL*8  X2(*),Y2(*), DXMIN, DXMAX, DYMIN, DYMAX
      REAL*4  X(NMAX), Y(NMAX), XMIN, XMAX, YMIN, YMAX
      INTEGER PGBEG, IUNIT, FLAG
C
      XMIN = DXMIN
      XMAX = DXMAX
      YMIN = DYMIN
      YMAX = DYMAX
      IF (N .GT. NMAX) THEN
	PRINT *, 'PGPLOT2D : dimension too large.'
	RETURN
      END IF
C
C Read input file
C
C     IUNIT = 10
C     OPEN(IUNIT,FILE="pgplot.dat")
C     READ(IUNIT,21) FLAG
C21    FORMAT(I1)
C
C Call PGBEG to initiate PGPLOT and open the output device; PGBEG
C will prompt the user to supply the device name and type. Always
C check the return code from PGBEG.
C
      IF (FLAG .EQ. 0) THEN
         IF (PGBEG(0,'/XWINDOW',1,1) .NE. 1) STOP
      END IF
C     READ(IUNIT,22) XMIN,XMAX
C     READ(IUNIT,22) YMIN,YMAX
C     CLOSE(IUNIT)
C22    FORMAT(1X,E20.10,1X,E20.10)
C
C Print information about device.
C
C     CALL PGEX0
C
C Override default colors.
C
      CALL PGSCRN(0, 'DarkSlateGray', IER)
      CALL PGSCRN(1, 'White', IER)
      CALL PGSCRN(2, 'Yellow', IER)
      CALL PGSCRN(3, 'Cyan', IER)
      CALL PGSCRN(4, 'SlateGray', IER)
C
C "Erase" the screen to fill with background color.
C
C     CALL PGERAS
C
C Set up window and viewport.
C
      CALL PGSCH(1.5)
      CALL PGVSTD
C
C Convert data from double precision to single precision
C
      DO I = 1, N
	X(I) = X2(I)
	Y(I) = Y2(I)
      END DO
      CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
C
C Fill the viewport in a different color.
C
      CALL PGSCI(4)
      CALL PGRECT(XMIN, XMAX, YMIN, YMAX)
C
C Annotation.
C
C     CALL PGLAB('(x)', '(i)', 'PGPLOT Data')
      CALL PGSCI(0)
      CALL PGBOX('G', 0.0, 0, 'G', 0.0, 0)
      CALL PGSCI(1)
      CALL PGSLW(3)
      CALL PGSCF(2)
      CALL PGBOX('BCNST', 0.0, 0, 'BCNSTV', 0.0, 0)
C     CALL PGLAB('x', ' ', 'PGPLOT Graph')
C
C Mark N points (coordinates in arrays X and Y), using symbol number 18.
C
C     CALL PGSCI(2)
C     CALL PGSLW(4)
C     CALL PGLINE(N,X,Y)
      CALL PGSCI(3)
      CALL PGSCH(3.0)
      DO I = 1, N
        J = NSTAT(I)
        IF (J .EQ. 0) THEN
          CALL PGPT(1,X(I),Y(I),17)
        ELSE IF (J .LT. 0) THEN
          J = J - J / 10 * 10
          CALL PGPT(1,X(I),Y(I),48-J)
        ELSE
          J = J - J / 10 * 10
          CALL PGPT(1,X(I),Y(I),48 + J)
          CALL PGPT(1,X(I),Y(I),25)
        END IF
      END DO
C     IF (FLAG .NE. 0) THEN
C       CALL PGLINE(N,X,YS)
C       CALL PGPT(N,X,YS,9)
C     END IF
C
C Finally write data back to file
C
C     IUNIT = 10
      FLAG = 1
C     OPEN(IUNIT,FILE="pgplot.dat")
C     WRITE(IUNIT,21) FLAG
C     WRITE(IUNIT,22) XMIN,XMAX
C     WRITE(IUNIT,22) YMIN,YMAX
C     DO I = 1, N
C        WRITE(IUNIT,23) X(I), Y(I) 
C     END DO
C     CLOSE(IUNIT)
      END
C-----------------------------------------------------------------------
C Plot 2D scatterplot
C-----------------------------------------------------------------------
      SUBROUTINE PLOTSCATTER2D(N,X2,Y2,DXMIN,DXMAX,DYMIN,DYMAX,FLAG) 
      INTEGER I, J, IER, N, NMAX, N2, NBINS
      PARAMETER (NMAX=20000)
      REAL*8  X2(*),Y2(*), DXMIN, DXMAX, DYMIN, DYMAX
      REAL*4  X(NMAX), Y(NMAX), XMIN, XMAX, YMIN, YMAX, XMEANS(NMAX)
      REAL*4  YMEANS(NMAX),  YMAXS(NMAX), YMINS(NMAX)
      INTEGER PGBEG, IUNIT, FLAG
C
      XMIN = DXMIN
      XMAX = DXMAX
      YMIN = DYMIN
      YMAX = DYMAX
      IF (N .GT. NMAX) THEN
	PRINT *, 'PGPLOT2D : dimension too large.'
	RETURN
      END IF
C
C Call PGBEG to initiate PGPLOT and open the output device; PGBEG
C will prompt the user to supply the device name and type. Always
C check the return code from PGBEG.
C
      IF (FLAG .EQ. 0) THEN
         IF (PGBEG(0,'/XWINDOW',1,1) .NE. 1) STOP
      END IF
C
C Print information about device.
C
C     CALL PGEX0
C
C Override default colors.
C
C     CALL PGSCRN(0, 'DarkSlateGray', IER)
      CALL PGSCRN(0, 'Black', IER)
      CALL PGSCRN(1, 'White', IER)
      CALL PGSCRN(2, 'Yellow', IER)
      CALL PGSCRN(3, 'Cyan', IER)
      CALL PGSCRN(4, 'SlateGray', IER)
C
C "Erase" the screen to fill with background color.
C
C     CALL PGERAS
C
C Set up window and viewport.
C
      CALL PGSCH(1.5)
      CALL PGSCI(4)
      CALL PGVSTD
C
C Convert data from double precision to single precision
C
      N2 = 1
      NBINS = 1
      XMEANS(1) = X2(1)
      YMEANS(1) = Y2(1)
      YMAXS(1)  = Y2(1)
      YMINS(1)  = Y2(1)
      DO I = 1, N
	 X(I) = X2(I)
	 Y(I) = Y2(I)
         IF (I .GT. 1) THEN
            IF (ABS(X(I) - X(I-1)) .LT. 1.0D-10) THEN
               N2 = N2 + 1
               YMEANS(NBINS) = YMEANS(NBINS) + Y(I)
               IF (Y(I) .GT. YMAXS(NBINS)) YMAXS(NBINS) = Y(I)
               IF (Y(I) .LT. YMINS(NBINS)) YMINS(NBINS) = Y(I)
            ELSE
               YMEANS(NBINS) = YMEANS(NBINS) / N2
               N2 = 1
               NBINS = NBINS + 1
               XMEANS(NBINS) = X(I)
               YMEANS(NBINS) = Y(I)
               YMAXS(NBINS)  = Y(I)
               YMINS(NBINS)  = Y(I)
            END IF
         END IF 
      END DO
      YMEANS(NBINS) = YMEANS(NBINS) / N2
      CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
C
C Fill the viewport in a different color.
C
      CALL PGSCI(1)
      CALL PGRECT(XMIN, XMAX, YMIN, YMAX)
C
C Annotation.
C
C grid line color (4 - black)
      CALL PGSCI(4)
      CALL PGBOX('G', 0.0, 0, 'G', 0.0, 0)
C tick mark color
      CALL PGSCI(4)
C character line thickness
      CALL PGSLW(5)
C font (2 - roman)
      CALL PGSCF(2)
      CALL PGSCI(4)
      CALL PGBOX('BCNST', 0.0, 0, 'BCNSTV', 0.0, 0)
      CALL PGSCI(2)
      CALL PGLAB('x', 'y', 'PGPLOT Scatter Plot')
C
C Mark N points (coordinates in arrays X and Y), using symbol number 18.
C
C point marker color (14 - black, 120 - x, size = 1.3)
      CALL PGSCI(14)
      CALL PGSCH(1.3)
      DO I = 1, N
         CALL PGPT(1,X(I),Y(I),120)
      END DO
C point marker color (12 - purple, 42 - asterisk)
      CALL PGSCI(11)
      CALL PGSCH(3.0)
      DO I = 1, NBINS
         CALL PGPT(1,XMEANS(I),YMEANS(I),42)
      END DO
C line (11 - blue)
      CALL PGSLW(6)
      CALL PGSCI(11)
      CALL PGLINE(NBINS,XMEANS,YMEANS)
C line (8 - red)
      CALL PGSCI(8)
      CALL PGLINE(NBINS,XMEANS,YMAXS)
      CALL PGLINE(NBINS,XMEANS,YMINS)
C     IF (FLAG .NE. 0) THEN
C       CALL PGLINE(N,X,YS)
C       CALL PGPT(N,X,YS,9)
C     END IF
C
C Finally write data back to file
C
C     IUNIT = 10
C     FLAG = 1
C     OPEN(IUNIT,FILE="pgplot.dat")
C     WRITE(IUNIT,21) FLAG
C     WRITE(IUNIT,22) XMIN,XMAX
C     WRITE(IUNIT,22) YMIN,YMAX
C     DO I = 1, N
C        WRITE(IUNIT,23) X(I), Y(I) 
C     END DO
C     CLOSE(IUNIT)
      END
C-----------------------------------------------------------------------
C Plot 2D multiple scatterplot
C X2 - the first variable that varies the slowest
C X3 - the second variable that varies faster 
C-----------------------------------------------------------------------
      SUBROUTINE PLOTSCATTERM2D(N,X2,X3,Y2,DXMIN,DXMAX,DYMIN,DYMAX,FLAG) 
      INTEGER I, J, IER, N, NMAX, NCOLOR, NCNT, NCUR
      PARAMETER (NMAX=20000)
      REAL*8  X2(*), X3(*), Y2(*), DXMIN, DXMAX, DYMIN, DYMAX
      REAL*4  X(NMAX), Y(NMAX), XC(NMAX), YMIN, YMAX
      REAL*4  XMIN, XMAX, YC(NMAX)
      INTEGER PGBEG, IUNIT, FLAG
C
C convert min and max to single precision 
C
      XMIN = DXMIN
      XMAX = DXMAX
      YMIN = DYMIN
      YMAX = DYMAX
      IF (N .GT. NMAX) THEN
	PRINT *, 'PGPLOTM2D : dimension too large.'
	RETURN
      END IF
C
C Call PGBEG to initiate PGPLOT and open the output device; PGBEG
C will prompt the user to supply the device name and type. Always
C check the return code from PGBEG.
C
      IF (FLAG .EQ. 0) THEN
         IF (PGBEG(0,'/XWINDOW',1,1) .NE. 1) STOP
      END IF
C
C Print information about device.
C
C     CALL PGEX0
C
C Override default colors.
C
C     CALL PGSCRN(0, 'DarkSlateGray', IER)
      CALL PGSCRN(0, 'Black', IER)
      CALL PGSCRN(1, 'White', IER)
      CALL PGSCRN(2, 'Yellow', IER)
      CALL PGSCRN(3, 'Cyan', IER)
      CALL PGSCRN(4, 'SlateGray', IER)
C
C "Erase" the screen to fill with background color.
C
C     CALL PGERAS
C
C Set up window and viewport.
C
      CALL PGSCH(1.5)
      CALL PGSCI(4)
      CALL PGVSTD
C
C Fill the viewport in a different color.
C
      CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
      CALL PGSCI(1)
      CALL PGRECT(XMIN, XMAX, YMIN, YMAX)
C
C Annotation.
C
C grid line color (4 - black)
      CALL PGSCI(4)
      CALL PGBOX('G', 0.0, 0, 'G', 0.0, 0)
C tick mark color
      CALL PGSCI(4)
C character line thickness
      CALL PGSLW(5)
C font (2 - roman)
      CALL PGSCF(2)
      CALL PGSCI(4)
      CALL PGBOX('BCNST', 0.0, 0, 'BCNSTV', 0.0, 0)
      CALL PGSCI(2)
      CALL PGLAB('x', 'y', 'PGPLOT Scatter2D Plot')
C
C point marker color (14 - black, 120 - x, size = 1.3)
      NCOLOR = 4
      CALL PGSCI(NCOLOR)
      CALL PGSCH(1.3)
C
C Convert data from double precision to single precision
C X has the faster varying variable
C
      DO I = 1, N
	 X(I) = X3(I)
	 Y(I) = Y2(I)
      END DO
C
C Convert data from double precision to single precision
C
      NCNT = 1
      YC(NCNT) = Y(1)
      XC(NCNT) = X(1)
      NCUR = 1
      DO I = 1, N
         IF (I .GT. 1) THEN
            IF (ABS(X2(I) - X2(I-1)) .LT. 1.0D-10) THEN
               IF (ABS(X(I) - XC(NCNT)) .LT. 1.0D-10) THEN
                  YC(NCNT) = YC(NCNT) + Y(I) 
                  NCUR = NCUR + 1
               ELSE
                  YC(NCNT) = YC(NCNT) / NCUR
                  CALL PGPT(1,XC(NCNT),YC(NCNT),120)
                  NCNT = NCNT + 1
                  XC(NCNT) = X(I)
                  YC(NCNT) = Y(I)
                  NCUR = 1
               END IF
            ELSE
               YC(NCNT) = YC(NCNT) / NCUR
               CALL PGPT(1,XC(NCNT),YC(NCNT),120)
               CALL PGLINE(NCNT,XC,YC)
               NCNT = 1
               YC(NCNT) = Y(I)
               XC(NCNT) = X(I)
               NCUR = 1
               NCOLOR = NCOLOR + 1
               IF (NCOLOR .GT. 14) THEN
                  NCOLOR = 4
               END IF
               CALL PGSCI(NCOLOR)
            END IF
         END IF 
      END DO 
               
C Finally write data back to file
C
C     IUNIT = 10
C     FLAG = 1
C     OPEN(IUNIT,FILE="pgplot.dat")
C     WRITE(IUNIT,21) FLAG
C     WRITE(IUNIT,22) XMIN,XMAX
C     WRITE(IUNIT,22) YMIN,YMAX
C     DO I = 1, N
C        WRITE(IUNIT,23) X(I), Y(I) 
C     END DO
C     CLOSE(IUNIT)
      END
C-----------------------------------------------------------------------
C call PGEND to terminate things properly.
C-----------------------------------------------------------------------
      SUBROUTINE PLOT2DIRIX(N,X2,Y2,NSTAT,DXMIN,DXMAX,DYMIN,DYMAX,FLAG) 
      INTEGER I, J, IER, N, NMAX, NSTAT(*)
      PARAMETER (NMAX=20000)
      REAL*8  X2(*),Y2(*), DXMIN, DXMAX, DYMIN, DYMAX
      REAL*4  X(NMAX), Y(NMAX), XMIN, XMAX, YMIN, YMAX
      INTEGER PGBEG, IUNIT, FLAG
C
      XMIN = DXMIN
      XMAX = DXMAX
      YMIN = DYMIN
      YMAX = DYMAX
      IF (N .GT. NMAX) THEN
	PRINT *, 'PGPLOT2D : dimension too large.'
	RETURN
      END IF
C
C Read input file
C
C     IUNIT = 10
C     OPEN(IUNIT,FILE="pgplot.dat")
C     READ(IUNIT,21) FLAG
C21    FORMAT(I1)
C
C Call PGBEG to initiate PGPLOT and open the output device; PGBEG
C will prompt the user to supply the device name and type. Always
C check the return code from PGBEG.
C
      IF (FLAG .EQ. 0) THEN
         IF (PGBEG(0,'/XWINDOW',1,1) .NE. 1) STOP
      END IF
C     READ(IUNIT,22) XMIN,XMAX
C     READ(IUNIT,22) YMIN,YMAX
C     CLOSE(IUNIT)
C22    FORMAT(1X,E20.10,1X,E20.10)
C
C Print information about device.
C
C     CALL PGEX0
C
C Override default colors.
C
      CALL PGSCRN(0, 'DarkSlateGray', IER)
      CALL PGSCRN(1, 'White', IER)
      CALL PGSCRN(2, 'Yellow', IER)
      CALL PGSCRN(3, 'Cyan', IER)
      CALL PGSCRN(4, 'SlateGray', IER)
C
C "Erase" the screen to fill with background color.
C
C     CALL PGERAS
C
C Set up window and viewport.
C
      CALL PGSCH(1.5)
      CALL PGVSTD
C
C Convert data from double precision to single precision
C
      DO I = 1, N
	X(I) = X2(I)
	Y(I) = Y2(I)
      END DO
      CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
C
C Fill the viewport in a different color.
C
      CALL PGSCI(4)
      CALL PGRECT(XMIN, XMAX, YMIN, YMAX)
C
C Annotation.
C
C     CALL PGLAB('(x)', '(i)', 'PGPLOT Data')
      CALL PGSCI(0)
      CALL PGBOX('G', 0.0, 0, 'G', 0.0, 0)
      CALL PGSCI(1)
      CALL PGSLW(3)
      CALL PGSCF(2)
      CALL PGBOX('BCNST', 0.0, 0, 'BCNSTV', 0.0, 0)
C     CALL PGLAB('x', ' ', 'PGPLOT Graph')
C
C Mark N points (coordinates in arrays X and Y), using symbol number 18.
C
C     CALL PGSCI(2)
C     CALL PGSLW(4)
C     CALL PGLINE(N,X,Y)
      CALL PGSCI(3)
      CALL PGSCH(3.0)
      DO I = 1, N
        J = NSTAT(I)
        CALL PGERR1(5,X(I),Y(I),0.0,1.0)
        CALL PGERR1(6,X(I),Y(I),0.0,1.0)
      END DO
C     IF (FLAG .NE. 0) THEN
C       CALL PGLINE(N,X,YS)
C       CALL PGPT(N,X,YS,9)
C     END IF
C
C Finally write data back to file
C
C     IUNIT = 10
      FLAG = 1
C     OPEN(IUNIT,FILE="pgplot.dat")
C     WRITE(IUNIT,21) FLAG
C     WRITE(IUNIT,22) XMIN,XMAX
C     WRITE(IUNIT,22) YMIN,YMAX
C     DO I = 1, N
C        WRITE(IUNIT,23) X(I), Y(I) 
C     END DO
C     CLOSE(IUNIT)
      END
C-----------------------------------------------------------------------
C call PGEND to terminate things properly.
C-----------------------------------------------------------------------
      SUBROUTINE PLOTEND
      PRINT *, 'type enter to continue'
      CALL PGEND
      END
C-----------------------------------------------------------------------

      SUBROUTINE PGEX0
C-----------------------------------------------------------------------
C This subroutine tests PGQINF and displays the information returned on
C the standard output.
C-----------------------------------------------------------------------
      CHARACTER*64 VALUE
      INTEGER LENGTH
      REAL X, Y, X1, X2, Y1, Y2
C
C Information available from PGQINF:
C
      CALL PGQINF('version',  VALUE, LENGTH)
      WRITE (*,*) 'version=', VALUE(:LENGTH)
      CALL PGQINF('state',    VALUE, LENGTH)
      WRITE (*,*) 'state=',   VALUE(:LENGTH)
      CALL PGQINF('user',     VALUE, LENGTH)
      WRITE (*,*) 'user=',    VALUE(:LENGTH)
      CALL PGQINF('now',      VALUE, LENGTH)
      WRITE (*,*) 'now=',     VALUE(:LENGTH)
      CALL PGQINF('device',   VALUE, LENGTH)
      WRITE (*,*) 'device=',  VALUE(:LENGTH)
      CALL PGQINF('file',     VALUE, LENGTH)
      WRITE (*,*) 'file=',    VALUE(:LENGTH)
      CALL PGQINF('type',     VALUE, LENGTH)
      WRITE (*,*) 'type=',    VALUE(:LENGTH)
      CALL PGQINF('dev/type', VALUE, LENGTH)
      WRITE (*,*) 'dev/type=',VALUE(:LENGTH)
      CALL PGQINF('hardcopy', VALUE, LENGTH)
      WRITE (*,*) 'hardcopy=',VALUE(:LENGTH)
      CALL PGQINF('terminal', VALUE, LENGTH)
      WRITE (*,*) 'terminal=',VALUE(:LENGTH)
      CALL PGQINF('cursor',   VALUE, LENGTH)
      WRITE (*,*) 'cursor=',  VALUE(:LENGTH)
C
C Get view surface dimensions:
C
      CALL PGQVSZ(1, X1, X2, Y1, Y2)
      X = X2-X1
      Y = Y2-Y1
      WRITE (*,100) X, Y, X*25.4, Y*25.4
  100 FORMAT (' Plot dimensions (x,y; inches): ',F9.2,', ',F9.2/
     1        '                          (mm): ',F9.2,', ',F9.2)
C-----------------------------------------------------------------------
      END

C-----------------------------------------------------------------------
      SUBROUTINE PLOT3D(N,Z2,FLAG) 
      INTEGER N, FLAG, NMAX
      PARAMETER (NMAX=20000)
      REAL*8  Z2(*), DSQRT, DTMP
      REAL*4  ZS(NMAX)
      REAL*4  ZMIN, ZMAX, SIZE
      INTEGER I, IER, PGBEG, IX, NUMCOL
C
C Check input size 
C
      IF (N .GT. NMAX) THEN
	 PRINT *, 'PGPLOT : buffer size not large enough, quit.'
         RETURN
      END IF
C
C check input
C
      DTMP = 1.0 * N
      DTMP = DSQRT(DTMP)
      IX   = INT(DTMP)
      IF (IX * IX .NE. N) THEN
	 PRINT *, 'PGPLOT : N should be a perfect square.'
         RETURN
      END IF
C
C Call PGBEG to initiate PGPLOT and open the output device; PGBEG
C will prompt the user to supply the device name and type. Always
C check the return code from PGBEG.
C
      IF (FLAG .EQ. 0) THEN
         IF (PGBEG(0,'/XWINDOW',1,1) .NE. 1) STOP
	 SIZE = 1.0
         CALL PGENV(0.,SIZE,0.,SIZE,1, -2)
      END IF
C
C set color
      NUMCOL = 20
      CALL PGSCR(2, 0.0, 0.0, 1.0)
      CALL PGSCR(3, 0.0, 0.1, 1.0)
      CALL PGSCR(4, 0.0, 0.2, 1.0)
      CALL PGSCR(5, 0.0, 0.3, 1.0)
      CALL PGSCR(6, 0.0, 0.4, 1.0)
      CALL PGSCR(7, 0.0, 0.5,  1.0)
      CALL PGSCR(8, 0.0, 0.6,  1.0)
      CALL PGSCR(9, 0.0, 0.7, 1.0)
      CALL PGSCR(10, 0.0, 0.8, 1.0)
      CALL PGSCR(11, 0.0, 0.9, 1.0)
      CALL PGSCR(12, 0.0, 1.0, 1.0)
      CALL PGSCR(13, 0.0, 1.0, 0.9)
      CALL PGSCR(14, 0.0, 1.0, 0.8)
      CALL PGSCR(15, 0.0, 1.0, 0.7)
      CALL PGSCR(16, 0.0, 1.0, 0.6)
      CALL PGSCR(17, 0.0, 1.0, 0.5)
      CALL PGSCR(18, 0.0, 1.0, 0.4)
      CALL PGSCR(19, 0.0, 1.0, 0.3)
      CALL PGSCR(20, 0.0, 1.0, 0.2)
      CALL PGSCR(21, 0.0, 1.0, 0.1)
      CALL PGSCR(22, 0.0, 1.0, 0.0)
      CALL PGSCR(23, 0.2, 1.0, 0.0)
      CALL PGSCR(24, 0.3, 1.0, 0.0)
      CALL PGSCR(25, 0.4, 1.0, 0.0)
      CALL PGSCR(26, 0.5, 1.0, 0.0)
      CALL PGSCR(27, 0.6, 1.0, 0.0)
      CALL PGSCR(28, 0.7, 1.0, 0.0)
      CALL PGSCR(29, 0.8, 1.0, 0.0)
      CALL PGSCR(30, 0.9, 1.0, 0.0)
      CALL PGSCR(31, 1.0, 1.0, 0.0)
      CALL PGSCR(32, 1.0, 0.9, 0.0)
      CALL PGSCR(33, 1.0, 0.8, 0.0)
      CALL PGSCR(34, 1.0, 0.7, 0.0)
      CALL PGSCR(35, 1.0, 0.6, 0.0)
      CALL PGSCR(36, 1.0, 0.5, 0.0)
      CALL PGSCR(37, 1.0, 0.4, 0.0)
      CALL PGSCR(38, 1.0, 0.3, 0.0)
      CALL PGSCR(39, 1.0, 0.2, 0.0)
      CALL PGSCR(40, 1.0, 0.1, 0.0)
      CALL PGSCR(41, 1.0, 0.0, 0.0)
      CALL PGSCI(3)
C     CALL PGLAB('Test X', 'Test Y', 'Test Title')
C
C Convert data from double precision to single precision
C
      ZMIN =  1.0E14
      ZMAX = -1.0E14
      DO I = 1, N
        IF (Z2(I) .GT. ZMAX) ZMAX = Z2(I)
        IF (Z2(I) .LT. ZMIN) ZMIN = Z2(I)
      END DO
      IF ((ZMAX - ZMIN) .LT. 1.0E-10) THEN
         ZMIN = ZMIN - 1.0E-10
         ZMAX = ZMAX + 1.0E-10
      END IF
      DTMP = ZMAX - ZMIN
      DTMP = 1.0 / DTMP
      DO I = 1, N
         ZS(I) = Z2(I) * DTMP
      END DO
C
C For the mesh option
C
C     CALL PGERAS
      SIZE = 1.0
C     CALL PGENV(0.,SIZE,0.,SIZE,1, -2)
      CALL FREDDY(ZS,IX,IX,SIZE,40.0)
C
C update PGPLOT variable write data back to file
C
      FLAG = 1
      END

C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE FREDDY(ARRAY,KX,NY,SIZE,ANGLE)
      INTEGER KX, NY
      REAL    ARRAY(KX,NY), SIZE, ANGLE
      INTEGER IFLAG(200, 200), FLAG
      REAL    HWIDTH
C
C Draws isometric plot of array
C
      REAL FMAX,FMIN,DELTAX,DELTAY,DELTAV,SINE,PEAK,X,DX,HEIGHT
      INTEGER I,J,KI,KJ,NX,MX,MY,STEP,LEFT,RIGHT,IT,MN,INCX
      LOGICAL VISBLE
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
      IF (KX .GT. 200 .or. NY .GT. 200) THEN
         FLAG = 0
      ELSE
         FLAG = 1
      END IF
      MN = KX*NY
      NX = KX
C     Check array size:
      IF(NX.LT.2 .OR. NY.LT.2) RETURN
      FMAX = ARRAY(1,1)
      FMIN = FMAX
      DO 20 J=1,NY
          DO 10 I=1,NX
              FMIN = AMIN1(ARRAY(I,J),FMIN)
              FMAX = AMAX1(ARRAY(I,J),FMAX)
   10     CONTINUE
   20 CONTINUE
      DELTAX = SIZE/(NX+NY)
      SINE = SIN(ANGLE/58.)
      DELTAY = DELTAX*SINE
      HEIGHT = SIZE*(1.-ABS(SINE))
      DELTAV = HEIGHT
      FMAX = FMAX-FMIN
      IF(FMAX.LT.0.0001) FMAX = DELTAV
      HWIDTH = FMAX / 39
      IF (FLAG .EQ. 1) THEN
         DO I = 1, KX
            DO J = 1, NY
               IFLAG(I,J) = (ARRAY(I,J) - FMIN) / HWIDTH
               IFLAG(I,J) = IFLAG(I,J) + 2
               IF (IFLAG(I,J) .LT. 2)  IFLAG(I,J) = 2
               IF (IFLAG(I,J) .GT. 41) IFLAG(I,J) = 41
            END DO
         END DO
      ELSE
         DO I = 1, KX
            DO J = 1, NY
               IFLAG(I,J) = 1
            END DO
         END DO
      END IF
      DELTAV = DELTAV/FMAX
      MX = NX+1
      MY = NY+1
      STEP = MX
C
C Start PGPLOT buffering.
C
      CALL PGBBUF
      CALL PGERAS
      CALL PGSCIR(1,32)
      CALL PGSCI(1)
C
C Work our way down the Y axis, then up the X axis,
C calculating the Y plotter coordinates for each
C column of the plot, doing the hidden-line suppression
C at the same time.
C
      DO 50 J=1,NY
          KJ = MY-J
          KI = 1
C               ( KI,KJ are coordinates of bottom of column)
          ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
   30     PEAK = ARRAY(KI,KJ)
   40     KI = KI+1
          KJ = KJ+1
          IF(KI.GT.NX .OR. KJ.GT.NY) GOTO 50
          ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(ARRAY(KI,KJ).GT.PEAK) GOTO 30
          IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          GOTO 40
   50 CONTINUE
C
C Now to work our way up the X axis
C
      DO 80 I=2,NX
          KI = I
          KJ = 1
          ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
   60     PEAK = ARRAY(KI,KJ)
   70     KI = KI+1
          KJ = KJ+1
          IF(KI.GT.NX .OR. KJ.GT.NY) GOTO 80
          ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(ARRAY(KI,KJ).GT.PEAK) GOTO 60
          IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          GOTO 70
   80 CONTINUE
C
C Draw a line along the bottom of the vertical faces
C
      CALL PGMOVE(DELTAX*(NX+NY-2), DELTAY*(MX))
      CALL PGDRAW(DELTAX*(NY-1),    DELTAY*2)
      CALL PGDRAW(0.0,              DELTAY*MY)
C
C Array is now ready for plotting.  If a point is
C positive, then it is to be plotted at that Y
C coordinate; if it is negative, then it is
C invisible, but at minus that Y coordinate (the point
C where the line heading towards it disappears has to
C be determined by finding the intersection of it and
C the cresting line).
C
C Plot rows:
C
      CALL PGSCI(1)
      DO 110 J=1,NY,2
          KJ = MY-J
          DX = DELTAX*(J-2)
          X = DX+DELTAX
          CALL PGMOVE(X,DELTAY*(KJ+1))
          CALL PGDRAW(X,ARRAY(1,KJ))
          VISBLE = .TRUE.
          DO 90 I=2,NX
              RIGHT = I+NX*(KJ-1)
              LEFT = RIGHT-1
              IT = RIGHT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN,IFLAG(I,KJ))
              CALL PGSCI(1)
   90     CONTINUE
C
C Now at far end of row so come back
C
          KJ = KJ-1
          IF(KJ.LE.0) GOTO 170
          VISBLE = ARRAY(NX,KJ).GE.0.0
          DX = DELTAX*(NX+J)
          IF(VISBLE) CALL PGMOVE(DX-DELTAX,ARRAY(NX,KJ))
          DELTAX = -DELTAX
          DO 100 I=2,NX
              KI = MX-I
              LEFT = KI+NX*(KJ-1)
              RIGHT = LEFT+1
              IT = LEFT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN,IFLAG(KI,KJ))
              CALL PGSCI(1)
  100     CONTINUE
C
          X = DX+DELTAX*NX
          IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(1,KJ))
          CALL PGDRAW(X,DELTAY*(KJ+1))
C               (set DELTAX positive for return trip)
          DELTAX = -DELTAX
  110 CONTINUE
C
C Now do the columns:
C as we fell out of the last DO-loop we do the
C columns in ascending-X order
C
      INCX = 1
      KI = 1
C               (set DELTAX -ve since scanning R to L)
      CALL PGSCI(2)
  120 DX = DELTAX*(KI+NY-1)
      DELTAX = -DELTAX
      X = DX+DELTAX
      CALL PGMOVE(X,ARRAY(1,1))
  130 VISBLE = .TRUE.
      DO 140 J=2,NY
          LEFT = KI+NX*(J-1)
          RIGHT = LEFT-NX
          IT = LEFT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN,IFLAG(KI,J))
          CALL PGSCI(1)
  140 CONTINUE
C
C At far end, increment X and check still inside array
C
      CALL PGSCI(6)
      KI = KI+INCX
      IF(KI.LE.0 .OR. KI.GT.NX) GOTO 180
      VISBLE = ARRAY(KI,NY).GE.0.0
      DELTAX = -DELTAX
      DX = DELTAX*(KI-2)
      X = DX+DELTAX
      IF(VISBLE) CALL PGMOVE(X,ARRAY(KI,NY))
      DO 150 J=2,NY
          KJ = MY-J
          RIGHT = KI+NX*(KJ-1)
          LEFT = RIGHT+NX
          IT = RIGHT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN,IFLAG(KI,KJ))
          CALL PGSCI(1)
  150 CONTINUE
C
      X = DX+DELTAX*NY
      IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(KI,1))
      IF(KI.EQ.1) GOTO 180
      CALL PGDRAW(X,DELTAY*(KI+1))
      KI = KI+INCX
      IF(KI.GT.NX) GOTO 180
      IF(KI.EQ.1) GOTO 120
  160 DELTAX = -DELTAX
      DX = DELTAX*(1-KI-NY)
      X = DX+DELTAX
      CALL PGMOVE(X,DELTAY*(KI+1))
      CALL PGDRAW(X,ARRAY(KI,1))
      GOTO 130
C
C Do columns backwards because ended rows at far end of X
C
  170 KI = NX
      INCX = -1
      DX = DELTAX*(KI+NY)
      GOTO 160
C
C
  180 CALL PGEBUF
      END
C-----------------------------------------------------------------------
      SUBROUTINE FREDGO(ARRAY,MN,COLOR)
      INTEGER MN, COLOR
      REAL ARRAY(MN)
C
      INTEGER STEP,LEFT,RIGHT,IT,NX
      LOGICAL VISBLE
      REAL AL,AR,BL,EM,XX,X,Y,DELTAX
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
C Test visibility
C
      CALL PGSCI(COLOR)
      IF(ARRAY(IT).LT.0.0) GOTO 80
C
C This point is visible - was last?
C
      IF(VISBLE) GOTO 50
C
C No: calculate point where this line vanishes
C
   10 IF(LEFT.LE.NX .OR. MOD(LEFT-1,NX).EQ.0 .OR.
     1     RIGHT.LE.NX .OR. MOD(RIGHT-1,NX).EQ.0) GOTO 100
      AL = ABS(ARRAY(LEFT))
      AR = ABS(ARRAY(RIGHT))
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C               Right-hand point is crested
   20 RIGHT = RIGHT-STEP
      IF(ARRAY(RIGHT).LT.0.0) GOTO 20
C               Left-hand end of cresting line is either
C               RIGHT+NX or RIGHT-1
      LEFT = RIGHT+NX
      IF(ARRAY(LEFT).LT.0.0) LEFT = RIGHT-1
C
C               RIGHT and LEFT index into the endpoints of the
C               cresting line
   30 BL = ABS(ARRAY(LEFT))
      EM = ABS(ARRAY(RIGHT))-BL
      XX = EM-AR+AL
      IF(ABS(XX).LT.0.0001) GOTO 60
      XX = (AL-BL)/XX
   40 Y = EM*XX+BL
      IF(DELTAX.GT.0.0) XX = 1.0-XX
      XX = X-XX*DELTAX
      IF(VISBLE) GOTO 90
C               Drawing a line from an invisible point
C               to a visible one
      CALL PGMOVE(XX,Y)
      VISBLE = .TRUE.
   50 CALL PGDRAW(X,ARRAY(IT))
      RETURN
C
   60 XX = 0.5
      GOTO 40
C
C Left-hand point crested
C
   70 LEFT = LEFT-STEP
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C
C Right-hand end of cresting line is either LEFT+1 or LEFT-NX
C
      RIGHT = LEFT+1
      IF(ARRAY(RIGHT).LT.0.0) RIGHT = LEFT-NX
      GOTO 30
C
C This point is invisible; if last one was too, then forget it;
C else draw a line towards it
C
   80 IF(.NOT.VISBLE) RETURN
      GOTO 10
C
   90 CALL PGDRAW(XX,Y)
  100 VISBLE = .FALSE.
      RETURN
      END

      SUBROUTINE HYBEG
C--BEGIN HYPOINVERSE BY INITIALIZING CERTAIN STRINGS.
C--THIS IS EQUIVALENT TO INITIALIZATION IN A DATA STATEMENT.
      INCLUDE 'common.inc'

C--THINGS FORMERLY SET IN BLOCK DATA
      LCOMP1=.FALSE.
      SUBMOD=.FALSE.
      CDOMAN='NC'
      CPVERS='01'

C--SUN/UNIX
       TERMIN='/dev/tty'
       INFILE(1)='hypinst'
       INFILE(0)='/home/calnet/klein/hypfiles/cal2000.hyp'
       STAFIL='/home/calnet/klein/hypfiles/all2.sta'
C       ATNFIL='/home/calnet/klein/hypfiles/all2.atn'
C       FMCFIL='/home/calnet/klein/hypfiles/all2.fmc'
C       XMCFIL='/home/calnet/klein/hypfiles/all2.xmc'
       BSTAFL='/home/calnet/klein/hypfiles/allsta2.bin'
       BCRUFL='/home/calnet/klein/hypfiles/multmod2.bin'
       ATNFIL=' '
       FMCFIL=' '
       XMCFIL=' '

C--VAX
C      TERMIN='TT:'
C      INFILE(1)='HYPINST.'
C      INFILE(0)='HOME:[KLEIN.HYPFILES]CAL2000.HYP'
CC      INFILE(0)='HYPOINV$MODELS'            !CUSP LOGICAL NAME

C--STATION DATA
      DO J=1,MAXSTA
        JCEXP(J)=0
        JCAL(J)=0
        JLMOD(J)=.FALSE.
      END DO

C--PHASE DATA
      DO K=1,MAXPHS
        KRMK6(K)=' '
        KCAL(K)=0
      END DO

C--NODE DATA
      DO I=1,NODMAX
        MODH(I)=0
      END DO

C--CRUST MODEL DATA
      DO I=1,LH
        POSM(I)=1.75
        POSB(I)=0.
        MODTYP(I)=-1
        MODALT(I)=0
        MODSAL(I)=0
        ELEVMX(I)=0.
        LELEV(I)=.FALSE.
      END DO
      DO I=1,LM
        CRODE(I)='   '
      END DO
      RETURN
      END

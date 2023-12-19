      SUBROUTINE HYDATUM
C--COMPUTES DEPTH DATUM FOR CRT AND CRH MODELS

      INCLUDE 'common.inc'
C--TEMPORARY STATION ARRAYS FOR SORTING & 5 CLOSEST UNIQUE STATIONS
      CHARACTER STASO*9,STASU*9
      DIMENSION STASO(MAXPHS), DISSO(MAXPHS), IELSO(MAXPHS)
      DIMENSION STASU(5), DISSU(5), IELSU(5)

C--DECIDE OF DOMINANT MODEL TYPE (CRT,CRH) NEEDS A DD CALCULATION
      IF (MODTYP(MODS(1)) .GT.1) THEN
C--CRE,CRV,CRL MODELS HAVE A DEPTH DATUM OF 0 (SEA LEVEL OR GEOID)
C--THE DOMINANT MODEL NUMBER IS MODS(1)
        IDEPDAT=0
        JDSTA=0
        IMODG=1
        GOTO 50
      END IF
      IMODG=0

C--GET AVERAGE ELEVATIONS OF THE 5 CLOSEST STATIONS
C--GET STATIONS AND DISTANCES FROM HI ARRAYS. THEY ARE NOT IN DISTANCE ORDER.
C--FILL ARRAYS AND SORT THEM BY DISTANCE, ARRAYS WILL BE REARRANGED IN SORTING
C--USE EACH STATION SNL ONLY ONCE, REGARDLESS OF COMPONENT C
      DO K=1,KSTA
        JINDXX=KINDX(K)
        STASO(K)=(STANAM(JINDXX)//JNET(JINDXX)//JSLOC(JINDXX))
        DISSO(K)=DIS(K)
        IELSO(K)=JELEV(JINDXX)
      END DO
      KKSTA=KSTA
      CALL SORT3 (KKSTA,DISSO,STASO,IELSO)

c--we now have a sorted station list, but there may be duplicate stations
c  or ones with 0 (unknown) elevations
c--put the 5 "unique" stations into the u arrays
c--KS is the index of sorted stations, J is the index of the 5 unique stations
c--every earthquake must have 3 stations
      J=0
      DO KSZ=1,KSTA
c--SKIP STATIONS WITH UNKNOWN ELEVATIONS
        IF (IELSO(KSZ).EQ.0) GOTO 30
        IF (J.EQ.0) THEN
          STASU(1)=STASO(KSZ)
          DISSU(1)=DISSO(KSZ)
          IELSU(1)=IELSO(KSZ)
          J=1
          JDSTA=1
          GOTO 30
        END IF

C--WE HAVE THE CLOSEST STATION AND J.GE.1
C--SEARCH FOR CURRENT STASO STATION IN STASU CLOSE LIST, SKIP IF WE HAVE IT
        DO JJ=1,J
          IF (STASO(KSZ) .EQ. STASU(JJ)) GOTO 30
        END DO

C--WE HAVE A NEW STATION STASO, SO ADD IT TO CLOSEST LIST
        J=J+1
        JDSTA=J
        STASU(J)=STASO(KSZ)
        DISSU(J)=DISSO(KSZ)
        IELSU(J)=IELSO(KSZ)
        IF (J.EQ.5) GOTO 40 !WE NOW HAVE 5 CLOSE, UNIQUE STATIONS, NEED NO MORE
        
30    CONTINUE
      END DO	!END OF LIST OF STATION PHASES

C--PROGRAM MAY REACH THIS POINT BEFORE J = 5
40    CONTINUE

C--FOR DEBUGGING, LIST OUT THE 5 CLOSEST STATIONS & THEIR ELEVATIONS
C      WRITE (6,*) KYEAR2,KMONTH,KDAY,KHOUR,KMIN
C      DO J=1,JDSTA
C        WRITE (6,1030) STASU(J),DISSU(J),IELSU(J)
C1030    FORMAT (A9,' DIS=',F7.2,' EL=',I5)
C      END DO


C--DEPTH DATAUM IS AVERAGE OF JDSTA (MAX 5) CLOSEST STATIONS
      IDEPDAT=0
      DO J=1,JDSTA
        IDEPDAT=IDEPDAT+IELSU(J)
      END DO
      IDEPDAT=NINT(1.*IDEPDAT/JDSTA)
      
C--GET GEOID DEPTH & SET FLAG OF WHAT REPORTED DEPTH IS
50    ZGEOID=Z1-IDEPDAT/1000.
      IF (IMODG.EQ.1) THEN	!CRE,CRV,CRL MODELS
        CZFLAG='G'
        ZREP=Z1
      ELSE IF (IMODG.EQ.0) THEN	!CRT,CRH MODELS
        IF (LGEOID) THEN
          CZFLAG='G'
          ZREP=ZGEOID	!FOR TRADITIONAL MODELS, MAKE DEPTH DATUM CORRECTION
        ELSE
          CZFLAG='M'
          ZREP=Z1	!FOR HYPOELLIPSE MODELS, DEPTH IS ALREADY TO GEOID
        END IF
      END IF
      
      RETURN
      END

	subroutine sort3 (n,ra,rb,irc)
c--heapsort subroutine from numerical recipes book
C--SORTS ENTRIES OF RA (NUMBERS) AND RB (9 CHARACTERS) IN ASCENDING ORDER
C  OF RA. RA, RB & IRC ARE REPLACED WITH SORTED VALUES. N VALUES ARE SORTED.
	character rb*9, rrb*9
	dimension ra(n), rb(n), irc(n)
	l=n/2+1
	ir=n
10	continue
	if (l.gt.1) then
	  l=l-1
	  rra=ra(l)
	  rrb=rb(l)
	  iirc=irc(l)
	else
	  rra=ra(ir)
	  rrb=rb(ir)
	  iirc=irc(ir)
	  ra(ir)=ra(1)
	  rb(ir)=rb(1)
	  irc(ir)=irc(1)
	  ir=ir-1
	  if (ir.eq.1) then
	    ra(1)=rra
	    rb(1)=rrb
	    irc(1)=iirc
	    return
	  end if
	end if
	i=l
	j=l+l

20	if (j.le.ir) then
	  if (j.lt.ir) then
	    if (ra(j) .lt. ra(j+1)) j=j+1
	  end if
	  if (rra.lt.ra(j)) then
	    ra(i)=ra(j)
	    rb(i)=rb(j)
	    irc(i)=irc(j)
	    i=j
	    j=j+j
	  else
	    j=ir+1
	  end if
	  goto 20
	end if
	ra(i)=rra
	rb(i)=rrb
	irc(i)=iirc
	goto 10
	end

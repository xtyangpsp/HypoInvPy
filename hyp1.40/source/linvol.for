c linvol.for    []
      subroutine linvol(delta, depthsv, stzsv, ak, voa, aha, vi, t,
     *  ain, dtdd, dtdh)
c model:  linearly increasing velocity over a halfspace
c bill gawthrop  (modified by j.c. lahr  3/6/91)
c     delta  epicentral distance (km)
c     depthsv eq depth (km)
c     stzsv   station depth (km)
c     ak     velocity gradient
c     voa    v at surface
c     aha    thickness of linear layer (depth to half space)
c     vi     velocity of half space
c     t      travel time
c     ain    angle of incidence
c     dtdd   partial derivative of tt wrt distance
c     dtdh   partial derivative of tt wrt depth
      integer punt
      common /punt/ punt 
      parameter (pi = 3.14159265)
      parameter (rad = pi/180.)
      parameter (deg = 1.0/rad)
      logical flip
      dist = sqrt(delta**2 + (depthsv - stzsv)**2)
      if (dist .lt. .0001) then
c       distance less than 10 cm
        t = 0.0
        dtdd = 1./(voa + ak*depthsv)
        dtdh = 1./(voa + ak*depthsv)
        ain = 90.
        return
      endif
      if((stzsv .ge. aha - .0001) .and.
     *  (depthsv .ge. aha - .0001)) then
c both eq and station are within (or within 10 cm) of the
c constant velocity halfspace
c       write(punt, *) 'eq and station both in halfspace'
        t = dist/vi
        dtdd = delta/(dist*vi)
        dtdh = (depthsv - stzsv)/(dist*vi)
        if (stzsv .eq. depthsv) then
          ain = 90.0
        else
          if(stzsv .lt. depthsv) then
c           upward
            ain = 180 + deg*atan(delta/(stzsv - depthsv))
          else
c           downward
            ain = deg*atan(delta/(stzsv - depthsv))
          endif
        endif
        return
      endif
      if(depthsv .ge. stzsv) then
        flip = .false.
        deptha = depthsv
        stza = stzsv
      else
        flip = .true.
        deptha = stzsv
        stza = depthsv
      endif
c if near boundary, then move to boundary
      if(abs(deptha - aha) .le. .0001) deptha = aha
      if(abs(stza - aha) .le. .0001) stza = aha
c set vsource
      if(flip) then
        if(stza .le. aha) then
          vsource = voa + stza*ak
        else
          vsource = vi
        endif
      else
        if(deptha .le. aha) then
	  vsource = voa + deptha*ak
	else
	  vsource = vi
	endif
      endif
c create modified model with station at surface
      ah = aha - stza
      vo = voa + ak*stza
      depth = deptha - stza
      stz = 0.0
      if (depth .le. ah) then
c source: within the layer
c       write(punt, *) ' source: within the layer'
        vz = vo + depth*ak
c       compute time of direct wave
        call tatime(delta,depth,vo,ak,ah,tat,aina)
c       compute time of refracted wave
        call tbtime(delta,depth,vo,ak,ah,vi,tbt,ainb)
        if (tbt .lt. tat) then
c         refracted wave is faster
          t = tbt
          if(flip) then
c           downward ray
            anin = vsource*sin(ainb)/vz
            ain = asin(anin)
          else
c           downward ray
            anin = sin(ainb)
            ain = ainb
          endif
        else if (tat .eq. 900000.) then
c         neither direct nor refracted wave could be computed!
          print *,
     *      'neither direct nor refracted wave could be computed!'
          print *, 'so stop'
          stop '1: abort from linvol'
        else
c         direct wave is faster
          t = tat
          if(flip) then
c           downward ray
            anin = vsource*sin(aina)/vz
            ain = asin(anin)
          else
c           upward or downward ray
            anin = sin(aina)
            ain = aina
          endif
        endif
      else
c source: within the halfspace
c       write(punt, *) ' source: within the halfspace'
c	write(punt, *) 'delta,depth,vo,ak,ah,vi'
c	write(punt, *)  delta,depth,vo,ak,ah,vi 
        call tdtime(delta,depth,vo,ak,ah,vi,td,aina)
        if (td.eq.900000.) then
c         direct wave could not be computed!
          print *,
     *      'direct wave could not be computed, so stop!'
          stop '2: abort from linvol'
        else
c         direct wave from halfspace
          t = td
          if(flip) then
c           downward ray
            anin = vsource*sin(aina)/vi
            ain = asin(anin)
          else
c           upward ray
            anin = sin(aina)
            ain = aina
          endif
        endif
      endif
      dtdd = anin/vsource
      dtdh = -cos(ain)/vsource
      ain = deg*ain
      return
      end
c end linvol
c----------------------------
c tatime.for    []
      subroutine tatime (delta,depth,vo,ak,ah,tat,toa1)
c model:  linearly increasing velocity over a halfspace
c source: within the layer
c path:   direct wave
c   bill gawthrop  (modified by j.c. lahr  3/6/91)
c     vo   v at surface
c     ak   velocity gradient
c     ah   thickness of linear layer (depth to half space)
      if (delta.ge. 0.1) then
        a = (vo/ak)**2
        b = (vo/ak + depth)**2
        xc = (delta*delta + a - b)/(2.*delta)
        r = sqrt(xc*xc + b)
        toa1 = acos(xc/r)
        if ((r - vo/ak).gt.ah.and.toa1.le.1.57) then
          tat = 900000.
          return
        endif
        c = r*r + delta*xc - xc*xc
        tat = alog((c + r*delta)/(c - r*delta))/(2.*ak)
      else
        tat = alog((vo + depth*ak)/(vo))/ak
        toa1 = 3.1415926
      endif
      return
      end
c end tatime
c-------------------
c tbtime.for    []
      subroutine tbtime (delta,depth,vo,ak,ah,vi,tbt,toa2)
c model:  linearly increasing velocity over a halfspace
c source: within the layer
c path:   refracted wave
c   bill gawthrop  (modified by j.c. lahr  3/6/91)
c     vo   v at surface
c     ak   velocity gradient
c     ah   thickness of linear layer (depth to half space)
c     vi   velocity of half space
      character*55 mess(3)
      integer punt
      common /punt/ punt
      common /logfil/ logfil
      a = (depth + vo/ak)**2
      b = (ah + vo/ak)**2
      c = (vo/ak)**2
      r = vi/ak
      r2 = r*r
      if(b .le. r2) then
        xc1 = sqrt(r2 - a)
        xi = xc1 - sqrt(r2 - b)
        xc2 = sqrt(r2 - c)
        x2 = xc2 - xc1 + xi
        if ((delta - xi - x2) .gt. 0.0) then
          d = r2 + xi*xc1 - xc1*xc1
          e = r2 + x2*xc2 - xc2*xc2
          tbt = alog(
     *     (d + r*xi)*(e + r*x2)/((d - r*xi)*(e - r*x2))
     *          )/(2.*ak) + (delta - xi - x2)/vi
          toa2 = acos(xc1/r)
        else
          tbt = 900000.
        endif
      else
        mess(1) =  ' for a linear increase over halfspace model the'
        mess(2) =  ' velocity at the base of the linear model is '
        mess(3) =  ' greater than the half space velocity, so stop.'
        write(punt, '(a, /, a, /, a)') mess
        write(logfil, '(a, /, a, /, a)') mess
        stop 'abort from tbtime'
      endif
      return
      end
c end tbtime
c--------------------
c tdtime.for    []
      subroutine tdtime (delta, depth, vo, ak, ah, vi, tdt, toa1)
c
c model:  linearly increasing velocity over a halfspace
c source: within the halfspace
c path:   direct wave
c bill gawthrop  (modified by j.c. lahr  3/6/91)
c
c   this sub determines the traveltime of a ray from an origin
c   below the moho in constant velocity material passing thru a crust
c   with velocity increasing linearly with depth
c   bill gawthrop
c     xo   delta
c     yo   depth
c     vo   v at surface
c     ak   velocity gradient
c     ah   thickness of linear layer (depth to half space)
c     vi   velocity of half space
c     tdt
c     toa1
      double precision xp, xp2, xp3, a, b, c, d, xt1, xt2, el
      double precision xpold, xt, xr, chx, denom, xpa, xpb
      double precision r2, r
      integer punt
      common /punt/ punt
c     write(punt, *) 'in tdtime: delta, depth, vo, ak, ah, vi'
c     write(punt, *)  delta, depth, vo, ak, ah, vi 
      if ((delta .lt. .001) .or. (delta/depth .lt. .001)) then
        tdt = (depth - ah)/vi + alog((vo + ak*ah)/(vo))/ak
        toa1 = 3.1415926
c	write(punt, *) 'tdt, toa1 ', tdt, toa1
        return
      endif
      a = (vi/ak)**2
      b = a - (vo/ak)**2
      c = a - (ah + vo/ak)**2
      d = a*(depth - ah)**2
ctest
      iwrt = 0
      chx = 0.d0
    6 continue
ctest
      xp = .8d0*delta
      do 15 i = 1, 25
        xpold = xp
        xp2 = xp*xp
        xp3 = xp2*xp
        xt1 = dsqrt(d/xp2 + b)
        xt2 = dsqrt(d/xp2 + c)
        xt = xp + xt1 - xt2
        xr = xt - delta
ctest
        if (iwrt .eq. 1) write(punt, 100)
     *    i, delta, depth, xr, xpold, chx,
     *    xp, xt1, xt2, xp3
100     format(i5, 8f13.7, f13.2)
ctest
        if (dabs(xr).lt..0001d0) then
          xc = delta - xt1
          r2 = a + d/xp2
          r = dsqrt(r2)
          ta = r2 - delta*xp + delta*xc + xp*xc - xc*xc
          tb = r*delta - r*xp
          el = dsqrt((xp2 + (depth - ah)**2)*1.d0)
c	  write(punt, *) 'xc, r2, r, ta, tb, el ', xc, r2, r, ta, tb, el
          tdt = alog((ta + tb)/(ta - tb))/(2.*ak) + el/vi
          toa1 = 1.5707963 + acos(xp/el)
c	  write(punt, *) 'tdt toa1 ', tdt, toa1
c
          return
        endif
        denom = xp3*xt1
        if (denom .eq. 0.d0) then
          xp = -1.d0
        else
          xpa = d/(xp3*xt1)
          xpb = d/(xp3*xt2)
          denom = 1.d0 - xpa + xpb
          if (denom .eq. 0.) then
            xp = -1.d0
          else
            chx = xr/(1.d0 - xpa + xpb)
            xp = xp - chx
          endif
        endif
   10   continue
        if (xp .lt. 0.000001) xp = xpold/4.
   15 continue
ctest
      tdt = 900000.
      if(iwrt .eq. 1) return
      write(punt, 200) delta, depth, vo, ak, ah, vi, tdt, toa1,
     * xr, xpold, chx, xp, xt1, xt2, xp3
  200 format(' sub. tdtime convergence problem', /,
     * '         delta           z          vo           ak
     *             ah            vi         tdt         toa1', /,
     * 5x, 8f13.7, /,
     *  '    xr          xpold          chx          ',
     *  ' xp          xt1          xt2          xp3', /, 5x, 7f13.7)
      iwrt = 1
      goto 6
ctest
      end
c end tdtime

c linv.for    []
c test driver for subroutine linv
c      data elevmx/3.0/, ielest/2000/, x/10./, zeq/3.0/
c      data v/3.0/, ak/1.0/
c      print *, 'welcome to linv!  this program tests subroutine'
c      print *, 'linv, which computes travel times in a half space'
c      print *, 'with linearly increasing velocity.'
c
c20    elevmx = raskk(
c     *'maximum topographic elevation (reference elevation in km)',
c     *          elevmx)
c      ielest = iaskk(
c     *'station elevation above sea level', ielest)
c      stz = elevmx - ielest*.001
c      v = raskk('velocity at reference elevation (km/s)', v)
c      ak = raskk('velocity gradient (km/s per km)', ak)
c30    x = raskk('epicentral distance (km)', x)
c      zeq = raskk('depth of eq beneath reference elevation (km)', zeq)
c
c      call linv(x, zeq, v, ak, t, ain, dtdd, dtdh, stz)
c      vsta = v + ak*(elevmx - ielest/1000.)
c      zeff = zeq - (elevmx - ielest/1000.)
c      call linvjl(x, zeff, vsta, ak, t, ain, dtdd, dtdh, stz,
c    *  vst, veq)
c
c      print *, 'linvjl results to test linv when eq is deeper',
c     * ' than station:'
c      print *, 'travel time (s) = ', t
c      ain =  asin(c      if(ain .lt. 0.) ain = 180. + ain
c      ain = 180. - ain
c      print *, 'ain in degrees = ', ain
c      print *, 'partial of tt wrt dist and depth = ', dtdd, dtdh
c
c      if(x .gt. 0.) goto 20
c      stop
c      end
      subroutine linv(x, zeq, v, ak, t, ain, dtdd, dtdh, stz,
     *  vst, veq)
c compute travel time and partial derivatives for a linear increasing
c velocity model in a region of significant topography.
      real x
c          x      epicentral distance in km
      real zeq
c          zeq    depth in km of earthquake beneath reference elevation
      real v
c          v      velocity at reference elevation in km/sec
      real ak
c          ak     velocity gradient km/sec per km
      real t
c          t      computed travel time in sec
      real ain
c          ain    angle-of-incidence
      real dtdd
c          dtdd   partial derivative of tt wrt distance
      real dtdh
c          dtdh   partial derivative of tt wrt depth
      real stz
c          stz    depth in km of station beneath reference elevation
      real vst
c          vst    velocity at station
      real veq
c          veq    velocity at earthquake
      real dz
c          dz     depth of eq beneath station
      real dc
c          dc     distance in km from eq or station, whichever is
c                 at higher elevation, to center of raypath circle
      parameter (pi = 3.14159265)
      parameter (rad = pi/180.)
      parameter (deg = 1.0/rad)
c      print *, 'station elevation (m)    = ', ielest
c      print *, 'eq depth (km)            = ', zeq
c      print *, 'epicentral distance      = ', x
c      print *, 'velocity at ref elev and grad = ', v, ak
      vst = v + ak*stz
      veq = v + ak*zeq
      dz = (zeq - stz)
      dist = sqrt(dz*dz + x*x)
      if (dist .lt. .0001) then
c       distance less than 10 cm
        t = 0.0
        dtdd = 1./veq
        dtdh = 1./veq
        anin = 1.
        ain = 90.
        return
      endif
      p = ak*ak*(x*x + dz*dz)/(2.*vst*veq) + 1.
      t = coshi(p)/ak
      if (dz .lt. 0) then
c earthqauke is above station (ray is down-going)
c        print *, 'eq is above station, dz = ', dz
        if (x .eq. 0) then
          anin = 0.0
          ain = 0.0
        else
          dc = (x*x + 2.*abs(dz)*veq/ak + dz*dz)/(2.*x)
          anin = (veq/ak) / sqrt(dc*dc + (veq/ak)**2)
          ain = deg*asin(anin)
c          print *, 'anin = ', anin
        endif
      else
c station is above earthquake (ray may be up- or down-going)
c        print *, 'station is above eq, dz = ', dz
        alam = vst*(cosh(ak*t) - 1.)/ak
        if(ak .le. .0001) alam = 0.0
        dif = dz - alam
        anin = x/sqrt(dif*dif + x*x)
      endif
      dtdd = x*ak/( vst*(vst + ak*dz)*sqrt(p*p - 1.) )
      dtdh = (dz*ak/(vst*(vst+ak*dz)) -
     *         ak*ak*(x*x + dz*dz)/(2.*vst*((vst + ak*dz)**2.)))/
     *         ((p*p - 1.)**0.5)
      if(dtdh .gt. 0.0) then
c       upgoing
        ain = 180. - deg*asin(anin)
      else
c       downgoing
        ain = deg*asin(anin)
      endif
      return
      end
c end linv
c coshi.for    []
      function coshi(x)
      p = x + sqrt(x*x - 1.)
      coshi = alog(p)
      return
      end
c end coshi

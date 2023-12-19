c trvcon.for    []
      subroutine trvcon(delta, zsv, t, ain, dtdd, dtdh,
     *  lbeg, lend, lbegm1, nlayers,
     *  ivlr, ivl, thk, nlay, ldx, wti, wtk,
     *  tid, did, jref, vsq, vsqd, v, vi, f, vs, vt, msta, stzsv)
c compute travel time and derivatives for model with constant
c    velocity layers
c indecies for this crustal structure, model number imod.  (lbegm1 = lbeg-1)
c
c    interface    index            velocity       thickness
c
c----surface------(lbeg)--------------------------------------(lstm1)---
c                                  v(lbeg)        thk(lbeg)
c
c-----------------(lbeg+1)------------------------------------(lst)-----
c                                    station is in this layer
c                                  v(lbeg+1)=vt   thk(lbeg+1)
c
c-----------------(lbeg+2) = (leqm1)--------------------------(lstp1)---
c                                  v(leqmq)       thk(leqm1)
c
c-----------------(lbeg+3) = (leq)--------------------------------------
c
c    epicenter is in this layer    v(leq)=vs       thk(leq)
c
c-----------------(lbeg+4) = (leqp1)------------------------------------
c                                  v(jj)          thk(jj)
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c-----------------(lend-1) = -------------------------------------------
c
c-top of half space-(lend)----------------------------------------------
c                                  v(lend)        infinite thickness
c
c
C FWK UNNEEDED PARAMETERS COMMENTED OUT 1/29/2011
C OTHERWISE CALCULATION CODE UNCHANGED
C      include 'params.inc'
c     HYPOELLIPSE params.inc for the SUN version
c
c     (npa - 1) is the maximum number of phases (P plus S plus S-P) for
c               a single earthquake.
C FWK      parameter (npa = 1501)

c     nsn is the maximum number of stations that can be included in the
c               station list.
c     parameter (nsn = 1024)
C FWK      parameter (nsn = 2048)

c     lmax is the maximum number of VELOCITY records allowed to define
c               the velocity models.
c     parameter (lmax = 96)
c      parameter (lmax = 128)
C FWK       parameter (lmax = 160)
      PARAMETER (LMAX=20)	!SHOULD BE THE MAX # LAYERS PER MODEL NLYR=20

c     mmax is the maximum number of velocity models that can be specified.
c     remember that a model with variable vp/vs ratio counts as two.
C FWK      parameter (mmax = 25)
C FWK      parameter (mmaxp3 = mmax+3)

c     lmmax is the maximum number of records that can be used to define
c     a single velocity model.  A model defined by N records will have
c     N-1 layers over a half space.
      parameter (lmmax = 20) 

c     lmax2 and narray are other array dimensions computed from lmax
C FWK      parameter (lmax2 = lmax + 2)
C FWK      parameter (narray = 2*lmmax*lmax)


      parameter (pi = 3.14159265)
      parameter (rad = pi/180.)
      parameter (deg = 1.0/rad)
      dimension jref(lmax)
      character*5 msta
      double precision umax, utemp, div, u1, u2, tkj, tkjsq, det, c1,
     *  c2, a, b, ac4, dmb, u0, vsqu, xtest, dlimit, d0, d1, d2, 
     *  ssqu
c     double precision top, bot
      dimension vsq(lmax), vsqu(lmax), vsqd(lmmax, lmax), v(lmax+2)
      dimension vi(lmax), thk(lmax+1)
      dimension tid(lmax, lmmax), did(lmax, lmmax), f(lmax, lmmax)
      dimension didj(lmmax), tidj(lmmax), tr(lmmax), div(lmax)
      integer punt
      DATA PUNT /6/
C FWK      common /punt/ punt
      logical flip
      idrt = 0
      tmin = 99999.
      z = zsv
      stz = stzsv
      if(abs(z - stz) .lt. 0.001) then
        stz = z
      endif
c
c decide which is deeper, station or eq
c and if station is deeper, then reverse roles for now so that
c "eq" is always deeper
      if (z .ge. stz) then
        flip = .false.
      else
        flip = .true.
        temp = z
        z = stz
        stz = temp
      endif
cd    print *, 'z = ', z, ', stz = ', stz
c determine which layer contains the hypocenter.
      sum = 0.0
      do 19 leq = lbeg, lend - 1
        sum = sum + thk(leq)
        if (z .le. sum) then
c         the hypocenter is not in the halfspace
c         leq is the layer containing the hypocenter, leqp1 is the one below,
c         and leqm1 is the one above.
          leqp1 = leq + 1
          leqm1 = leq - 1
          leqr = leq - lbegm1
c         tkj is the distance to the boundary above the eq
          tkj = z - sum + thk(leq)
          tkjsq = tkj**2
cd        print *, 'eq is in layer # ', leq, ', ', tkj, ' km below top'
          goto 23
        endif
19    continue
c the hypocenter is in the halfspace
      leq = lend
cd    print *, 'eq is in halfspace, layer ', leq
      leqm1 = leq - 1
      leqr = leq - lbegm1
      tkj = z - sum
      tkjsq = tkj**2
c set weight if this ray was supposed to be refracted
      if (nlay .ne. 0) wti = 0.0
      idrt = 1
c determine which layer contains the station
23    sum = 0.0
      do 24 lst = lbeg, lend - 1
        sum = sum + thk(lst)
        if (stz .le. sum) then
c         the station is not in the halfspace
c         lst is the layer containing the station, lstp1 is the one below,
c         and lstm1 is the one above.
          lstp1 = lst + 1
          lstm1 = lst - 1
          lstr = lst - lbegm1
c         tkjst is the distance to the boundary above the station
          tkjst = stz - sum + thk(lst)
          tkjstsq = tkj**2
cd        print *,
cd   *      'sta is in layer # ', lst, ', ', tkjst, ' km below top'
          goto 25
        endif
24    continue
c the station is in the halfspace
      lst = lend
cd    print *, 'station is in halfspace, layer ', lst
      lstm1 = lst - 1
      lstr = lst - lbegm1
      tkjst = stz - sum
      tkjstsq = tkj**2
      idrt = 1
      if (leq .eq. lst) then
c eq and station in halfspace
        vs = v(leq)
        vt = vs
        goto 43
      endif
c set the velocity at the source to vs and receiver to vt
25    if(flip) then
        vs = v(lst)
        vt = v(leq)
      else
        vs = v(leq)
        vt = v(lst)
      endif
      tjljll = tid(leq, leqr)
      djljll = did(leq, leqr)
c     correct for variable layers
      if ((ivlr .gt. 0) .and. (leq .gt. ivl) .and.
     *    (jref(leq) .ne. 0)) then
        tjljll = tjljll + thk(ivl)*vsqd(leqr, ivl)*vi(leq)
        djljll = djljll + thk(ivl)/vsqd(leqr, ivl)
        if (leq .gt. ivl + 1) then
          tjljll = tjljll + thk(ivl+1)*vsqd(leqr, ivl+1)*vi(leq)
          djljll = djljll + thk(ivl+1)/vsqd(leqr, ivl+1)
        endif
      endif
c     correct for portion of layer above station
      if ((leq .gt. lst) .and. (vsqd(leqr, lst) .ne. 0.0)) then
        tjljll = tjljll - tkjst*vsqd(leqr, lst)*vi(leq)
cd      print *, 'vsqd(leqr, lst), leqr, lst ', vsqd(leqr, lst), leqr,
cd   *    lst
        djljll = djljll - tkjst/vsqd(leqr, lst)
      endif
c     correct for full layers above station
      if(lst .gt. lbeg) then
        do 30 layer = lbeg, lst - 1
          if (vsqd(leqr, layer) .eq. 0.0) goto 30
          tjljll = tjljll - thk(layer)*vsqd(leqr, layer)*vi(leq)
          djljll = djljll - thk(layer)/vsqd(leqr, layer)
30      continue
      endif
      if ((idrt .eq. 1) .or. (delta .lt. .0001)) goto 43
c***************************************************************************
c calculations for refracted waves.
c***************************************************************************
c if this arrival is to be refracted along a certain boundary
c   calculate travel time or set weight = 0 if this refracted wave
c   can not exist.
      if (nlay .ne. 0) then
        k = nlay + 1
c       k at non-refraction boundary
        if (jref(lbeg + nlay) .eq. 0) goto 31
c       k below(?) halfspace or k above hypocenter layer
        if ((k .gt. nlayers) .or. (k .lt. leqr)) goto 31
c       refraction in layer of hypocenter (compute for direct wave)
        if (k .eq. leqr) then
          idrt = 1
          goto 43
        endif
        didj(k) = did(leq, k)-tkj/vsqd(k, leq)
        if (ivlr .gt. 0) then
          if (k .gt. ivlr) didj(k) = didj(k) +
     *      f(ivl, leqr)*thk(ivl)/vsqd(k, ivl)
          if (k .gt. ivlr+1)
     *      didj(k) = didj(k) + f(ivl+1, leqr)*thk(ivl+1)/
     *      vsqd(k, ivl+1)
        endif
        if (didj(k) .gt. delta) goto 31
        kk = nlay + lbegm1 + 1
        tmin = tid(leq, k )-tkj*vsqd(k, leq)*vi(kk) + delta*vi(kk)
        if (ivlr .le. 0) goto 42
        if (k .gt. ivlr)
     *    tmin = tmin+f(ivl, leqr)*thk(ivl)*vsqd(k, ivl)*vi(k+lbegm1)
        if (k .gt. ivlr+1) tmin = tmin + f(ivl+1, leqr)*thk(ivl+1)*
     *    vsqd(k, ivl+1)*vi(k+lbegm1)
        goto 42
c       if not possible, set p & s weights to zero.
31      wti = 0.0
        if (ldx .ne. 0) wtk = 0.0
      endif
c calculate intercept times and critical distances for all possible
c refracted waves.
      do 34 m = leqp1, lend
        mm = m - lbegm1
        tidj(mm) = 0.
        didj(mm) = 0.
        if (jref(m) .eq. 0) goto 34
        tidj(mm) = tid(leq, mm) - tkj*vsqd(mm, leq)*vi(m)
        didj(mm) = did(leq, mm) - tkj/vsqd(mm, leq)
c       correct for portion of layer above station
        tidj(mm) = tidj(mm) - tkjst*vsqd(mm, lst)*vi(m)
        didj(mm) = didj(mm) - tkjst/vsqd(mm, lst)
c       correct for full layers above station
        if(lst .gt. lbeg) then
          do 32 layer = lbeg, lst - 1
            tidj(mm) = tidj(mm) - thk(layer)*vsqd(mm, layer)*vi(m)
            didj(mm) = didj(mm) - thk(layer)/vsqd(mm, layer)
32        continue
        endif
        if ((ivlr .le. 0) .or. (m .le. ivl)) goto 34
        xtemp = f(ivl, leqr)*thk(ivl)
        tidj(mm) = tidj(mm)+xtemp*vsqd(mm, ivl)*vi(m)
        didj(mm) = didj(mm)+xtemp/vsqd(mm, ivl)
        if (m .le. ivl+1) goto 34
        xtemp = f(ivl+1, leqr)*thk(ivl+1)
        tidj(mm) = tidj(mm) + xtemp*vsqd(mm, ivl+1)*vi(m)
        didj(mm) = didj(mm) + xtemp/vsqd(mm, ivl+1)
34    continue
c xovmax is the maximum cross over distance.  at this or greater
c distance, refracted wave must be first.
      xovmax = (tidj(leqr + 1)-tjljll)/(vi(leq)-vi(leqp1))
cd    print *, 'xovmax = ', xovmax
      if (jref(leqp1) .eq. 0 .or. jref(leq) .eq. 0) xovmax = 9999.
c calculate all possible refracted wave travel-times and find least
      do 40 m = leqr + 1, nlayers
        if (delta .lt. didj(m)) goto 40
        if (jref(m+lbegm1) .eq. 0) goto 40
        tr(m) = tidj(m)+delta*vi(m+lbegm1)
        if (tr(m) .ge. tmin) goto 40
        k = m
        kk = k + lbegm1
        tmin = tr(m)
40    continue
cd    print *, 'min refracted tt and layer = ', tmin, k
c     if (delta .ge. xovmax) then
      if (delta .lt. xovmax) goto 43
c refracted wave is fastest
c relabel travel-time and find derivatives
42      t = tmin
cd      print *, 'refracted wave is fastest'
        dtdd = vi(kk)
        if(flip) then
          dtdh = -vsqd(k, lst)*dtdd
          anin = v(lst)*vi(kk)
        else
          dtdh = -vsqd(k, leq)*dtdd
          anin = v(leq)*vi(kk)
        endif
        ain = deg*asin(anin)
        return
c     endif
c*************************************************************************
c calculations for the direct waves
c*************************************************************************
43    if (leq .eq. lst) then
c for hypocenter and station in same layer
cd      print *, 'hypocenter & station in same layer'
        sqt = sqrt((stz - z)**2 + delta**2)
        tdj1 = sqt/v(leq)
        if (tdj1 .ge. tmin) goto 42
        t = tdj1
        if (sqt .lt. .0001) then
c         hypocentral distance (sqt) less than 0.1 meters
          dtdd = vi(leq)
          dtdh = vi(leq)
          ain = 90.
        else
          dtdd = delta*vi(leq)/sqt
          anin = delta/(sqt)
          if(flip) then
c           downgoing
            dtdh = -(z - stz)*vi(leq)/sqt
            ain = deg*asin(anin)
cd          print *, 'downgoing ain = ', ain
          else
c           upgoing
            dtdh = (z - stz)*vi(leq)/sqt
            ain = 180. - deg*asin(anin)
cd          print *, 'upgoing ain = ', ain
          endif
        endif
        return
      endif
c find parameters of direct wave from hypocenter not in station layer.
c invert the funct. delta(u) (where u is the sine of the angle
c    of incidence) by iterating to find u.
c    use the interpolating funct.:
c       d(u) = c/(umax - u) + a*(u - umax) + b
c    to find better values of u0 which satisfy
c    delta(u0) = given distance within some error, say .03 km.
c use a linear interpolating finction for large epicenteral
c    distances to save computing time.
c----
c find umax for this layer and velocity structure.
c (umax is less than 1 for a hypocenter in a low velocity layer)
      if (delta .lt. 0.0001) then
c vertical ray path
cd      print *, 'vertical ray path'
        u0 = 0.0
        ssqu = 1.d0
        umax = 1.0
        vsqu(lst) = (thk(lst)-tkjst)/(v(leq)/v(lst))
        if(lst+1 .le. leqm1) then
          do 52 l = lst + 1, leqm1
            vsqu(l) = thk(l)/(v(leq)/v(l))
52        continue
        endif
        goto 220
      endif
c non-vertical ray path
cd    print *, 'non-vertical ray path'
      umax = v(leq)/v(lst)
      do 54 l = lst, leqm1
        utemp = v(leq)/v(l)
        div(l) = utemp**2
        if (utemp .lt. umax) umax = utemp
54    continue
cd    print *, 'non-vertical ray.  umax = ', umax
      if (umax .gt. 1.d0) umax = 1.d0
c choose 2 initail points (u1, d1) and (u2, d2) for interpolation,
c    depending on epicentral distance.
      u1 = .75d0*umax
      u2 = delta/sqrt((z-stz)**2+delta**2)
      if (u2 .gt. umax-1.d-4) u2 = umax-1.d-4
      if (u2 .gt. u1) goto 56
      utemp = u1
      u1 = u2
      u2 = utemp
56    continue
c first take care of eq layer and station layer:
      d1 = tkj*u1/dsqrt(1.d0-u1**2) +
     *     (thk(lst) - tkjst)*u1/dsqrt(  vsq(leq)/vsq(lst) - u1**2 )
c    *     (thk(lst) - tkjst)*u1/dsqrt( (v(leq)/v(lst))**2 - u1**2 )
      d2 = tkj*u2/dsqrt(1.-u2**2) +
     *     (thk(lst) - tkjst)*u2/dsqrt(  vsq(leq)/vsq(lst) - u2**2 )
c    *     (thk(lst) - tkjst)*u2/dsqrt( (v(leq)/v(lst))**2 - u2**2 )
c then add terms for layers in between:
      if(lst+1 .le. leqm1) then
        do 58 l = lst+1, leqm1
          d1 = d1 + u1*thk(l)/dsqrt(  vsq(leq)/vsq(l) - u1**2 )
c         d1 = d1 + u1*thk(l)/dsqrt( (v(leq)/v(l))**2 - u1**2 )
          d2 = d2 + u2*thk(l)/dsqrt(  vsq(leq)/vsq(l) - u2**2 )
c         d2 = d2 + u2*thk(l)/dsqrt( (v(leq)/v(l))**2 - u2**2 )
58      continue
      endif
c begin the iteration loop.
      do 150 ll = 1, 25
cd      print *, 'u1, d1, u2, d2 ', u1, d1, u2, d2
c now find the constants a, b, c of d(u) subject to the conditions
c       d(u1) = d1, d(u2) = d2, and d(0) = 0
        det = u1*u2*(u2-u1)
        c1 = d1*u2*(umax-u1)
        c2 = d2*u1*(umax-u2)
        a = (c1-c2)/det
        b = (c1*(2.d0*umax-u2)-c2*(2.d0*umax-u1))/det
        ac4 = 4.d0*a*(a*umax**2-b*umax)
        dmb = delta - b
c invert d(u0) for u0 and find d0 = delta(u0)
        u0 = umax+(dmb-dsqrt(dmb**2+ac4))/(2.d0*a)
        ssqu = dsqrt(1. - u0**2)
        if (u0 .ge. .99999d0) then
c if wave leaves hypocenter nearly horizontally (ie u goes to 1)
c    then compute tdir as wave travelling along top of layer leq
cd        print *, 'compute tdir as wave along top of eq layer'
          tdir = tjljll + delta*vi(leq)
          if (tdir .lt. tmin) then
c           travel time derivatives for hypocenter at top of a layer
            t = tdir
            if(flip) then
c             downgoing
              anin = v(lst)*vi(leq)
              ain =  deg*asin(anin)
              dtdh = -sqrt(1. - anin**2)*vi(lst)
              dtdd = vi(leq)
            else
c             upgoing
              ain = 90.
              dtdh = 0.0
              dtdd = vi(leq)
            endif
            return
          endif
          goto 42
        endif
        if (u0 .ge. umax-1.d-6) then
          u0 = umax-1.d-6
          ssqu = dsqrt(1. - u0**2)
        endif
        d0 = tkj*u0/ssqu
        vsqu(lst) = (thk(lst) - tkjst)/dsqrt( vsq(leq)/vsq(lst)-u0**2)
c       vsqu(lst) = (thk(lst) - tkjst)/dsqrt((v(leq)/v(lst))**2-u0**2)
        d0 = d0 + u0*vsqu(lst)
        if(lst+1 .le. leqm1) then
          do 112 l = lst+1, leqm1
            vsqu(l) = thk(l)/dsqrt( vsq(leq)/vsq(l)-u0**2)
c           vsqu(l) = thk(l)/dsqrt((v(leq)/v(l))**2-u0**2)
            d0 = d0 + u0*vsqu(l)
  112     continue
        endif
c set the distance accuracy of iteration.
        dlimit = .04d0
        if (u0 .lt. .05d0 ) goto 113
        dlimit = .015d0*v(leq)/u0
c test to see if new ray is within required distance accuracy.
  113   xtest = dabs(delta-d0)
        if (xtest  .lt.  dlimit) then
cd        print *, 'trvcon converged direct ray in ', ll, ' iterations'
          goto 220
        endif
c prepare to iterate again
c replace (u1, d1) or (u2, d2) by the new point (u0, d0)
        if (delta .lt. d2) goto 114
        d1 = d2
        u1 = u2
        d2 = d0
        u2 = u0
        goto 150
  114   if (delta .gt. d1) goto 116
        d2 = d1
        u2 = u1
        d1 = d0
        u1 = u0
        goto 150
  116   if (delta .lt. d0) goto 118
        d1 = d0
        u1 = u0
        goto 150
  118   d2 = d0
        u2 = u0
  150 continue
c if solution of direct wave does not converge in 25 iterations
      uu = asin(u0)*deg
      if (uu .lt. 0.0) uu = uu + 180.
      uu = 180. - uu
      write(punt, 172) msta, delta, z, uu, xtest, dlimit
  172 format(' travel-time solution for a direct wave from station ',
     *  a5, ' at distance ', f6.2, ' from depth ', f6.2,
     *  / ' emerging at an angle of ', f10.4,
     *  ' iterated to the limit.', /,
     *  ' came within ', f10.4, ' km.  xlimit = ', f10.4/)
c travel time & derivatives for direct wave below first layer.
  220 continue
      tdir = tkj*vi(leq)/ssqu
      do 240 l = lst, leqm1
        tdir = tdir + v(leq)*vsqu(l)/vsq(l)
240   continue
cd    print *, 'compare tdir and tmin ', tdir, tmin
      if (tdir .ge. tmin) goto 42
c     if refracted wave is faster than direct wave goto 42
c direct wave is faster
cd    print *, 'direct wave is faster than refracted'
      t = tdir
cd    print *, 'umax, u0 ', umax, u0
      if (umax-u0 .lt. .00001) then
cd      print *, 'umax - u0 .lt. .00001 '
cd      print *,
cd   *    'asin(umax), asin(u0) ', deg*asin(umax), asin(u0)
c       horizontal ray at lower limit
        if(flip) then
c         downgoing
          anin = v(lst)*vi(leq)
          ain = deg*asin(anin)
          dtdh = -sqrt(1. - anin**2)*vi(lst)
          dtdd = vi(leq)
        else
c         upgoing
          ain = 90.
          dtdh = 0.0
          dtdd = vi(leq)
        endif
        return
      endif
c compute partial derivative wrt distance (is independent of flip)
c     bot = 0.d0
c     top = 0.d0
c     do 250 l = lst, leqm1
c       term = (vsq(leq)*vsqu(l)**3)/(vsq(l)*thk(l)**2)
c       bot = bot + term
c       top = top + term*u0/v(leq)
c50   continue
c     dtdd1 = top/bot
c     simplify the formula to:
      dtdd = abs(u0)*vi(leq)
cd    print *, 'flip = ', flip
      if(flip) then
c       downgoing
cd	print *, 'downgoing, lst = ', lst
        anin = abs(u0*v(lst)*vi(leq))
        ain = deg*asin(anin)
        ssqu = sqrt(1.0 - anin**2)
c       dtdh1 = (v(lst)*anin*dtdd - 1.0)/(ssqu*v(lst))
c       simplify the formula to:
        dtdh = -sqrt(1. - anin**2)*vi(lst)
c	print *, 'dtdh fancy, plain = ', dtdh1, dtdh
      else
c       upgoing
cd	print *, 'upgoing, leq = ', leq
        anin = abs(u0)
        ain = 180. - deg*asin(anin)
c       dtdh1 = (1.0 - v(leq)*u0*dtdd)/(ssqu*v(leq))
c       simplify the formula to:
        dtdh = sqrt(1. - anin**2)*vi(leq)
c	print *, 'dtdh fancy, plain = ', dtdh1, dtdh
      endif
      return
      end
c end trvcon

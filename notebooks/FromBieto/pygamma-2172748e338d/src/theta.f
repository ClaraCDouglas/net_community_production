      double precision function theta(s,t0,p0,pr)
c
c ***********************************
c to compute local potential temperature at pr
c using bryden 1973 polynomial for adiabatic lapse rate
c and runge-kutta 4-th order integration algorithm.
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c fofonoff,n.,1977,deep-sea res.,24,489-491
c units:
c       pressure        p0       decibars
c       temperature     t0       deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       reference prs   pr       decibars
c       potential tmp.  theta    deg celsius
c checkvalue: theta= 36.89073 c,s=40 (ipss-78),t0=40 deg c,
c p0=10000 decibars,pr=0 decibars
c
c      set-up intermediate temperature and pressure variables
c
      implicit double precision (a-h,o-z)
      double precision h
c      real p,p0,t,t0,h,pr,xk,s,q,theta,atg
      p=p0
      t=t0
c**************
      h = pr - p
      xk = h*atg(s,t,p)
      t = t + 0.5*xk
      q = xk
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      t = t + 0.29289322*(xk-q)
      q = 0.58578644*xk + 0.121320344*q
      xk = h*atg(s,t,p)
      t = t + 1.707106781*(xk-q)
      q = 3.414213562*xk - 4.121320344*q
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      theta = t + (xk-2.0*q)/6.0
      return
      end

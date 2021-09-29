      double precision function svan(s,t,p0,sigma)
c
c  modified rcm
c ******************************************************
c specific volume anomaly (steric anomaly) based on 1980 equation
c of state for seawater and 1978 practical salinity scale.
c references
c millero, et al (1980) deep-sea res.,27a,255-264
c millero and poisson 1981,deep-sea res.,28a pp 625-629.
c both above references are also found in unesco report 38 (1981)
c units:
c       pressure        p0       decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       spec. vol. ana. svan     m**3/kg *1.0e-8
c       density ana.    sigma    kg/m**3
c ******************************************************************
c check value: svan=981.3021 e-8 m**3/kg.  for s = 40 (ipss-78) ,
c t = 40 deg c, p0= 10000 decibars.
c check value: sigma = 59.82037  kg/m**3 for s = 40 (ipss-78) ,
c t = 40 deg c, p0= 10000 decibars.
c *******************************************************
c
      implicit double precision (a-h,o-z)
      double precision k,k0,kw,k35
c      real p,t,s,sig,sr,r1,r2,r3,r4
c      real a,b,c,d,e,a1,b1,aw,bw,k,k0,kw,k35
c      real svan,p0,sigma,r3500,dr350,v350p,sva,dk,gam,pk,dr35p
c      real dvan
c equiv
      equivalence (e,d,b1),(bw,b,r3),(c,a1,r2)
      equivalence (aw,a,r1),(kw,k0,k)
c ********************
c data
      data r3500,r4/1028.1063,4.8314d-4/
      data dr350/28.106331/
c   r4 is refered to as  c  in millero and poisson 1981
c convert pressure to bars and take square root salinity.
      p=p0/10.
      sr = sqrt(abs(s))
c *********************************************************
c pure water density at atmospheric pressure
c   bigg p.h.,(1967) br. j. applied physics 8 pp 521-537.
c
      r1 = ((((6.536332d-9*t-1.120083d-6)*t+1.001685d-4)*t
     x-9.095290d-3)*t+6.793952d-2)*t-28.263737
c seawater density atm press.
c  coefficients involving salinity
c  r2 = a   in notation of millero and poisson 1981
      r2 = (((5.3875d-9*t-8.2467d-7)*t+7.6438d-5)*t-4.0899d-3)*t
     x+8.24493d-1
c  r3 = b  in notation of millero and poisson 1981
      r3 = (-1.6546d-6*t+1.0227d-4)*t-5.72466d-3
c  international one-atmosphere equation of state of seawater
      sig = (r4*s + r3*sr + r2)*s + r1
c specific volume at atmospheric pressure
      v350p = 1.0/r3500
      sva = -sig*v350p/(r3500+sig)
      sigma=sig+dr350
c  scale specific vol. anamoly to normally reported units
      svan=sva*1.0d+8
      if(p.eq.0.0) return
c ******************************************************************
c ******  new high pressure equation of state for seawater ********
c ******************************************************************
c        millero, et al , 1980 dsr 27a, pp 255-264
c               constant notation follows article
c********************************************************
c compute compression terms
      e = (9.1697d-10*t+2.0816d-8)*t-9.9348d-7
      bw = (5.2787d-8*t-6.12293d-6)*t+3.47718d-5
      b = bw + e*s
c
      d = 1.91075d-4
      c = (-1.6078d-6*t-1.0981d-5)*t+2.2838d-3
      aw = ((-5.77905d-7*t+1.16092d-4)*t+1.43713d-3)*t
     x-0.1194975
      a = (d*sr + c)*s + aw
c
      b1 = (-5.3009d-4*t+1.6483d-2)*t+7.944d-2
      a1 = ((-6.1670d-5*t+1.09987d-2)*t-0.603459)*t+54.6746
      kw = (((-5.155288d-5*t+1.360477d-2)*t-2.327105)*t
     x+148.4206)*t-1930.06
      k0 = (b1*sr + a1)*s + kw
c evaluate pressure polynomial
c ***********************************************
c   k equals the secant bulk modulus of seawater
c   dk=k(s,t,p)-k(35,0,p)
c  k35=k(35,0,p)
c ***********************************************
      dk = (b*p + a)*p + k0
      k35  = (5.03217d-5*p+3.359406)*p+21582.27
      gam=p/k35
      pk = 1.0 - gam
      sva = sva*pk + (v350p+sva)*p*dk/(k35*(k35+dk))
c  scale specific vol. anamoly to normally reported units
      svan=sva*1.0d+8
      v350p = v350p*pk
c  ****************************************************
c compute density anamoly with respect to 1000.0 kg/m**3
c  1) dr350: density anamoly at 35 (ipss-78), 0 deg. c and 0 decibars
c  2) dr35p: density anamoly 35 (ipss-78), 0 deg. c ,  pres. variation
c  3) dvan : density anamoly variations involving specfic vol. anamoly
c ********************************************************************
c check value: sigma = 59.82037  kg/m**3 for s = 40 (ipss-78),
c t = 40 deg c, p0= 10000 decibars.
c *******************************************************
      dr35p=gam/v350p
      dvan=sva/(v350p*(v350p+sva))
      sigma=dr350+dr35p-dvan
      return
      end

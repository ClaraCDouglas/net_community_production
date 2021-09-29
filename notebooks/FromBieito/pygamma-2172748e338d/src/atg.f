      double precision function atg(s,t,p)
*
c ****************************
c adiabatic temperature gradient deg c per decibar
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c units:
c       pressure        p        decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       adiabatic       atg      deg. c/decibar
c checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
c t=40 deg c,p0=10000 decibars
*
      implicit double precision (a-h,o-z)
*

      ds = s - 35.0

      atg = (((-2.1687d-16*t+1.8676d-14)*t-4.6206d-13)*p
     x+((2.7759d-12*t-1.1351d-10)*ds+((-5.4481d-14*t
     x+8.733d-12)*t-6.7795d-10)*t+1.8741d-8))*p
     x+(-4.2393d-8*t+1.8932d-6)*ds
     x+((6.6228d-10*t-6.836d-8)*t+8.5258d-6)*t+3.5803d-5


      return
      end

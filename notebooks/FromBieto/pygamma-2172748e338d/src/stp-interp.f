        subroutine stp_interp(s,t,p,n,s0,t0,p0)
ccc
ccc
ccc
ccc     DESCRIPTION :   Interpolate salinity and in situ temperature
ccc                     on a cast by linearly interpolating salinity
ccc                     and potential temperature
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s(n)            array of cast salinities
ccc                     t(n)            array of cast in situ temperatures
ccc                     p(n)            array of cast pressures
ccc                     n               length of cast
ccc                     p0              pressure for which salinity and
ccc                                     in situ temperature are required
ccc
ccc     OUTPUT :        s0              interpolated value of salinity
ccc                     t0              interpolated value of situ temperature
ccc
ccc     UNITS :         salinities      psu (IPSS-78)
ccc                     temperatures    degrees C (IPTS-68)
ccc                     pressures       db
ccc
ccc
ccc     AUTHOR :        David Jackett
ccc
ccc     CREATED :       June 1993
ccc
ccc     REVISION :      1.1             30/6/93
ccc
ccc
ccc
        implicit double precision (a-h,o-z)

        dimension s(n),t(n),p(n)

        external indx

        data pr0/0.0/



        call indx(p,n,p0,k)

        r = (p0-p(k))/(p(k+1)-p(k))

        s0 = s(k) + r*(s(k+1)-s(k))

        thk = theta(s(k),t(k),p(k),pr0)

        th0 = thk + r*(theta(s(k+1),t(k+1),p(k+1),pr0)-thk)

        t0 = theta(s0,th0,pr0,p0)



        return
        end


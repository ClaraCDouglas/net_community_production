        subroutine scv_solve(s,t,p,e,n,k,s0,t0,p0,sscv,tscv,pscv,iter)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the zero of the v function using a
ccc                     bisection method
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s(n)            array of cast salinities
ccc                     t(n)            array of cast in situ temperatures
ccc                     p(n)            array of cast pressures
ccc                     e(n)            array of cast e values
ccc                     n               length of cast
ccc                     k               interval (k-1,k) contains the zero
ccc                     s0              the bottle salinity
ccc                     t0              the bottle in situ temperature
ccc                     p0              the bottle pressure
ccc
ccc     OUTPUT :        sscv            salinity of the scv surface
ccc                                     intersection with the cast
ccc                     tscv            in situ temperature of the intersection
ccc                     pscv            pressure of the intersection
ccc
ccc
ccc     UNITS :         salinities      psu (IPSS-78)
ccc                     temperatures    degrees C (IPTS-68)
ccc                     pressures       db
ccc
ccc
ccc     AUTHOR :        David Jackett
ccc
ccc     CREATED :       February 1995
ccc
ccc     REVISION :      1.1             1/2/95
ccc
ccc
ccc
        implicit double precision (a-h,o-z)

        dimension s(n),t(n),p(n),e(n)

        data n2/2/



        pl = p(k-1)
        el = e(k-1)
        pu = p(k)
        eu = e(k)

        iter = 0
        isuccess = 0

        do while(isuccess.eq.0)

          iter = iter+1

          pm = (pl+pu)/2

          call stp_interp(s(k-1),t(k-1),p(k-1),n2,sm,tm,pm)

          sdum = svan(s0,theta(s0,t0,p0,pm),pm,sigl)
          sdum = svan(sm,tm,pm,sigu)
          em = sigu-sigl

          if(el*em.lt.0.) then
            pu = pm
            eu = em
          elseif(em*eu.lt.0.) then
            pl = pm
            el = em
          elseif(em.eq.0.) then
            sscv = sm
            tscv = tm
            pscv = pm
            isuccess = 1
          end if

          if(isuccess.eq.0) then
            if(abs(em).le.5.d-5.and.abs(pu-pl).le.5.d-3) then
              sscv = sm
              tscv = tm
              pscv = pm
              isuccess = 1
            elseif(iter.le.20) then
              isuccess = 0
            else
              print *, 'WARNING 1 in scv-solve.f'
              print *, iter,'  em',abs(em),'  dp',pl,pu,abs(pu-pl)
              sscv = -99.0
              tscv = -99.0
              pscv = -99.0
              isuccess = 1
            end if
          end if

        end do



        return

        end

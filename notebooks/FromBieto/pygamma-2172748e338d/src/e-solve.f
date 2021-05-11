        subroutine e_solve(s,t,p,e,n,k,s0,t0,p0,sns,tns,pns,iter)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the zero of the e function using a
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
ccc     OUTPUT :        sns             salinity of the neutral surface
ccc                                     intersection with the cast
ccc                     tns             in situ temperature of the intersection
ccc                     pns             pressure of the intersection
ccc
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

          call sig_vals(s0,t0,p0,sm,tm,pm,sigl,sigu)
          em = sigu-sigl

          if(el*em.lt.0.) then
            pu = pm
            eu = em
          elseif(em*eu.lt.0.) then
            pl = pm
            el = em
          elseif(em.eq.0.) then
            sns = sm
            tns = tm
            pns = pm
            isuccess = 1
          end if

          if(isuccess.eq.0) then
            if(abs(em).le.5.d-5.and.abs(pu-pl).le.5.d-3) then
              sns = sm
              tns = tm
              pns = pm
              isuccess = 1
            elseif(iter.le.20) then
              isuccess = 0
            else
              print *, 'WARNING 1 in e-solve.f'
              print *, iter,'  em',abs(em),'  dp',pl,pu,abs(pu-pl)
              sns = -99.0
              tns = -99.0
              pns = -99.0
              isuccess = 1
            end if
          end if

        end do



        return

        end

        subroutine goor_solve(sl,tl,el,su,tu,eu,p,s0,t0,p0,sigb,sns,tns)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the intersection of a potential density surface
ccc                     between two bottles using a bisection method
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         sl, su          bottle salinities
ccc                     tl, tu          bottle in situ temperatures
ccc                     el, eu          bottle e values
ccc                     p               bottle pressures (the same)
ccc                     s0              emanating bottle salinity
ccc                     t0              emanating bottle in situ temperature
ccc                     p0              emanating bottle pressure
ccc
ccc     OUTPUT :        sns             salinity of the neutral surface
ccc                                     intersection with the bottle pair
ccc                     tns             in situ temperature of the intersection
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



        rl = 0.
        ru = 1.

        pmid = (p+p0)/2.0

        thl = theta(sl,tl,p,pmid)
        thu = theta(su,tu,p,pmid)

        iter = 0
        isuccess = 0

        do while(isuccess.eq.0)

          iter = iter+1

          rm = (rl+ru)/2

          sm = sl+rm*(su-sl)
          thm = thl+rm*(thu-thl)

          tm = theta(sm,thm,pmid,p)

          sd = svan(sm,thm,pmid,sigma)

          em = sigma-sigb

          if(el*em.lt.0.) then
            ru = rm
            eu = em
          elseif(em*eu.lt.0.) then
            rl = rm
            el = em
          elseif(em.eq.0.) then
            sns = sm
            tns = tm
            isuccess = 1
          end if

          if(isuccess.eq.0) then
            if(abs(em).le.5.e-5.and.abs(ru-rl).le.5.e-3) then
              sns = sm
              tns = tm
              isuccess = 1
            elseif(iter.le.20) then
              isuccess = 0
            else
              print *, 'WARNING 1 in goor-solve.f'
              sns = sm
              tns = tm
              isuccess = 1
            end if
          end if

        end do



        return

        end

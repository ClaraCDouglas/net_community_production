        subroutine depth_scv(s,t,p,n,s0,t0,p0,sscv,tscv,pscv,nscv)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the position which the scv surface through a
ccc                     specified bottle intersects a neighbouring cast
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s(n)            array of cast salinities
ccc                     t(n)            array of cast in situ temperatures
ccc                     p(n)            array of cast pressures
ccc                     n               length of cast
ccc                     s0              the bottle salinity
ccc                     t0              the bottle in situ temperature
ccc                     p0              the bottle pressure
ccc
ccc     OUTPUT :        sscv            salinities of the scv surface
ccc                                     intersections with the cast
ccc                     tscv            temperatures of the intersections
ccc                     pscv            pressures of the intersections
ccc                     nscv            number of intersections
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
ccc     REVISION :      1.1             9/2/95
ccc
ccc
ccc
        implicit double precision (a-h,o-z)

        parameter(n_max=2000,nscv_max=50)

        dimension s(n),t(n),p(n),e(n_max)
        dimension sscv(nscv_max),tscv(nscv_max),pscv(nscv_max)

        data n2/2/


        ierr = 0
        if(n.gt.n_max) then
          print *, '\nparameter n_max in depth-scv.f < ',n,'\n'
ccc       stop
        ierr = 9
        return
        end if

ccc
ccc             find the bottle pairs containing a crossing
ccc

        ncr = 0
        nscv = 0

        do k = 1,n

          sdum = svan(s0,theta(s0,t0,p0,p(k)),p(k),sigl)
          sdum = svan(s(k),t(k),p(k),sigu)
          e(k) = sigu-sigl

          if(k.gt.1) then

ccc
ccc             an exact crossing at the k-1 bottle
ccc

            if(e(k-1).eq.0.) then

              ncr = ncr+1
              sscv_tmp = s(k-1)
              tscv_tmp = t(k-1)
              pscv_tmp = p(k-1)


ccc
ccc             a crossing between k-1 and k bottles
ccc

            elseif(e(k)*e(k-1).lt.0.0) then

              ncr = ncr+1

ccc
ccc             some Newton-Raphson iterations to find the crossing
ccc

              pc0 = p(k-1)-e(k-1)*(p(k)-p(k-1))/(e(k)-e(k-1))

              iter = 0
              isuccess = 0

              do while(isuccess.eq.0)

                iter = iter+1

                call stp_interp(s(k-1),t(k-1),p(k-1),n2,sc0,tc0,pc0)

                sdum = svan(s0,theta(s0,t0,p0,pc0),pc0,sigl)
                sdum = svan(sc0,tc0,pc0,sigu)
                ec0 = sigu-sigl

                p1 = (p(k-1)+pc0)/2
                ez1 = (e(k-1)-ec0)/(pc0-p(k-1))
                p2 = (pc0+p(k))/2
                ez2 = (ec0-e(k))/(p(k)-pc0)
                r = (pc0-p1)/(p2-p1)
                ecz_0 = ez1+r*(ez2-ez1)

                if(iter.eq.1) then
                  ecz0 = ecz_0
                else
                  ecz0 = -(ec0-ec_0)/(pc0-pc_0)
                  if(ecz0.eq.0) ecz0 = ecz_0
                end if

                pc1 = pc0+ec0/ecz0

ccc
ccc             strategy when the iteration jumps out of the inteval
ccc

                if(pc1.le.p(k-1).or.pc1.ge.p(k)) then
                  call scv_solve(s,t,p,e,n,k,s0,t0,p0,
     &                                  sscv_tmp,tscv_tmp,pscv_tmp,niter)
                  if(pscv_tmp.lt.p(k-1).or.pscv_tmp.gt.p(k)) then
                    ierr = 10
                    return
ccc                 stop 'ERROR 1 in depth-scv.f'
                  else
                    isuccess = 1
                  endif
                else

ccc
ccc             otherwise, test the accuracy of the iterate
ccc

                  eps = abs(pc1-pc0)

                  if(abs(ec0).le.5.e-5.and.eps.le.5.e-3) then
                    sscv_tmp = sc0
                    tscv_tmp = tc0
                    pscv_tmp = pc0
                    isuccess = 1
                    niter = iter
                  elseif(iter.gt.10) then
                    call scv_solve(s,t,p,e,n,k,s0,t0,p0,
     &                                  sscv_tmp,tscv_tmp,pscv_tmp,niter)
                    isuccess = 1
                  else
                    pc_0 = pc0
                    ec_0 = ec0
                    pc0 = pc1
                    isuccess = 0
                  end if

                end if

              end do

            end if

          end if

ccc
ccc             the last bottle
ccc

          if(k.eq.n.and.e(k).eq.0.0) then
            ncr = ncr+1
            sscv_tmp = s(k)
            tscv_tmp = t(k)
            pscv_tmp = p(k)
          end if


ccc
ccc             store multiples
ccc

          if(ncr.gt.nscv) then
            nscv = nscv+1
            if(nscv.gt.nscv_max) then
            ierr = 10
            return
ccc     stop 'ERROR 2 in depth-scv.f'
        end if


            sscv(nscv) = sscv_tmp
            tscv(nscv) = tscv_tmp
            pscv(nscv) = pscv_tmp
          end if

        end do


ccc
ccc             no crossings
ccc

        if(nscv.eq.0) then
          sscv(1) = -99.0
          tscv(1) = -99.0
          pscv(1) = -99.0
        end if



        return
        end

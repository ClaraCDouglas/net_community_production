        subroutine depth_ns(s,t,p,n,s0,t0,p0,sns,tns,pns, ierr)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the position which the neutral surface through a
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
ccc     OUTPUT :        sns             salinity of the neutral surface
ccc                                     intersection with the cast
ccc                     tns             in situ temperature of the intersection
ccc                     pns             pressure of the intersection
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

        parameter(nmax=100)

        dimension s(n),t(n),p(n),e(nmax)

        data n2/2/


        ierr = 0
        if(n.gt.nmax) then
        ierr = 7
          print *, '\nparameter nmax in depth-ns.f < ',n,'\n'
        return
ccc       stop
        end if

ccc
ccc             find the bottle pairs containing a crossing
ccc

        ncr = 0

        do k = 1,n

          call sig_vals(s0,t0,p0,s(k),t(k),p(k),sigl,sigu)
          e(k) = sigu-sigl

          if(k.gt.1) then

ccc
ccc             an exact crossing at the k-1 bottle
ccc

            if(e(k-1).eq.0.) then

              ncr = ncr+1
              sns = s(k-1)
              tns = t(k-1)
              pns = p(k-1)

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

                call sig_vals(s0,t0,p0,sc0,tc0,pc0,sigl,sigu)
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
                  call e_solve(s,t,p,e,n,k,s0,t0,p0,sns,tns,pns,niter)
                  if(pns.lt.p(k-1).or.pns.gt.p(k)) then
                ierr = 8
                return
ccc                 stop 'ERROR 1 in depth-ns.f'
                  else
                    isuccess = 1
                  endif
                else

ccc
ccc             otherwise, test the accuracy of the iterate
ccc

                  eps = abs(pc1-pc0)

                  if(abs(ec0).le.5.e-5.and.eps.le.5.e-3) then
                    sns = sc0
                    tns = tc0
                    pns = pc0
                    isuccess = 1
                    niter = iter
                  elseif(iter.gt.10) then
                    call e_solve(s,t,p,e,n,k,s0,t0,p0,sns,tns,pns,niter)
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
            sns = s(k)
            tns = t(k)
            pns = p(k)
          end if

        end do

ccc
ccc             multiple crossings
ccc

        if(ncr.eq.0) then
          sns = -99.0
          tns = -99.0
          pns = -99.0
        elseif(ncr.ge.2) then
c         print *, 'WARNING in depth-ns.f: multiple crossings'
          sns = -99.2
          tns = -99.2
          pns = -99.2
        end if



        return
        end

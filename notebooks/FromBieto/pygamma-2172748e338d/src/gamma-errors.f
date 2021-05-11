        subroutine gamma_errors(s,t,p,gamma,a,n,along,alat,
     &                                  s0,t0,p0,sns,tns,pns,kns,
     &                              gamma_ns,pth_error,scv_l_error,scv_h_error, ierr)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the p-theta and the scv errors associated
ccc                     with the basic neutral surface calculation
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s(n)            array of Levitus cast salinities
ccc                     t(n)            array of cast in situ temperatures
ccc                     p(n)            array of cast pressures
ccc                     gamma(n)        array of cast neutral densities
ccc                     a(n)            array of cast quadratic coefficients
ccc                     n               length of cast
ccc                     along           longitude of Levitus cast
ccc                     alat            latitude of Levitus cast
ccc                     s0              bottle salinity
ccc                     t0              bottle temperature
ccc                     p0              bottle pressure
ccc                     sns             salinity of neutral surface on cast
ccc                     tns             temperature of neutral surface on cast
ccc                     pns             pressure of neutral surface on cast
ccc                     kns             index of neutral surface on cast
ccc                     gamma_ns        gamma value of neutral surface on cast
ccc
ccc     OUTPUT :        pth_error       p-theta gamma error bar
ccc                     scv_l_error     lower scv gamma error bar
ccc                     scv_h_error     upper scv gamma error bar
ccc
ccc     UNITS :         salinity        psu (IPSS-78)
ccc                     temperature     degrees C (IPTS-68)
ccc                     pressure        db
ccc                     gamma           kg m-3
ccc
ccc
ccc     AUTHOR :        David Jackett
ccc
ccc     CREATED :       March 1995
ccc
ccc     REVISION :      1.1             9/3/95
ccc
ccc
ccc
        parameter(nscv_max=50)

        implicit double precision (a-h,o-z)

        dimension s(n),t(n),p(n),gamma(n),a(n)

        dimension sscv_m(nscv_max),tscv_m(nscv_max),pscv_m(nscv_max)


        data pr0/0.0/, Tb/2.7e-8/, gamma_limit/26.845/, test_limit/0.1/


        ierr = 0
ccc
ccc             p-theta error
ccc

        th0 = theta(s0,t0,p0,pr0)
        thns = theta(sns,tns,pns,pr0)

        sdum = svan(sns,tns,pns,sig_ns)
        rho_ns = 1000+sig_ns

        call sig_vals(s(kns),t(kns),p(kns),s(kns+1),t(kns+1),p(kns+1),
     &                                                          sig_l,sig_h)

        b = (gamma(kns+1)-gamma(kns))/(sig_h-sig_l)

        dp = pns-p0
        dth = thns-th0

        pth_error = rho_ns*b*Tb*abs(dp*dth)/6


ccc
ccc             scv error
ccc

        scv_l_error = 0.0
        scv_h_error = 0.0

        if(alat.le.-60.0.or.gamma(1).ge.gamma_limit) then

          drldp = (sig_h-sig_l)/(rho_ns*(p(kns+1)-p(kns)))

          test = Tb*dth/drldp


cc
cc              approximation
cc

          if(abs(test).le.test_limit) then

            if(dp*dth.ge.0.0) then
              scv_h_error = (3*pth_error)/(1.0-test)
            else
              scv_l_error = (3*pth_error)/(1.0-test)
            end if

          else


cc
cc              explicit scv solution, when necessary
cc

            call depth_scv(s,t,p,n,s0,t0,p0,sscv_m,tscv_m,pscv_m,nscv, ierr)
            if(ierr.ne.0) then
                return
            end if

            if(nscv.eq.0) then

              continue

            else

              if(nscv.eq.1) then

                pscv = pscv_m(1)

              else

                pscv_mid = pscv_m((1+nscv)/2)

                if(p0.le.pscv_mid) then
                  pscv = pscv_m(1)
                else
                  pscv = pscv_m(nscv)
                end if

              end if

              call indx(p,n,pscv,kscv)
              call gamma_qdr(p(kscv),gamma(kscv),a(kscv),
     &                          p(kscv+1),gamma(kscv+1),pscv,gamma_scv)

              if(pscv.le.pns) then
                scv_l_error = gamma_ns-gamma_scv
              else
                scv_h_error = gamma_scv-gamma_ns
              end if

            end if

          end if

        else

          continue

        end if


ccc
ccc             check for positive gamma errors
ccc

        if(pth_error.lt.0.0.or.
     &          scv_l_error.lt.0.0.or.
     &                  scv_h_error.lt.0.0) then
        ierr = 6
        return
ccc       stop 'ERROR 1 in gamma-errors: negative scv error'

        end if



        return
        end


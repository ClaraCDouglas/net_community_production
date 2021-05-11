        subroutine goor(s,t,p,gamma,n,sb,tb,pb,
     &                                  gammab,g1_err,g2_l_err,g2_h_err, ierr)
ccc
ccc
ccc
ccc     DESCRIPTION :   Extend a cast of hydrographic data so that
ccc                     a bottle outside the gamma range of the cast
ccc                     can be labelled with the neutral density variable
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s(n)            array of cast salinities
ccc                     t(n)            array of cast in situ temperatures
ccc                     p(n)            array of cast pressures
ccc                     gamma(n)        array of cast gammas
ccc                     n               length of cast
ccc                     sb              bottle salinity
ccc                     tb              bottle temperature
ccc                     pb              bottle pressure
ccc
ccc     OUTPUT :        gammab          bottle gamma value
ccc                     g1_err          bottle Type i error estimate
ccc                     g2_l_err        bottle Type ii lower error estimate
ccc                     g2_h_err        bottle Type ii upper error estimate
ccc
ccc     UNITS :         salinity        psu (IPSS-78)
ccc                     temperature     degrees C (IPTS-68)
ccc                     pressure        db
ccc                     gamma           kg m-3
ccc
ccc
ccc     AUTHOR :        David Jackett
ccc
ccc     CREATED :       June 1993
ccc
ccc     REVISION :      1.2             2/11/94
ccc
ccc
ccc
        implicit double precision (a-h,o-z)

        dimension s(n),t(n),p(n),gamma(n)

        data delt_b/-0.1/, delt_t/0.1/, slope/-0.14/, pr0/0.0/, Tbp/2.7e-8/

        ierr = 0

ccc
ccc             determine if its bottom data
ccc

        pmid = (p(n)+pb)/2.0

        sd = svan(s(n),theta(s(n),t(n),p(n),pmid),pmid,sigma)

        sd = svan(sb,theta(sb,tb,pb,pmid),pmid,sigb)

ccc
ccc             a bottom extension
ccc

        if(sigb.gt.sigma) then

cc
cc              extend the cast data till it is denser
cc

          n_sth = 0
          s_new = s(n)
          t_new = t(n)
          e_new = sigma-sigb

          do while (sigma.lt.sigb)
            s_old = s_new
            t_old = t_new
            e_old = e_new
            n_sth = n_sth+1
            s_new = s(n)+n_sth*delt_b*slope
            t_new = t(n)+n_sth*delt_b
            sd = svan(s_new,theta(s_new,t_new,p(n),pmid),pmid,sigma)
            e_new = sigma-sigb
          end do

cc
cc              find the salinity and temperature with
cc              the same neutral density
cc

          if(sigma.eq.sigb) then
            sns = s_new
            tns = t_new
          else
            call goor_solve(s_old,t_old,e_old,s_new,t_new,e_new,p(n),
     &                                              sb,tb,pb,sigb,sns,tns)
          end if

cc
cc              now compute the new gamma value
cc

          call sig_vals(s(n-1),t(n-1),p(n-1),s(n),t(n),p(n),sigl,sigu)
          bmid = (gamma(n)-gamma(n-1))/(sigu-sigl)

          sd = svan(s(n),t(n),p(n),sigl)
          sd = svan(sns,tns,p(n),sigu)

          gammab = gamma(n)+bmid*(sigu-sigl)

          pns = p(n)

        else

ccc
ccc             determine if its top data
ccc

          pmid = (p(1)+pb)/2.0

          sd = svan(s(1),theta(s(1),t(1),p(1),pmid),pmid,sigma)

          sd = svan(sb,theta(sb,tb,pb,pmid),pmid,sigb)

ccc
ccc             a top extension
ccc

          if(sigb.lt.sigma) then

cc
cc              extend the cast data till it is lighter
cc

            n_sth = 0
            s_new = s(1)
            t_new = t(1)
            e_new = sigma-sigb
            do while (sigma.gt.sigb)
              s_old = s_new
              t_old = t_new
              e_old = e_new
              n_sth = n_sth+1
              s_new = s(1)
              t_new = t(1)+n_sth*delt_t
              sd = svan(s_new,theta(s_new,t_new,p(1),pmid),pmid,sigma)
              e_new = sigma-sigb
            end do

cc
cc              find the salinity and temperature with
cc              the same neutral density
cc

            if(sigma.eq.sigb) then
              sns = s_new
              tns = t_new
            else
              call goor_solve(s_new,t_new,e_new,s_old,t_old,e_old,p(1),
     &                                              sb,tb,pb,sigb,sns,tns)
            end if

cc
cc              now compute the new gamma value
cc

            call sig_vals(s(1),t(1),p(1),s(2),t(2),p(2),sigl,sigu)
            bmid = (gamma(2)-gamma(1))/(sigu-sigl)

            sd = svan(sns,tns,p(1),sigl)
            sd = svan(s(1),t(1),p(1),sigu)

            gammab = gamma(1)-bmid*(sigu-sigl)

            pns = p(1)

ccc
ccc             neither top nor bottom extension
ccc

          else
        ierr = 5
        return
ccc         stop 'ERROR 1 in gamma-out-of-range.f'

          end if

        end if


ccc
ccc             error estimate
ccc


        thb = theta(sb,tb,pb,pr0)
        thns = theta(sns,tns,pns,pr0)

        sdum = svan(sns,tns,pns,sig_ns)
        rho_ns = 1000+sig_ns

        b = bmid

        dp = pns-pb
        dth = thns-thb

        g1_err = rho_ns*b*Tbp*abs(dp*dth)/6

        g2_err = rho_ns*b*Tbp*dp*dth/2

        if(g2_err.le.0.0) then
          g2_l_err = -g2_err
          g2_h_err = 0.0
        else
          g2_l_err = 0.0
          g2_h_err = g2_err
        end if



        return
        end

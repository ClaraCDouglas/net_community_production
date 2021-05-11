        subroutine neutral_surfaces(s,t,p,gamma,n,glevels,ng,
     &                                          sns,tns,pns,dsns,dtns,dpns, ierrout)
ccc
ccc
ccc
ccc     DESCRIPTION :   For a cast of hydrographic data which has been
ccc                     labelled with the neutral density variable gamma,
ccc                     find the salinities, temperatures and pressures
ccc                     on ng specified neutral density surfaces.
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s(n)            array of cast salinities
ccc                     t(n)            array of cast in situ temperatures
ccc                     p(n)            array of cast pressures
ccc                     gamma(n)        array of cast gamma values
ccc                     n               length of cast
ccc                     glevels(ng)     array of neutral density values
ccc                     ng              number of neutral density surfaces
ccc
ccc     OUTPUT :        sns(ng)         salinity on the neutral density surfaces
ccc                     tns(ng)         in situ temperature on the surfaces
ccc                     pns(ng)         pressure on the surfaces
ccc                     dsns(ng)        surface salinity errors
ccc                     dtns(ng)        surface temperature errors
ccc                     dpns(ng)        surface pressure errors
ccc
ccc                     NOTE:           sns, tns and pns values of -99.0
ccc                                     denotes under or outcropping
ccc
ccc                                     non-zero dsns, dtns and dpns values
ccc                                     indicates multiply defined surfaces,
ccc                                     and file 'ns-multiples.dat' contains
ccc                                     information on the multiple solutions
ccc
ccc     UNITS :         salinity        psu (IPSS-78)
ccc                     temperature     degrees C (IPTS-68)
ccc                     pressure        db
ccc                     gamma           kg m-3
ccc
ccc
ccc     AUTHOR :        David Jackett
ccc
ccc     CREATED :       July 1993
ccc
ccc     REVISION :      2.1             17/5/95
ccc
ccc
ccc
        parameter(nint_max=50)

        implicit double precision (a-h,o-z)

        integer int(nint_max)

        dimension s(n),t(n),p(n),gamma(n)
        dimension glevels(ng),sns(ng),tns(ng),pns(ng)
        dimension dsns(ng),dtns(ng),dpns(ng)

        double precision alfa_l,beta_l,alfa_u,beta_u,alfa_mid,beta_mid
        double precision rhomid,thl,thu,dels,delth,pl,pu,delp,delp2,bmid
        double precision a,b,c,q

        data n2/2/, pr0/0.0/, ptol/1.0e-3/




ccc
ccc             detect error condition
ccc

        in_error = 0
        do k = 1,n
          if(gamma(k).le.0.d0) in_error = 1
        end do

        ierrout = 0
        if(in_error.eq.1) then
            ierrout = 2
            return
        end if
ccc     &         stop '\nERROR 1 in neutral-surfaces.f : missing gamma value'


ccc
ccc             loop over the surfaces
ccc

c       call system('rm -f ns-multiples.dat')

        ierr = 0

        do ig = 1,ng


ccc
ccc             find the intervals of intersection
ccc

          nint = 0

          do k = 1,n-1

            gmin = min(gamma(k),gamma(k+1))
            gmax = max(gamma(k),gamma(k+1))

            if(gmin.le.glevels(ig).and.glevels(ig).le.gmax) then
              nint = nint+1
              if(nint.gt.nint_max) then
            ierrout = 3
            return
       end if
cccstop 'ERROR 2 in neutral-surfaces.f'
              int(nint) = k
            end if

          end do


ccc
ccc             find point(s) of intersection
ccc

          if(nint.eq.0) then

            sns(ig) = -99.0
            tns(ig) = -99.0
            pns(ig) = -99.0
            dsns(ig) = 0.0
            dtns(ig) = 0.0
            dpns(ig) = 0.0

          else


ccc
ccc             choose the central interval
ccc

            if(mod(nint,2).eq.0.and.int(1).gt.n/2) then
              int_middle = (nint+2)/2
            else
              int_middle = (nint+1)/2
            end if


ccc
ccc             loop over all intersections
ccc

            do i_int = 1,nint

              k = int(i_int)

ccc
ccc             coefficients of a quadratic for gamma
ccc

              pmid = (p(k)+p(k+1))/2.

              call alphabeta(s(k),t(k),p(k),alfa,beta)
              alfa_l = alfa
              beta_l = beta
              call alphabeta(s(k+1),t(k+1),p(k+1),alfa,beta)
              alfa_u = alfa
              beta_u = beta

              alfa_mid = (alfa_l+alfa_u)/2.0
              beta_mid = (beta_l+beta_u)/2.0

              call stp_interp(s(k),t(k),p(k),n2,smid,tmid,pmid)

              sd = svan(smid,tmid,pmid,sigmid)
              rhomid = 1000.+sigmid

              thl = theta(s(k),t(k),p(k),pr0)
              thu = theta(s(k+1),t(k+1),p(k+1),pr0)

              dels = s(k+1)-s(k)
              delth = thu-thl

              pl = p(k)
              pu = p(k+1)
              delp = pu-pl
              delp2 = delp*delp

              bden = rhomid*(beta_mid*dels-alfa_mid*delth)

              if(abs(bden).le.1.d-6) bden = 1.d-6

              bmid = (gamma(k+1)-gamma(k))/bden
     &

cc
cc              coefficients
cc

              a = dels*(beta_u-beta_l)-delth*(alfa_u-alfa_l)
              a = (a*bmid*rhomid)/(2*delp2)

              b = dels*(pu*beta_l-pl*beta_u) - delth*(pu*alfa_l-pl*alfa_u)
              b = (b*bmid*rhomid)/delp2

              c = dels*(beta_l*(pl-2.*pu)+beta_u*pl) -
     &                 delth*(alfa_l*(pl-2.*pu)+alfa_u*pl)
              c = gamma(k) + (bmid*rhomid*pl*c)/(2*delp2)
              c = c - glevels(ig)

ccc
ccc             solve the quadratic
ccc

              if(a.ne.0.d0.and.bden.ne.1.d-6) then

                q = -(b+sign(1.d0,b)*sqrt(b*b-4*a*c))/2.0

                pns1 = q/a
                pns2 = c/q

                if(pns1.ge.p(k)-ptol.and.pns1.le.p(k+1)+ptol) then
                  pns(ig) = min(p(k+1),max(pns1,p(k)))
                elseif(pns2.ge.p(k)-ptol.and.pns2.le.p(k+1)+ptol) then
                  pns(ig) = min(p(k+1),max(pns2,p(k)))
                else
        ierrout = 4
        return
ccc               stop 'ERROR 3 in neutral-surfaces.f'
                end if

              else

                rg = (glevels(ig)-gamma(k))/(gamma(k+1)-gamma(k))
                pns(ig) = p(k)+rg*(p(k+1)-p(k))

              end if

              call stp_interp(s,t,p,n,sns(ig),tns(ig),pns(ig))


ccc
ccc             write multiple values to file
ccc

              if(nint.gt.1) then

                if(ierr.eq.0) then
                  ierr = 1
                  call system('rm -f ns-multiples.dat')
                  call get_lunit(lun, ierrout)
                if (ierrout.ne.0) return
                  open(lun,file='ns-multiples.dat',status='unknown')
                end if

                if(i_int.eq.1) write(lun,*) ig,nint

                write(lun,*) sns(ig),tns(ig),pns(ig)


ccc
ccc             find median values and errors
ccc

                if(i_int.eq.1) then
                  sns_top = sns(ig)
                  tns_top = tns(ig)
                  pns_top = pns(ig)
                elseif(i_int.eq.int_middle) then
                  sns_middle = sns(ig)
                  tns_middle = tns(ig)
                  pns_middle = pns(ig)
                elseif(i_int.eq.nint) then
                  if((pns_middle-pns_top).gt.(pns(ig)-pns_middle)) then
                    dsns(ig) = sns_middle-sns_top
                    dtns(ig) = tns_middle-tns_top
                    dpns(ig) = pns_middle-pns_top
                  else
                    dsns(ig) = sns(ig)-sns_middle
                    dtns(ig) = tns(ig)-tns_middle
                    dpns(ig) = pns(ig)-pns_middle
                  end if
                  sns(ig) = sns_middle
                  tns(ig) = tns_middle
                  pns(ig) = pns_middle
                end if

              else

                dsns(ig) = 0.0
                dtns(ig) = 0.0
                dpns(ig) = 0.0

              end if

            end do

          end if


        end do


        if(ierr.eq.1) close(lun)




        return
        end

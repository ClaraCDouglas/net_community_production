        subroutine gamma_n(filename,s,t,p,n,along,alat,gamma,dg_lo,dg_hi, ierr)
ccc
ccc
ccc
ccc     DESCRIPTION :   Label a cast of hydrographic data at a specified
ccc                     location with neutral density
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s(n)            array of cast salinities
ccc                     t(n)            array of cast in situ temperatures
ccc                     p(n)            array of cast pressures
ccc                     n               length of cast (n=1 for single bottle)
ccc                     along           longitude of cast (0-360)
ccc                     alat            latitude of cast (-80,64)
ccc
ccc     OUTPUT :        gamma(n)        array of cast gamma values
ccc                     dg_lo(n)        array of gamma lower error estimates
ccc                     dg_hi(n)        array of gamma upper error estimates
ccc
ccc                     NOTE:           -99.0 denotes algorithm failed
ccc                                     -99.1 denotes input data is outside
ccc                                           the valid range of the present
ccc                                           equation of state
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
ccc     REVISION :      3.1             27/11/95
ccc
ccc
ccc
        parameter(nx=90,ny=45,nz=33,ndx=4,ndy=4)

        implicit double precision (a-h,o-z)

        integer*4 ioce,iocean0(2,2),n0(2,2)

      character *(*) filename

        dimension s(n),t(n),p(n),gamma(n),dg_lo(n),dg_hi(n)
        dimension along0(2),alat0(2)
        dimension s0(nz,2,2),t0(nz,2,2),p0(nz),gamma0(nz,2,2),a0(nz,2,2)

        dimension gwij(4),wtij(4)

        data pr0/0.0/, dgamma_0/0.0005/, dgw_max/0.3/

        save

        external indx

        ierr = 0

ccc
ccc             detect error conditions
ccc

        if(along.lt.0.0) then
          along = along+360.0
          ialtered = 1
        elseif(along.ge.360.0) then
          along = along-360.0
          ialtered = 2
        else
          ialtered = 0
        end if

        if(along.lt.0.0.or.along.gt.360.0.or.
     &                          alat.lt.-80.0.or.alat.gt.64.0) then
          print *, '\nERROR 1 in gamma-n.f : out of oceanographic range'
          print *, '\n',along,alat,'  is not contained in [0,360]x[-80,64]'
          print *
          ierr = 1
           return
        end if

        do k = 1,n
          if(s(k).lt.0.0.or.s(k).gt.42.0.or.
     &           t(k).lt.-2.5.or.t(k).gt.40.0.or.
     &               p(k).lt.0.0.or.p(k).gt.10000.0) then
            gamma(k) = -99.1
            dg_lo(k) = -99.1
            dg_hi(k) = -99.1
          else
            gamma(k) = 0.0
            dg_lo(k) = 0.0
            dg_hi(k) = 0.0
          end if
        end do


ccc
ccc             read records from the netCDF data file
ccc
        call read_nc(filename,along,alat,s0,t0,p0,gamma0,a0,n0,along0,alat0,iocean0, ierr)
        if (ierr.ne.0) return

ccc
ccc             find the closest cast
ccc

        dist2_min = 1.e10

        do j0 = 1,2
        do i0 = 1,2

          if(n0(i0,j0).ne.0) then
            dist2 = (along0(i0)-along)*(along0(i0)-along) +
     &                                  (alat0(j0)-alat)*(alat0(j0)-alat)
            if(dist2.lt.dist2_min) then
              i_min = i0
              j_min = j0
              dist2_min = dist2
            end if
          end if

        end do
        end do

        ioce = iocean0(i_min,j_min)


ccc
ccc             label the cast
ccc

        dx = abs(mod(along,dble(ndx)))
        dy = abs(mod(alat+80.0,dble(ndy)))
        rx = dx/dble(ndx)
        ry = dy/dble(ndy)

        do k = 1,n
        if(gamma(k).ne.-99.1) then

          thk = theta(s(k),t(k),p(k),pr0)

          dgamma_1 = 0.0
          dgamma_2_l = 0.0
          dgamma_2_h = 0.0

          wsum = 0.0

          nij = 0


ccc
ccc             average the gammas over the box
ccc

          do j0 = 1,2
          do i0 = 1,2
          if(n0(i0,j0).ne.0) then

            if(j0.eq.1) then
              if(i0.eq.1) then
                wt = (1.-rx)*(1.-ry)
              elseif(i0.eq.2) then
                wt = rx*(1-ry)
              end if
            elseif(j0.eq.2) then
              if(i0.eq.1) then
                wt = (1.-rx)*ry
              elseif(i0.eq.2) then
                wt = rx*ry
              end if
            end if

            wt = wt+1.e-6

            call ocean_test(along,alat,ioce,along0(i0),alat0(j0),
     &                                          iocean0(i0,j0),p(k),itest)

            if(itest.eq.0) wt = 0.0

            call depth_ns(s0(1,i0,j0),t0(1,i0,j0),p0,n0(i0,j0),
     &                                  s(k),t(k),p(k),sns,tns,pns, ierr)
            if(ierr.ne.0) then
                return
            end if

            if(pns.gt.-99.) then

              call indx(p0,n0(i0,j0),pns,kns)
              call gamma_qdr(p0(kns),gamma0(kns,i0,j0),a0(kns,i0,j0),
     &                              p0(kns+1),gamma0(kns+1,i0,j0),pns,gw)


ccc
ccc             error bars
ccc
              call gamma_errors(s0(1,i0,j0),t0(1,i0,j0),p0,gamma0(1,i0,j0),
     &                           a0(1,i0,j0),n0(i0,j0),along0(i0),alat0(j0),
     &                            s(k),t(k),p(k),sns,tns,pns,kns,
     &                                          gw,g1_err,g2_l_err,g2_h_err, ierr)
            if (ierr.ne.0) return

            elseif(pns.eq.-99.) then

              call goor(s0(1,i0,j0),t0(1,i0,j0),p0,
     &                  gamma0(1,i0,j0),n0(i0,j0),s(k),t(k),p(k),
     &                                  gw,g1_err,g2_l_err,g2_h_err, ierr)
            if (ierr.ne.0) return

ccc
ccc             adjust weight for gamma extrapolation
ccc

              if(gw.gt.gamma0(n0(i0,j0),i0,j0)) then
                rw = min(dgw_max,gw-gamma0(n0(i0,j0),i0,j0))/dgw_max
                wt = (1-rw)*wt
              end if


            else
              gw = 0.0
              g1_err = 0.0
              g2_l_err = 0.0
              g2_h_err = 0.0
            end if

            if(gw.gt.0.) then
              gamma(k) = gamma(k)+wt*gw
              dgamma_1 = dgamma_1+wt*g1_err
              dgamma_2_l = max(dgamma_2_l,g2_l_err)
              dgamma_2_h = max(dgamma_2_h,g2_h_err)
              wsum = wsum+wt
              nij = nij+1
              wtij(nij) = wt
              gwij(nij) = gw
            end if

          end if
          end do
          end do


ccc
ccc             the average
ccc

          if(wsum.ne.0.0) then

            gamma(k) = gamma(k)/wsum
            dgamma_1 = dgamma_1/wsum


ccc
ccc             the gamma errors
ccc

            dgamma_3 = 0.0
            do ij = 1,nij
              dgamma_3 = dgamma_3+wtij(ij)*abs(gwij(ij)-gamma(k))
            end do
            dgamma_3 = dgamma_3/wsum

            dg_lo(k) = max(dgamma_0,dgamma_1,dgamma_2_l,dgamma_3)
            dg_hi(k) = max(dgamma_0,dgamma_1,dgamma_2_h,dgamma_3)

          else

            gamma(k) = -99.0
            dg_lo(k) = -99.0
            dg_hi(k) = -99.0

          end if

        end if
        end do


        if(ialtered.eq.1) then
          along = along-360.0
        elseif(ialtered.eq.2) then
          along = along+360.0
        end if




        return
        end


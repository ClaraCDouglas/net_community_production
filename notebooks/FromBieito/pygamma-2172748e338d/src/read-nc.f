        subroutine read_nc(filename,along,alat,s0,t0,p0,gamma0,a0,n0,
     &                                          along0,alat0,iocean0, ierr)
ccc
ccc
ccc
ccc     DESCRIPTION :   Read variables from the netcdf labelled data file
ccc
ccc     PRECISION:      Double precision
ccc
ccc     INPUT :
ccc            filename    full path to nc file
ccc             along           longitude of record
ccc                     alat            latitude of record
ccc
ccc     OUTPUT :        s0(nz,2,2)      array of cast salinities
ccc                     t0(nz,2,2)      array of cast in situ temperatures
ccc                     p0(nz)          array of cast pressures
ccc                     gamma0(nz,2,2)  array of cast gamma values
ccc                     a0(nz,2,2)      array of cast a values
ccc                     n0(2,2)         length of casts
ccc                     along0(2)       array of cast longitudes
ccc                     alat0(2)        array of cast latitudes
ccc                     iocean0(2,2)    array of cast oceans
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
ccc     REVISION :      1.3             15/11/94
ccc
ccc
ccc
        parameter(nx=90,ny=45,nz=33,ndx=4,ndy=4)

        implicit double precision (a-h,o-z)

        integer*4 id_gnc,id_lon,id_lat,id_p,id_n
        integer*4 id_ocean,id_s,id_t,id_gamma,id_a
        integer*4 start(3),count(3),n0(2,2),n0_t(2)

        integer*4 iocean(nx,ny),iocean0(2,2)

        real*4 along_s(nx),alat_s(ny)
        real*4 s0_s(nz,2,2),t0_s(nz,2,2),p0_s(nz),gamma0_s(nz,2,2),a0_s(nz,2,2)
        real*4 s0_t(nz,2),t0_t(nz,2),gamma0_t(nz,2),a0_t(nz,2)

        dimension along_d(nx),alat_d(ny),along0(2),alat0(2)
        dimension s0(nz,2,2),t0(nz,2,2),p0(nz),gamma0(nz,2,2),a0(nz,2,2)

        character*(*) filename

        save along_d, alat_d, i0, j0, iocean, p0_s

        data start/1,1,1/, count/nz,nx,ny/, nw/0/, i0/1/, j0/1/

        ilong(alng) = int(alng/ndx + 1)
        jlat(alt) = int((88+alt)/ndy + 1)

      ierr = 0

ccc
ccc             only read when you have to
ccc

        dx = along-along_d(i0)
        dy = alat-alat_d(j0)

        if(dx.lt.0.0.or.dx.ge.4.0.or.dy.lt.0.0.or.dy.ge.4.0.or.
     &                                      (i0.eq.1.and.j0.eq.1)) then

ccc
ccc             open the netcdf file and get the variable id's
ccc

C         filename = PWD//'/gamma.nc'
          ierr = nf_open(filename, 0, id_gnc)
            if (ierr.ne.0) then
              print *, ierr
                ierr = 13
                return
            end if
          ! Error handling above looks wrong...
          if(i0.eq.1.and.j0.eq.1) then

            ! These look unnecessary, duplicating data statement.
            start(2) = 1
            count(2) = nx
            start(3) = 1
            count(3) = ny

            ierr = nf_inq_varid(id_gnc, 'lon', id_lon)
            ierr = nf_get_var_double(id_gnc, id_lon, along_d)
            ierr = nf_inq_varid(id_gnc, 'lat', id_lat)
            ierr = nf_get_var_double(id_gnc, id_lat, alat_d)
            ierr = nf_inq_varid(id_gnc, 'pressure', id_p)
            ierr = nf_get_var_double(id_gnc, id_p, p0)

            ierr = nf_inq_varid(id_gnc, 'iocean', id_ocean)
            ierr = nf_get_var_int(id_gnc, id_ocean, iocean)

          end if

          ierr = nf_inq_varid(id_gnc, 'n', id_n)
          ierr = nf_inq_varid(id_gnc, 's', id_s)
          ierr = nf_inq_varid(id_gnc, 't', id_t)
          ierr = nf_inq_varid(id_gnc, 'gamma', id_gamma)
          ierr = nf_inq_varid(id_gnc, 'a', id_a)

ccc
ccc             read the appropriate records
ccc

          i0 = ilong(along)
          j0 = jlat(alat)

          along0(1) = along_d(i0)
          alat0(1) = alat_d(j0)
          alat0(2) = alat0(1)+ndy

          if(i0.eq.nx+1) i0 = 1

          if(i0.lt.nx) then

            along0(2) = along0(1)+ndx

            start(2) = i0
            count(2) = 2
            start(3) = j0
            count(3) = 2

            ierr = nf_get_vara_int(id_gnc, id_n, start(2), count(2), n0)
            ierr = nf_get_vara_double(id_gnc, id_s, start, count, s0)
            ierr = nf_get_vara_double(id_gnc, id_t, start, count, t0)
            ierr = nf_get_vara_double(id_gnc, id_gamma, start, count, gamma0)
            ierr = nf_get_vara_double(id_gnc, id_a, start, count, a0)


          elseif(i0.eq.nx) then

            start(2) = i0
            count(2) = 1
            start(3) = j0
            count(3) = 2

            ierr = nf_get_vara_int(id_gnc, id_n, start(2), count(2), n0_t)
            ierr = nf_get_vara_real(id_gnc, id_s, start, count, s0_t)
            ierr = nf_get_vara_real(id_gnc, id_t, start, count, t0_t)
            ierr = nf_get_vara_real(id_gnc, id_gamma, start, count, gamma0_t)
            ierr = nf_get_vara_real(id_gnc, id_a, start, count, a0_t)

            do j = 1,2
              n0(1,j) = n0_t(j)
              do k = 1,nz
                s0_s(k,1,j) = s0_t(k,j)
                t0_s(k,1,j) = t0_t(k,j)
                gamma0_s(k,1,j) = gamma0_t(k,j)
                a0_s(k,1,j) = a0_t(k,j)
              end do
            end do

            along0(2) = 0.0

            start(2) = 1
            count(2) = 1
            start(3) = j0
            count(3) = 2

            ierr = nf_get_vara_int(id_gnc, id_n, start(2), count(2), n0_t)
            ierr = nf_get_vara_real(id_gnc, id_s, start, count, s0_t)
            ierr = nf_get_vara_real(id_gnc, id_t, start, count, t0_t)
            ierr = nf_get_vara_real(id_gnc, id_gamma, start, count, gamma0_t)
            ierr = nf_get_vara_real(id_gnc, id_a, start, count, a0_t)


            do j = 1,2
              n0(2,j) = n0_t(j)
              do k = 1,nz
                s0_s(k,2,j) = s0_t(k,j)
                t0_s(k,2,j) = t0_t(k,j)
                gamma0_s(k,2,j) = gamma0_t(k,j)
                a0_s(k,2,j) = a0_t(k,j)
              end do
            end do

          end if

          ! call ncclos(id_gnc,ierr)
          ierr = nf_close(id_gnc)


ccc
ccc             get the ocean information
ccc

          do j = 1,2
          do i = 1,2
            iocean0(i,j) = iocean(ilong(along0(i)),jlat(alat0(j)))
          end do
          end do


        end if

        return
        end

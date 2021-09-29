ccc     Note the insertion of loops around the main function calls
ccc     so that we can get accurate profiling.

        parameter(nz=80,nlevels=3)

        implicit double precision (a-h,o-z)


        dimension s(nz),t(nz),p(nz),gamma(nz),dgl(nz),dgh(nz)
        dimension glevels(nlevels)
        dimension sns(nlevels),tns(nlevels),pns(nlevels)
        dimension dsns(nlevels),dtns(nlevels),dpns(nlevels)

ccc     Note the modification from the original: use double precision.

        data glevels/26.8_8,27.9_8,28.1_8/

        character(50) :: filename

ccc
ccc             an example of labelling data
ccc
        filename='../data/gamma.nc'

        open(10,file='../data/example.dat',status='old')

        read(10,*) along,alat,n
        do k = 1,n
          read(10,*) s(k),t(k),p(k)
        end do
ccc
ccc             label
ccc
        ! do k = 1,10000
        call gamma_n(filename,s,t,p,n,along,alat,gamma,dgl,dgh,ierr)
        ! end do
        print *
        print *, 'location'
        print *
        write(6,'(2f12.4,i8)') along,alat

        print *
        print *, 'labels'
        print *
        do k = 1,n
          write(6,'(f8.2,3f12.6)') p(k),gamma(k),dgl(k),dgh(k)
        end do


ccc
ccc             fit some surfaces
ccc
        ! do k = 1,10000
        call neutral_surfaces(s,t,p,gamma,n,glevels,nlevels,
     &                        sns,tns,pns,dsns,dtns,dpns, ierr)
        ! end do
        print *
        print *, 'surfaces'
        print *

        do k = 1,nlevels
          print '(f8.2,2f12.6,f14.6,f6.1)',
     &                          glevels(k),sns(k),tns(k),pns(k),dpns(k)
        end do

        print *




        end

        subroutine indx(x,n,z,k)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the index of a real number in a
ccc                     monotonically increasing real array
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         x               array of increasing values
ccc                     n               length of array
ccc                     z               real number
ccc
ccc     OUTPUT :        k               index k - if x(k) <= z < x(k+1), or
ccc                                     n-1     - if z = x(n)
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

        dimension x(n)


        if(x(1).lt.z.and.z.lt.x(n)) then

          kl=1
          ku=n

          do while (ku-kl.gt.1)
            km=(ku+kl)/2
            if(z.gt.x(km))then
              kl=km
            else
              ku=km
            endif
          end do

          k=kl

          if(z.eq.x(k+1)) k = k+1

        else

          if(z.eq.x(1)) then
            k = 1
          elseif(z.eq.x(n)) then
            k = n-1
          else
            print *, 'ERROR 1 in indx.f : out of range'
            print *, z,n,x
          end if

        end if



        return
        end

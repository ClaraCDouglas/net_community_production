        subroutine ocean_test(x1,y1,io1,x2,y2,io2,z,itest)
ccc
ccc
ccc
ccc     DESCRIPTION :   Test whether two locations are connected by ocean
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         x1              longitude of first location
ccc                     y1              latitude of first location
ccc                     io1             ocean of first location
ccc                     x2              longitude of second location
ccc                     y2              latitude of second location
ccc                     io2             ocean of second location
ccc                     z               depth of connection
ccc
ccc     OUTPUT :        itest           success of connection
ccc
ccc
ccc     AUTHOR :        David Jackett
ccc
ccc     CREATED :       June 1994
ccc
ccc     REVISION :      1.1             7/7/94
ccc
ccc
ccc
        implicit double precision (a-h,o-z)

        integer*4 io1,io2

        dimension x_js(3),y_js(3)

        data x_js/129.87, 140.37, 142.83/
        data y_js/ 32.75,  37.38,  53.58/



        y = (y1+y2)/2



ccc
ccc             same ocean talks
ccc

        if(io1.eq.io2) then

          itest = 1
          return

        elseif(y.le.-20.) then

ccc
ccc             land of South America doesn't talk
ccc

          if(y.ge.-48..and.(io1*io2).eq.12) then
            itest = 0

ccc
ccc             everything else south of -20 talks
ccc

          else
            itest = 1
          end if

ccc
ccc             Pacific talks
ccc

        elseif((io1.eq.1.or.io1.eq.2).and.
     &                    (io2.eq.1.or.io2.eq.2)) then
          itest = 1

ccc
ccc             Indian talks
ccc

        elseif((io1.eq.3.or.io1.eq.4).and.
     &                    (io2.eq.3.or.io2.eq.4)) then
          itest = 1

ccc
ccc             Atlantic talks
ccc

        elseif((io1.eq.5.or.io1.eq.6).and.
     &                    (io2.eq.5.or.io2.eq.6)) then
          itest = 1

ccc
ccc             Indonesian throughflow
ccc

        elseif(io1*io2.eq.8.and.z.le.1200..and.
     &                  x1.ge.124..and.x1.le.132..and.
     &                          x2.ge.124..and.x2.le.132.) then
          itest = 1

ccc
ccc             anything else doesn't
ccc
        else
          itest = 0
        end if


ccc
ccc             exclude Japan Sea from talking
ccc

        if( (x_js(1).le.x1.and.x1.le.x_js(3).and.
     &                  y_js(1).le.y1.and.y1.le.y_js(3)) .or.

     &      (x_js(1).le.x2.and.x2.le.x_js(3).and.
     &                  y_js(1).le.y2.and.y2.le.y_js(3)) ) then

          em1 = (y_js(2)-y_js(1))/(x_js(2)-x_js(1))
          c1 = y_js(1)-em1*x_js(1)

          em2 = (y_js(3)-y_js(2))/(x_js(3)-x_js(2))
          c2 = y_js(2)-em2*x_js(2)

          if((y1-em1*x1-c1).ge.0.0.and.(y1-em2*x1-c2).ge.0.0) then
            isj1 = 1
          else
            isj1 = 0
          end if

          if((y2-em1*x2-c1).ge.0.0.and.(y2-em2*x2-c2).ge.0.0) then
            isj2 = 1
          else
            isj2 = 0
          end if

          if(isj1.eq.isj2) then
            itest = 1
          else
            itest = 0
          end if

        end if

ccc
ccc             exclude Antarctic tip
ccc

        if(io1*io2.eq.12.and.y.lt.-60.) itest = 0



        return
        end

        subroutine gamma_qdr(pl,gl,a,pu,gu,p,gamma)
ccc
ccc
ccc
ccc     DESCRIPTION :   Evaluate the quadratic gamma profile at a pressure
ccc                     between two bottles
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         pl, pu          bottle pressures
ccc                     gl, gu          bottle gamma values
ccc                     a               quadratic coefficient
ccc                     p               pressure for gamma value
ccc
ccc     OUTPUT :        gamma           gamma value at p
ccc
ccc     UNITS :         pressure        db
ccc                     gamma           kg m-3
ccc                     a               kg m-3
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



        p1 = (p-pu)/(pu-pl)
        p2 = (p-pl)/(pu-pl)

        gamma = (a*p1+(gu-gl))*p2+gl



        return
        end





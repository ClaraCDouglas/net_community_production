        subroutine sig_vals(s1,t1,p1,s2,t2,p2,sig1,sig2)
ccc
ccc
ccc
ccc     DESCRIPTION :   Computes the sigma values of two neighbouring
ccc                     bottles w.r.t. the mid pressure
ccc
ccc     PRECISION :     Double
ccc
ccc     INPUT :         s1,s2           bottle salinities
ccc                     t1,t2           bottle in situ temperatures
ccc                     p1,p2           bottle pressures
ccc
ccc     OUTPUT :        sig1,sig2       bottle potential density values
ccc
ccc     UNITS :         salinity        psu (IPSS-78)
ccc                     temperature     degrees C (IPTS-68)
ccc                     pressure        db
ccc                     density         kg m-3
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



        pmid = (p1+p2)/2.0

        sd = svan(s1,theta(s1,t1,p1,pmid),pmid,sig1)

        sd = svan(s2,theta(s2,t2,p2,pmid),pmid,sig2)



        return
        end


        subroutine get_lunit(lun, ierr)
ccc
ccc
ccc
ccc     DESCRIPTION :   Find the first FORTRAN logical unit (>=20)
ccc                     which is available for writing
ccc
ccc     PRECISION :     Double
ccc
ccc     OUTPUT :        lun             available logical unit
ccc
ccc
ccc     AUTHOR :        David Jackett
ccc
ccc     CREATED :       October 1994
ccc
ccc     REVISION :      1.2             2/12/94
ccc
ccc
ccc

        implicit double precision (a-h,o-z)

        integer lun

        logical lv

        data lun0/20/, lun1/70/


        ierr = 0

        ifound = 0
        lun = lun0

        do while (lun.le.lun1.and.ifound.eq.0)

          inquire(unit=lun,opened=lv)

          if(lv) then
            lun = lun+1
          else
            ifound = 1
          end if

        end do

        if(ifound.eq.0) then
            ierr = 12
            return
        end if

cccstop 'ERROR 1 in get-lun.f'



        return
        end


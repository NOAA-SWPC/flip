!dbg20120330: copied from flip:RSSUBC.FOR
! called from module_field_line_grid.f90
!..  Evaluate the ExB term for the continuity equation
!        CALL GET_DIVDRIFT_VELOCITY(RE,PCO,X(J), DIV_VEM(J))
!RE:earth radius?
!PCO = lshell
!X(J)??? calculated in subroutine FIELD
!UNDERCONSTRUCTION!!!

C:::::::::::::::::::: GET_DIVDRIFT_VELOCITY ::::::::::::::::::::::::::::::
      SUBROUTINE GET_DIVDRIFT_VELOCITY (RE,PCO,X, DIV_VEM)
C
C.... This subroutine calculates the ExB term in the continuty equation.
C.... This expression comes from Bailey, PSS 1983, page 392. Note that Re in
C.... his equation (5) is the equatorial radius to the flux tube
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL RE
 
      NUMERATOR = 6.0000 * ((DSIN(X))**2.0)*(1.0+(DCOS(X))**2.0)
      DENOMINATOR = RE*PCO*1.0E5*((1.0+3.0*(DCOS(X))**2.0))**2.0
      DIV_VEM = NUMERATOR/DENOMINATOR
      RETURN
      END

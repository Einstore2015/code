#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic velocity variances \label{sec:variances}
!
! !INTERFACE:
   subroutine waveinduced_turbulence(nlev)
!
! !DESCRIPTION:

!  Using \eq{bijVertical} and the solution shown in \eq{b13} and
!  the variances of the turbulent velocity fluctations can be

!
! !USES:
  use meanflow,       only: z,u_taus
  use observations,   only: beta,wave_turbulence,Hs_,WTp_
  use turbulence,     only: diss_w
  use airsea_driver,  only: u10,v10
  use stokes_drift,   only: us0,uwnd,vwnd

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!

!  number of vertical layers
  integer,  intent(in)                :: nlev
  
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:

    integer                   :: i
    REALTYPE                  :: WA,WTS,WC,WK,WUS,WL,WSTAR,WUSURF,WVSURF,WUSTAR

    WUSURF = uwnd%value/1025.
    WVSURF = vwnd%value/1025.
    WUSTAR  = SQRT(WUSURF*WUSURF+WVSURF*WVSURF)
    WA  = 0.5 * Hs_%value
    WTS = 0.91 * WTp_%value
    WC  = 1.56 * WTS
    WL  = WC * WTS
    WK  = 2*3.1415/WL
    WUS = WC*WA*WA*WK*WK
    do i=1,nlev-1
        diss_w(i) = 148 * beta * sqrt(Hs_%value / WL) * WUS * u_taus**2 * &
                    exp(2.*WK*z(i)) / WL
    end do
    
   return
   end subroutine waveinduced_turbulence

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

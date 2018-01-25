module resources

use parameters_site
use parameters_plant

implicit none

! Resource variables
real :: DTRINT    ! = MJ GR m-2 d-1 Interception of global radiation
real :: PARAV     ! = mumol PAR m-2 s-1 Average PAR during the photoperiod
real :: PARINT    ! = mol PAR m-2 d-1 PAR interception
real :: TRANRF    ! = Transpiration realisation factor

contains

! Calculate DTRINT,PARAV,PARINT = light interception variables
Subroutine Light(DAYL,DTR,LAI,PAR)
  real :: DAYL,DTR,LAI,PAR
  if (DAYL > 0) then
    PARAV = PAR * (1E6/(24*3600)) / DAYL
  else
    PARAV = 0.
  end if
  PARINT = PAR * (1 - exp(-1.0*KLAI*LAI))  ! PAR extinction, Simon renamed K to KLAI
  DTRINT = DTR * (1 - exp(-0.75*KLAI*LAI)) ! GR has different extinction, Simon renamed K to KLAI
end Subroutine Light

! Calculate EVAP,TRAN,TRANRF
! See equations in Marcel van Oijen and Peter Leffelaar Crop Ecology 2010
! Chapter 10(B): Lintul-2: water limited crop growth
Subroutine EVAPTRTRF(Fdepth,PEVAP,PTRAN,CRT,ROOTD,WAL,WCL, EVAP,TRAN)
  real :: Fdepth, PEVAP, PTRAN, CRT,ROOTD, WAL,WCL,  EVAP, TRAN
  real :: AVAILF, FR, WAAD, WCCR
!  real :: WCL ! Simon use previously calculated WCL
!  if (Fdepth < ROOTD) then
!    WCL = WAL*0.001 / (ROOTD-Fdepth)
!  else
!    WCL = 0
!  end if                                                       ! (m3 m-3)
  EVAP = PEVAP * max(0., min(1., (WCL-WCAD)/(WCFC-WCAD) ))      ! = mm d-1 Evaporation of water from soil surface
  WCCR = WCWP + (WCFC - WCWP) * max(0.0, PTRAN/(PTRAN+TRANCO))  ! = m3 m-3 Critical water content below which transpiration is reduced (Eqn 1)
!  WCCR = WCWP + (WCFC - WCWP) * max(0.0, PTRAN/(PTRAN+TRANCO*1000)) ! Simon made TRANCO 1000x smaller
!  WCCR = WCWP + (WCFC - WCWP) * max(0.0, PTRAN/(PTRAN+TRANCO*CRT)) ! Simon made this a function of CRT and removed lower bound
  if (WCL > WCCR) then                                          ! Transpiraiton reduction factor (Fig 4)
    FR = max(0., min(1., (WCST-WCL)/(WCST-WCWET) ))             ! Transpiration reduction in wet conditions
  else if (WCCR > WCWP) then                                    ! Gradual reduction due to dryness
    FR = max(0., min(1., (WCL-WCWP)/(WCCR-WCWP)  ))             ! Transpiration reduction in dry conditions
  else
    FR = 0.0                                                     ! Simon added this case explicitly instead of setting lower bound to WCCR
  end if
  TRAN = PTRAN * FR                                             ! = mm d-1 Transpiration reduction due to WCL
  WAAD = 1000. * WCAD * (ROOTD-Fdepth)                          ! = mm Water in non-frozen soil at air dryness
  if (EVAP+TRAN > 0.) then
    AVAILF = min( 1., ((WAL-WAAD)/DELT) / (EVAP+TRAN) )         ! = Reduction when near air dryness to prevent WAL<WAAD
  else
    AVAILF = 0
  end if
  EVAP = EVAP * AVAILF                                          ! (mm d-1)
  TRAN = TRAN * AVAILF                                          ! (mm d-1)
  if (PTRAN > 0.) then
    TRANRF = TRAN / PTRAN                                       ! (-) Water restriction on plant processes
  else
    TRANRF = 1                                                  ! (-) Assume no restriction when PTRAN==0
  end if
end Subroutine EVAPTRTRF

! Calculate root depth growth rate RROOTD,EXPLOR
Subroutine ROOTDG(Fdepth,ROOTD,WAL,WCL,FAGE, EXPLOR,RROOTD)
  real :: Fdepth,ROOTD,WAL,WCL,FAGE
  real :: EXPLOR,RROOTD
!  real :: WCL ! Simon use previously calculated WCL
!  if (Fdepth < ROOTD) then
!    WCL = WAL*0.001 / (ROOTD-Fdepth)
!  else
!    WCL = 0
!  end if                                                        ! (m3 m-3)
  if ( (ROOTD<ROOTDM) .and. (WCL>WCWP) ) then
     RROOTD = min( RRDMAX, (ROOTDM-ROOTD)/DELT ) ! = m d-1 Root depth growth rate (basically constant = RRDMAX)
  else
     RROOTD = 0.
  end if
  EXPLOR = 1000. * RROOTD * WCFC                 ! = mm d-1 Increased access to water by root depth growth
end Subroutine ROOTDG

end module resources
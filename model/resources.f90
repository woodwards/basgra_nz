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
Subroutine EVAPTRTRF(Fdepth,PEVAP,PTRAN,ROOTD,WAL,WCL, EVAP,TRAN)
  real :: Fdepth, PEVAP, PTRAN, ROOTD, WAL,WCL,  EVAP, TRAN
  real :: AVAILF, FR, WAAD, WCCR
!  real :: WCL
!  if (Fdepth < ROOTD) then
!    WCL = WAL*0.001 / (ROOTD-Fdepth)
!  else
!    WCL = 0
!  end if                                                        ! (m3 m-3)
  WAAD = 1000. * WCAD * (ROOTD-Fdepth)                          ! = mm Water in non-frozen soil at air dryness
  EVAP = PEVAP * max(0., min(1., (WCL-WCAD)/(WCFC-WCAD) ))      ! = mm d-1 Evaporation of water from soil surface
  WCCR = WCWP + max( 0.01, PTRAN/(PTRAN+TRANCO) * (WCFC-WCWP) ) ! = m3 m-3 Water content "critical"?? FIXME
  if (WCL > WCCR) then
    FR = max(0., min(1., (WCST-WCL)/(WCST-WCWET) ))             ! Transpiration realisation in wet conditions? FIXME
  else
    FR = max(0., min(1., (WCL-WCWP)/(WCCR-WCWP)  ))             ! Transpiration realisation in dry conditions? FIXME looks wrong
  end if                                                        !
  TRAN = PTRAN * FR                                             ! = mm d-1 Transpiration reduction due to WCL?
  if (EVAP+TRAN > 0.) then
    AVAILF = min( 1., ((WAL-WAAD)/DELT) / (EVAP+TRAN) )         ! = Availability of water for evapotranspiration
  else
    AVAILF = 0
  end if                                                        ! (mm mm-1)
  EVAP = EVAP * AVAILF                                          ! (mm d-1)
  TRAN = TRAN * AVAILF                                          ! (mm d-1)
  if (PTRAN > 0.) then
    TRANRF = TRAN / PTRAN                                       ! (-)
  else
    TRANRF = 1                                                  ! (-)
  end if
end Subroutine EVAPTRTRF

! Calculate root depth growth rate RROOTD,EXPLOR
Subroutine ROOTDG(Fdepth,ROOTD,WAL,WCL, EXPLOR,RROOTD)
  real :: Fdepth,ROOTD,WAL,WCL
  real :: EXPLOR,RROOTD
!  real :: WCL
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

module soil

use parameters_site
use parameters_plant

implicit none

! Soil variables
real :: FO2     ! = mol O2 mol-1 gas Soil oxygen as a fraction of total gas
real :: fPerm   ! = not used
real :: Tsurf   ! = soil surface temperature, and fPerm = not used
real :: WCL     ! = Liquid soil water content between frost depth and root depth

contains

! Calculate WCL = Liquid soil water content between frost depth and root depth
Subroutine SoilWaterContent(Fdepth,ROOTD,WAL)
  real :: Fdepth,ROOTD,WAL
  if (Fdepth < ROOTD) then
    WCL = WAL*0.001 / (ROOTD-Fdepth)
  else
    WCL = 0
  end if
end Subroutine SoilWaterContent

! Calculate Tsurf = soil surface temperature, and fPerm = not used
Subroutine Physics(DAVTMP,Fdepth,ROOTD,Sdepth,WAS, Frate)
  real :: DAVTMP,Fdepth,ROOTD,Sdepth,WAS
  real :: Frate
  if (Fdepth > 0.) then
    Tsurf = DAVTMP / (1. + 10. * (Sdepth / Fdepth) ) ! Temperature extinction under snow when soil is frozen? FIXME
    fPerm = 0. ! Not used
  else
    Tsurf = DAVTMP * exp(-KTSNOW*Sdepth) ! Temperature extinction under snow (KTSNOW = gamma ~ 65 m-1)
    fPerm = 1. ! Not used
  end if
  call Frozensoil(Fdepth,ROOTD,WAS, Frate)
end Subroutine Physics

   ! Calculate Frate = m d-1 Rate of increase of frost layer depth
   Subroutine FrozenSoil(Fdepth,ROOTD,WAS, Frate)
     real :: Fdepth,ROOTD,WAS
     real :: Frate
     real :: alpha, PFrate, WCeff
     ! Determining the amount of solid water that contributes in transportation of heat to surface 'WCeff'
     if (Fdepth > ROOTD) then
       WCeff = WCFC
     else if (Fdepth > 0.) then
       WCeff = (0.001*WAS) / Fdepth
     else
       WCeff = WCL
     end if
     ! Calculating potential frost rate 'PFrate'
!     if ((Fdepth == 0.).and.(Tsurf>0.)) then ! No soil frost present AND no frost starting
     if (((Fdepth == 0.).and.(Tsurf>0.)).or.(WCeff == 0.)) then ! No soil frost present AND no frost starting
       PFrate = 0.
     else
       alpha  = LAMBDAsoil / ( RHOwater * WCeff * LatentHeat )
       PFrate = sqrt( max(0.,Fdepth**2 - 2.*alpha*Tsurf) ) - Fdepth
     end if
     if ((PFrate >= 0.).and.(Fdepth > 0.).and.(Fdepth < ROOTD)) then
       Frate = PFrate * (0.001*WAS/Fdepth) / WCFC ! Soil frost increasing
     else if ((PFrate+Fdepth/DELT) < 0.) then
       Frate = -Fdepth / DELT                     ! Remaining soil frost thaws away
     else
       Frate = PFrate
     end if
   end Subroutine FrozenSoil

! Calculate DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS
Subroutine FRDRUNIR(EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS, &
                                               DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS)
  real :: EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS
  real :: DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS
  real :: INFILTOT,WAFC,WAST
  WAFC   = 1000. * WCFC * max(0.,(ROOTD-Fdepth))                      ! (mm)
  WAST   = 1000. * WCST * max(0.,(ROOTD-Fdepth))                      ! (mm)
  INFILTOT = INFIL + poolDrain
  if (Fdepth < ROOTD) then
    FREEZEL = max(0., min( WAL/DELT + (INFILTOT - EVAP - TRAN), &
                         (Frate/(ROOTD-Fdepth))*WAL))                 ! = mm d-1 Freezing of soil water
  else
    FREEZEL = 0.
  end if
  if ((Fdepth > 0.) .and. (Fdepth <= ROOTD)) then
    THAWS   = max(0.,min( WAS/DELT, -Frate*WAS/Fdepth ))              ! = mm d-1 Thawing of soil frost
  else
    THAWS   = 0.
  end if
  DRAIN  = max(0.,min( DRATE, (WAL-WAFC)/DELT + &
         (INFILTOT - EVAP - TRAN - FREEZEL + THAWS) ))                ! = mm d-1 Drainage
  RUNOFF = max(0.,            (WAL-WAST)/DELT + &
         (INFILTOT - EVAP - TRAN - FREEZEL + THAWS - DRAIN) )         ! = mm d-1 Runoff
  IRRIG  = IRRIGF *  (        (WAFC-WAL)/DELT - &
         (INFILTOT - EVAP - TRAN - FREEZEL + THAWS - DRAIN - RUNOFF)) ! = mm d-1 Irrigation
end Subroutine FRDRUNIR

! Calculate FO2 = mol O2 mol-1 gas	Soil oxygen as a fraction of total gas
Subroutine O2status(O2,ROOTD)
  real :: O2,ROOTD
  FO2 = O2 / (ROOTD * FGAS * 1000./22.4) ! FGAS is a parameter
end Subroutine O2status

Subroutine O2fluxes(O2,PERMgas,ROOTD,RplantAer, O2IN,O2OUT)
  real :: O2,PERMgas,ROOTD,RplantAer
  real :: O2IN,O2OUT
  real :: O2MX
  O2OUT = RplantAer * KRTOTAER * 1./12. * 1.
  O2MX  = FO2MX * ROOTD * FGAS * 1000./22.4
  O2IN  = PERMgas * ( (O2MX-O2) + O2OUT*DELT )
end Subroutine O2fluxes

end module soil

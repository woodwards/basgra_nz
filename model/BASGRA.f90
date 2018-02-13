subroutine BASGRA(PARAMS,MATRIX_WEATHER,DAYS_HARVEST,NDAYS,NOUT,y)
!-------------------------------------------------------------------------------
! This is the BASic GRAss model originally written in MATLAB/Simulink by Marcel
! van Oijen, Mats Hoglind, Stig Morten Thorsen and Ad Schapendonk.
! 2011-07-13: Translation to FORTRAN by David Cameron and Marcel van Oijen.
! 2014-03-17: Extra category of tillers added
! 2014-04-03: Vernalization added
! 2014-04-03: Lower limit of temperature-driven leaf senescence no longer zero
!-------------------------------------------------------------------------------

! Allows access to all public objects in the other modules
use parameters_site
use parameters_plant
use environment
use resources
use soil
use plant

implicit none

! Define model inputs
integer               :: NDAYS, NOUT
integer, dimension(100,3) :: DAYS_HARVEST     ! Simon added third column (= percent harvested)
integer, parameter    :: NPAR     = 88        ! Note: NPAR is also hardwired in set_params.f90
! BASGRA handles two types of weather files with different data columns
#ifdef weathergen
  integer, parameter  :: NWEATHER =  7
#else
  integer, parameter  :: NWEATHER =  8
#endif
real                  :: PARAMS(NPAR)
real                  :: MATRIX_WEATHER(NMAXDAYS,NWEATHER)
real                  :: y(NDAYS,NOUT)

! Define time variables
integer               :: day, doy, i, year

! Define state variables
real :: CLV, CLVD, YIELD, CRES, CRT, CST, CSTUB, DRYSTOR, Fdepth, LAI, LT50, O2, PHEN, AGE
real ::            YIELD_LAST
real :: ROOTD, Sdepth, TILG1, TILG2, TILV, TANAER, WAL, WAPL, WAPS, WAS, WETSTOR
integer :: VERN
real :: VERND, DVERND

! Define intermediate and rate variables
real :: DeHardRate, DLAI, DLV, DLVD, DPHEN, DRAIN, DRT, DSTUB, dTANAER, DTILV, EVAP, EXPLOR, FAGE
real :: Frate, FREEZEL, FREEZEPL, GLAI, GLV, GPHEN, GRES, GRT, GST, GSTUB, GTILV, HardRate
real :: HARVFR, HARVLA, HARVLV, HARVPH, HARVRE, HARVST, HARVTILG2, INFIL, IRRIG, O2IN
real :: O2OUT, PackMelt, poolDrain, poolInfil, Psnow, reFreeze, RGRTV
real :: RGRTVG1, RROOTD, RUNOFF, SnowMelt, THAWPS, THAWS, TILVG1, TILG1G2, TRAN, Wremain
integer :: HARV

! Extract parameters
call set_params(PARAMS)

! Extract calendar and weather data
YEARI  = MATRIX_WEATHER(:,1)
DOYI   = MATRIX_WEATHER(:,2)
GRI    = MATRIX_WEATHER(:,3)
TMMNI  = MATRIX_WEATHER(:,4)
TMMXI  = MATRIX_WEATHER(:,5)
#ifdef weathergen
  RAINI = MATRIX_WEATHER(:,6)
  PETI  = MATRIX_WEATHER(:,7)
#else
  VPI   = MATRIX_WEATHER(:,6)
  RAINI = MATRIX_WEATHER(:,7)
  WNI   = MATRIX_WEATHER(:,8)
#endif

! Initialise state variables
AGE     = 0.0
CLV     = CLVI
CLVD    = CLVDI
CRES    = CRESI
CRT     = CRTI
CST     = CSTI
CSTUB   = CSTUBI
DAYL    = 0.5        ! Simon used to initialise YDAYL
DRYSTOR = DRYSTORI
Fdepth  = FdepthI
LAI     = max(SLAMAX * FSLAMIN * CLV, min(LAII, SLAMAX * CLV)) ! Simon constrain initial LAI (note SLAMAX is in gC units)
LT50    = LT50I
O2      = FGAS * ROOTDM * FO2MX * 1000./22.4
PHEN    = PHENI
ROOTD   = ROOTDM
Sdepth  = SDEPTHI
TANAER  = TANAERI
TILG1   = TILTOTI *       FRTILGI *    FRTILGG1I
TILG2   = TILTOTI *       FRTILGI * (1-FRTILGG1I)
TILV    = TILTOTI * (1. - FRTILGI)
VERN    = 0
VERND   = 0.0        ! Simon initialise count of cold days
YIELD   = YIELDI
YIELD_LAST = YIELDI
WAL     = 1000. * ROOTDM * WCI
WAPL    = WAPLI
WAPS    = WAPSI
WAS     = WASI
WETSTOR = WETSTORI

! Loop through days
do day = 1, NDAYS

  ! Calculate intermediate and rate variables
  !    SUBROUTINE      INPUTS                          OUTPUTS

  ! Environment
  call set_weather_day(day,DRYSTOR,                    year,doy) ! set weather for the day
  call SoilWaterContent(Fdepth,ROOTD,WAL)                        ! calculate WCL
  call Physics        (DAVTMP,Fdepth,ROOTD,Sdepth,WAS, Frate)    ! calculate Tsurf
  call MicroClimate   (doy,DRYSTOR,Fdepth,Frate,LAI,Sdepth,Tsurf,WAPL,WAPS,WETSTOR, &
                                                       FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil, &
                                                       pSnow,reFreeze,SnowMelt,THAWPS,wRemain)
                                                                 ! calculate water, snow and ice
  call DDAYL          (doy)                                      ! calculate DAYL
#ifdef weathergen
  call PEVAPINPUT     (LAI)                                      ! calculate PEVAP, PTRAN
#else
  call PENMAN         (LAI)                                      ! calculate PEVAP, PTRAN
#endif

  ! Resources
  call Light          (DAYL,DTR,LAI,PAR)                             ! calculate Light interception DTRINT,PARINT,PARAV
  call EVAPTRTRF      (Fdepth,PEVAP,PTRAN,CRT,ROOTD,WAL,WCL,EVAP,TRAN)   ! calculate EVAP,TRAN,TRANRF
  call ROOTDG         (Fdepth,ROOTD,WAL,WCL,FAGE,      EXPLOR,RROOTD)! calculate root depth increase rate RROOTD,EXPLOR

  ! Soil
  call FRDRUNIR       (EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS, &
                                                       DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS)
                                              ! calculate water movement etc DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS
  call O2status       (O2,ROOTD)              ! calculate FO2

  ! Plant
  call Harvest        (CLV,CRES,CST,year,doy,DAYS_HARVEST,LAI,PHEN,TILG2,TILG1,TILV, &
                                                       GSTUB,HARVLA,HARVLV,HARVPH,HARVRE,HARVST,HARVTILG2,HARVFR,HARV)
  call Biomass        (CLV,CRES,CST,CSTUB)
  call Phenology      (DAYL,PHEN,AGE,                  DPHEN,GPHEN,HARVPH,FAGE)
  call Vernalisation  (DAYL,PHEN,YDAYL,TMMN,TMMX,VERN,VERND,DVERND)       ! Simon vernalisation function
  call SLA
  call LUECO2TM       (PARAV)
  call HardeningSink  (CLV,DAYL,doy,LT50,Tsurf)
  call Growth         (CLV,CRES,CST,PARINT,TILG2,TILG1,TILV,TRANRF,GLV,GRES,GRT,GST)
  call PlantRespiration(FO2,RESPHARD)
  call Senescence     (CLV,CRT,CSTUB,doy,LAI,LT50,PERMgas,TANAER,TILV,Tsurf, &
                                                       DeHardRate,DLAI,DLV,DRT,DSTUB,dTANAER,DTILV,HardRate)
  call Decomposition  (CLVD,DAVTMP,WCL,                DLVD,RDLVD)    ! Simon decomposition function
  call Tillering       (DAYL,GLV,LAI,TILV,TILG1,TRANRF,Tsurf,VERN,FAGE, &
                                                       GLAI,GTILV,TILVG1,TILG1G2)
  ! Soil 2
  call O2fluxes       (O2,PERMgas,ROOTD,RplantAer,     O2IN,O2OUT)

  !================
  ! Outputs
  !================
  y(day, 1) = year + (doy-0.5)/366 ! "Time" = Decimal year (approximation)
  y(day, 2) = year
  y(day, 3) = doy
  y(day, 4) = DAVTMP

  y(day, 5) = CLV
  y(day, 6) = CLVD
  y(day, 7) = TRANRF
  y(day, 8) = CRES
  y(day, 9) = CRT
  y(day,10) = CST
  y(day,11) = CSTUB
  y(day,12) = VERND        ! (Simon changed)
  y(day,13) = PHOT         ! (Simon changed)
  y(day,14) = LAI
  y(day,15) = RESMOB       ! (Simon changed)
  y(day,16) = RAIN         ! mm Daily rainfall (Simon)
  y(day,17) = PHEN
  y(day,18) = ROOTD        ! m Root depth
  y(day,19) = DAYL         ! (Simon changed)
  y(day,20) = TILG2        ! (Simon changed)
  y(day,21) = TILG1        ! (Simon changed)
  y(day,22) = TILV
  y(day,23) = WAL          ! mm Soil water amount liquid
  y(day,24) = WCL * 100.0  ! WCL = volumetric water content pecentage (Simon changed)
  y(day,25) = DAYLGE       ! (Simon changed)
  y(day,26) = RDLVD        ! (Simon changed)
  y(day,27) = HARVFR * HARV! (Simon changed)

  ! Extra derived variables for calibration
  y(day,28) = (CLV+CST+CSTUB)/0.45 + CRES/0.40 + CLVD/0.45     ! "DM"  = Aboveground dry matter in g m-2 (Simon included CLVD)
  y(day,29) = (CRES/0.40) / ((CLV+CST+CSTUB)/0.45 + CRES/0.40) ! "RES" = Reserves in gDM gDM-1 aboveground green matter
  y(day,30) = LERG                               ! = m d-1 Leaf elongation rate per leaf for generative tillers
  y(day,31) = NELLVG                             ! = tiller-1 Number of growing leaves per elongating tiller
  y(day,32) = RLEAF                              ! = leaves tiller-1 d-1 Leaf appearance rate per tiller
  y(day,33) = LAI / CLV                          ! SLA     = m2 leaf area gC-1 dry matter vegetative tillers (RES not included?) Note in gC units
!  y(day,33) = LAI / (CLV/0.45)                   ! "SLA"     = m2 leaf area g-1 dry matter vegetative tillers (RES not included?) Note SLANEW internally is gC units
  y(day,34) = TILG1 + TILG2 + TILV               ! "TILTOT"  = Total tiller number in # m-2
  y(day,35) = (TILG1+TILG2) / (TILG1+TILG2+TILV) ! "FRTILG"  = Fraction of tillers that is generative
  y(day,36) =  TILG1        / (TILG1+TILG2+TILV) ! "FRTILG1" = Fraction of tillers that is in TILG1
  y(day,37) =        TILG2  / (TILG1+TILG2+TILV) ! "FRTILG2" = Fraction of tillers that is in TILG2
  y(day,38) = RDRT                               ! = d-1 Relative leaf death rate due to high temperature
  y(day,39) = VERN                               ! = Vernalisation flag

  ! Simon added additional output variables
  y(day,40) = DRAIN
  y(day,41) = RUNOFF
  y(day,42) = EVAP
  y(day,43) = TRAN
  y(day,44) = PARINT / PAR                       ! = Percentage light interception
  y(day,45) = RESPHARD                                ! Simon put variables here for debugging

  ! Update state variables
  AGE     = AGE     + 1.0
  CLV     = CLV     + GLV   - DLV    - HARVLV
  CLVD    = CLVD            + DLV             - DLVD       ! Simon included decomposition of dead material
  YIELD   = (HARVLV + HARVST) / 0.45 + HARVRE / 0.40       ! Simon separated HARVST and GSTUB
  if (YIELD>0) YIELD_LAST = YIELD
  CRES    = CRES    + GRES  - RESMOB - HARVRE
  CRT     = CRT     + GRT   - DRT
  CST     = CST     + GST            - HARVST - GSTUB      ! Simon separated HARVST and GSTUB
  CSTUB   = CSTUB   + GSTUB - DSTUB
  DRYSTOR = DRYSTOR + reFreeze + Psnow - SnowMelt
  Fdepth  = Fdepth  + Frate
  LAI     = LAI     + GLAI - DLAI   - HARVLA
  LT50    = LT50    + DeHardRate - HardRate
  O2      = O2      + O2IN - O2OUT
  PHEN    = max(0.0, min(1.0, PHEN + GPHEN - DPHEN - HARVPH))
  ROOTD   = ROOTD   + RROOTD
  Sdepth  = Sdepth  + Psnow/RHOnewSnow - PackMelt
  TANAER  = TANAER  + dTANAER
  TILG1   = TILG1           + TILVG1 - TILG1G2
  TILG2   = TILG2                    + TILG1G2 - HARVTILG2
  TILV    = TILV    + GTILV - TILVG1           - DTILV
  VERN    = VERN
  VERND   = VERND   + DVERND
  WAL     = WAL  + THAWS  - FREEZEL  + poolDrain + INFIL + EXPLOR + IRRIG - DRAIN - RUNOFF - EVAP - TRAN
  WAPL    = WAPL + THAWPS - FREEZEPL + poolInfil - poolDrain
  WAPS    = WAPS - THAWPS + FREEZEPL
  WAS     = WAS  - THAWS  + FREEZEL
  WETSTOR = WETSTOR + Wremain - WETSTOR

enddo

end

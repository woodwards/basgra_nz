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
integer, dimension(100,3) :: DAYS_HARVEST     ! Simon added third column (percent leaf removed)
integer, parameter    :: NPAR     = 81
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
real :: CLV, CLVD, YIELD, CRES, CRT, CST, CSTUB, DRYSTOR, Fdepth, LAI, LT50, O2, PHEN
real ::            YIELD_LAST
real :: ROOTD, Sdepth, TILG1, TILG2, TILV, TANAER, WAL, WAPL, WAPS, WAS, WETSTOR
integer :: VERN, NEWVERN

! Define intermediate and rate variables
real :: DeHardRate, DLAI, DLV, DLVD, DPHEN, DRAIN, DRT, DSTUB, dTANAER, DTILV, EVAP, EXPLOR ! Simon added DLVD
real :: Frate, FREEZEL, FREEZEPL, GLAI, GLV, GPHEN, GRES, GRT, GST, GSTUB, GTILV, HardRate
real :: HARVFR, HARVLA, HARVLV, HARVPH, HARVRE, HARVST, HARVTILG2, INFIL, IRRIG, O2IN ! Simon exposed HARVFR
real :: O2OUT, PackMelt, poolDrain, poolInfil, Psnow, reFreeze, RESMOB, RGRTV
real :: RGRTVG1, RROOTD, RUNOFF, SnowMelt, THAWPS, THAWS, TILVG1, TILG1G2, TRAN, Wremain
integer :: HARV

! Extract parameters
call set_params(PARAMS)  ! Note: NPAR is hardwired in set_params()

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
CLV     = CLVI
CLVD    = CLVDI
YIELD   = YIELDI
YIELD_LAST = YIELDI
CRES    = CRESI
CRT     = CRTI
CST     = CSTI
CSTUB   = CSTUBI
DRYSTOR = DRYSTORI
Fdepth  = FdepthI
LAI     = LAII
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
  call EVAPTRTRF      (Fdepth,PEVAP,PTRAN,ROOTD,WAL,WCL,EVAP,TRAN)   ! calculate EVAP,TRAN,TRANRF (Simon passed WCL)
  call ROOTDG         (Fdepth,ROOTD,WAL,WCL,           EXPLOR,RROOTD)! calculate root depth increase rate RROOTD,EXPLOR (Simon passed WCL)

  ! Soil
  call FRDRUNIR       (EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS, &
                                                       DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS)
                                              ! calculate water movement etc DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS
  call O2status       (O2,ROOTD)              ! calculate FO2

  ! Plant
  call Harvest        (CLV,CRES,CST,year,doy,DAYS_HARVEST,LAI,PHEN,TILG2,TILG1,TILV, &     ! Simon added TILG1
                                                       GSTUB,HARVLA,HARVLV,HARVPH,HARVRE,HARVST,HARVTILG2,HARVFR,HARV) ! Simon added HARVFR
  call Biomass        (CLV,CRES,CST)
  call Phenology      (DAYL,PHEN,                      DPHEN,GPHEN)
  call Vernalisation  (DAYL,DAVTMP,VERN,               NEWVERN)       ! Simon vernalisation function
  call Foliage1
  call LUECO2TM       (PARAV)
  call HardeningSink  (CLV,DAYL,doy,LT50,Tsurf)
  call Growth         (CLV,CRES,CST,PARINT,TILG2,TILV,TRANRF, &
                                                       GLV,GRES,GRT,GST,RESMOB)
  call PlantRespiration(FO2,RESPHARD)
  call Senescence     (CLV,CRT,CSTUB,doy,LAI,LT50,PERMgas,TANAER,TILV,Tsurf, &
                                                       DeHardRate,DLAI,DLV,DRT,DSTUB,dTANAER,DTILV,HardRate)
  call Decomposition  (CLVD,DAVTMP,WCL,                DLVD,DELTA)    ! Simon decomposition function
  call Foliage2       (DAYL,GLV,LAI,TILV,TILG1,TRANRF,Tsurf,VERN, &
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
  y(day, 7) = YIELD_LAST
  y(day, 8) = CRES
  y(day, 9) = CRT
  y(day,10) = CST
  y(day,11) = CSTUB
  y(day,12) = DRYSTOR      ! mm Snow amount as SWE (Soil Water Equivalent)
  y(day,13) = Fdepth       ! m Soil frost layer depth
  y(day,14) = LAI
  y(day,15) = LT50         ! deg C Temperature that kills half the plants in a day
  y(day,16) = RAIN         ! mm Daily rainfall (Simon)
  y(day,17) = PHEN
  y(day,18) = ROOTD        ! m Root depth
  y(day,19) = Sdepth       ! m Snow depth
  y(day,20) = TILG2        ! (Simon changed)
  y(day,21) = TILG1        ! (Simon changed)
  y(day,22) = TILV
  y(day,23) = WAL          ! mm Soil water amount liquid
  y(day,24) = WCL * 100    ! WCL = volumetric water content pecentage (Simon changed)
  y(day,25) = WAPS         ! mm Pool water amount solid (ice)
  y(day,26) = DELTA        ! (Simon changed)
  y(day,27) = HARVFR * 100 * HARV ! (Simon changed)

  ! Extra derived variables for calibration
  y(day,28) = (CLV+CST+CSTUB)/0.45 + CRES/0.40   ! "DM"      = Aboveground dry matter in g m-2 (includes CSTUB but excludes CLVD)
  y(day,29) = (CRES/0.40) / y(day,28)            ! "RES"     = Reserves in g g-1 aboveground dry matter
  y(day,30) = LERG                               ! = m d-1 Leaf elongation rate per leaf for generative tillers
  y(day,31) = NELLVG                             ! = tiller-1 Number of growing leaves per elongating tiller
  y(day,32) = RLEAF                              ! = leaves tiller-1 d-1 Leaf appearance rate per tiller
  y(day,33) = LAI / (CLV/0.45)                   ! "SLA"     = m2 leaf area g-1 dry matter vegetative tillers
  y(day,34) = TILG1 + TILG2 + TILV               ! "TILTOT"  = Total tiller number in # m-2
  y(day,35) = (TILG1+TILG2) / (TILG1+TILG2+TILV) ! "FRTILG"  = Fraction of tillers that is generative
  y(day,36) =  TILG1        / (TILG1+TILG2+TILV) ! "FRTILG1" = Fraction of tillers that is in TILG1
  y(day,37) =        TILG2  / (TILG1+TILG2+TILV) ! "FRTILG2" = Fraction of tillers that is in TILG2
  y(day,38) = RDRT                               ! = d-1 Relative leaf death rate due to high temperature
  y(day,39) = VERN                               ! = Vernalisation flag

  ! Extra variables added by Simon
  y(day,40) = DRAIN
  y(day,41) = RUNOFF
  y(day,42) = EVAP
  y(day,43) = TRAN
  y(day,44) = PARINT / PAR * 100                 ! = Percentage light interception

  ! Update state variables
  CLV     = CLV     + GLV   - DLV    - HARVLV
  CLVD    = CLVD            + DLV    - DLVD                ! Simon, no decomposition of dead material FIXME
  YIELD   = (HARVLV + HARVST) / 0.45 + HARVRE / 0.40       ! Simon, separated HARVST and GSTUB
  if (YIELD>0) YIELD_LAST = YIELD
  CRES    = CRES    + GRES  - RESMOB - HARVRE
  CRT     = CRT     + GRT   - DRT
  CST     = CST     + GST            - HARVST - GSTUB      ! Simon, separated HARVST and GSTUB
  CSTUB   = CSTUB   + GSTUB - DSTUB
  DRYSTOR = DRYSTOR + reFreeze + Psnow - SnowMelt
  Fdepth  = Fdepth  + Frate
  LAI     = LAI     + GLAI - DLAI   - HARVLA
  LT50    = LT50    + DeHardRate - HardRate
  O2      = O2      + O2IN - O2OUT
  PHEN    = min(1., PHEN + GPHEN - DPHEN - HARVPH)
  ROOTD   = ROOTD   + RROOTD
  Sdepth  = Sdepth  + Psnow/RHOnewSnow - PackMelt
  TANAER  = TANAER  + dTANAER
  TILG1   = TILG1           + TILVG1 - TILG1G2
  TILG2   = TILG2                    + TILG1G2 - HARVTILG2
  TILV    = TILV    + GTILV - TILVG1           - DTILV
  VERN    = NEWVERN
  WAL     = WAL  + THAWS  - FREEZEL  + poolDrain + INFIL +EXPLOR+IRRIG-DRAIN-RUNOFF-EVAP-TRAN
  WAPL    = WAPL + THAWPS - FREEZEPL + poolInfil - poolDrain
  WAPS    = WAPS - THAWPS + FREEZEPL
  WAS     = WAS  - THAWS  + FREEZEL
  WETSTOR = WETSTOR + Wremain - WETSTOR                    ! Simon this looks weird

enddo

end

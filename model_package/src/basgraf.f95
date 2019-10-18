module basgramodule
    use, intrinsic :: iso_c_binding

    implicit none
    private
    public :: BASGRA

contains

subroutine BASGRA(PARAMS,MATRIX_WEATHER,DAYS_HARVEST,NDAYS,NOUT,y) bind(C, name = "BASGRA_")
!-------------------------------------------------------------------------------
! This is the BASic GRAss model originally written in MATLAB/Simulink by Marcel
! van Oijen, Mats Hoglind, Stig Morten Thorsen and Ad Schapendonk.
! 2011-07-13: Translation to FORTRAN by David Cameron and Marcel van Oijen.
! 2014-03-17: Extra category of tillers added
! 2014-04-03: Vernalization added
! 2014-04-03: Lower limit of temperature-driven leaf senescence no longer zero
! 2018-08-01: Modified by Simon Woodward for New Zealand ryegrass simulations
! 2019-06-09: Added C wrapper to allow model to be compiled into an R package
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
integer(kind = c_int), intent(in), value            :: NDAYS
integer(kind = c_int), intent(in), value            :: NOUT
integer(kind = c_int), intent(in), dimension(100,3) :: DAYS_HARVEST ! Simon added third column (= pc harvested)
integer, parameter                                  :: NPAR     = 108 ! NPAR also hardwired in set_params.f90
! BASGRA handles two types of weather files with different data columns
#ifdef weathergen
  integer, parameter                                :: NWEATHER =  7
#else
  integer, parameter                                :: NWEATHER =  8
#endif
real(kind = c_double), intent(in), dimension(NPAR)              :: PARAMS
real(kind = c_double), intent(in), dimension(NMAXDAYS,NWEATHER) :: MATRIX_WEATHER
real(kind = c_double), intent(out), dimension(NDAYS,NOUT)       :: y

! Define time variables
integer               :: day, doy, i, year

! Define state variables
real :: CLV, CLVD, YIELD, CRES, CRT, CST, CSTUB, DRYSTOR, Fdepth, LAI, LT50, O2, PHEN, AGE
real :: ROOTD, Sdepth, TILG1, TILG2, TILV, TANAER, WAL, WAPL, WAPS, WAS, WETSTOR
!integer :: VERN
real :: VERN                                  ! Simon made VERN a continuous function of VERND
real :: VERND, DVERND, WALS, BASAL

! Define intermediate and rate variables
real :: DeHardRate, DLAI, DLV, DLVD, DPHEN, DRAIN, DRT, DSTUB, dTANAER, DTILV, EVAP, EXPLOR
real :: Frate, FREEZEL, FREEZEPL, GLAI, GLV, GPHEN, GRES, GRT, GST, GSTUB, GTILV, HardRate
real :: HARVFR, HARVFRIN, HARVLA, HARVLV, HARVLVD, HARVPH, HARVRE, HARVST, HARVTILG2, INFIL, IRRIG, O2IN
real :: O2OUT, PackMelt, poolDrain, poolInfil, Psnow, reFreeze, RGRTV, RDRHARV
real :: RGRTVG1, RROOTD, RUNOFF, SnowMelt, THAWPS, THAWS, TILVG1, TILG1G2, TRAN, Wremain, SP
integer :: HARV

! Extra output variables (Simon)
real :: Time, DM, RES, SLA, TILTOT, FRTILG, FRTILG1, FRTILG2, LINT, DEBUG, TSIZE

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

! Initialise harvest array index
HARVI   = 1
do while ( (DAYS_HARVEST(HARVI,1)<YEARI(1)) .or. ((DAYS_HARVEST(HARVI,1)==YEARI(1)).and.(DAYS_HARVEST(HARVI,2)<DOYI(1))) )
  HARVI = HARVI + 1
end do

! Extract parameters
call set_params(PARAMS)

! Initial value transformations, Simon moved to here
CLVI  = 10**LOG10CLVI
!CLVDI =
!CRESI = 10**LOG10CRESI
CRTI  = 10**LOG10CRTI
!LAII  = 10**LOG10LAII

! Soil water parameter scaling, Simon moved to here
WCAD  = FWCAD  * WCST
WCWP  = FWCWP  * WCST
WCFC  = FWCFC  * WCST
WCWET = FWCWET * WCST

! Initialise state variables
AGE     = 0.0
CLV     = CLVI
CLVD    = CLVI * 0.2                             ! Simon initially set equal to leaf mass * 0.2
CRES    = (FCOCRESMN*0.5+0.5) * COCRESMX * (CLVI + CSTI)   ! Simon start at av
CRT     = CRTI
CST     = CSTI
CSTUB   = CSTUBI                                 ! Currently constant 0
DAYL    = 0.5                                    ! Simon used to initialise YDAYL
DRYSTOR = DRYSTORI
Fdepth  = FdepthI
LAI     = (FSLAMIN*0.5+0.5) * SLAMAX * CLV       ! Simon start at av, LAI is defined over entire area
LT50    = LT50I
O2      = FGAS * ROOTDM * FO2MX * 1000./22.4
PHEN    = PHENI
Sdepth  = SDEPTHI
TANAER  = TANAERI
TILG1   = TILTOTI *       FRTILGI *    FRTILGG1I
TILG2   = TILTOTI *       FRTILGI * (1-FRTILGG1I)
TILV    = TILTOTI * (1. - FRTILGI)
BASAL   = BASALI
ROOTD   = ROOTDM * CRT/BASAL / (CRT/BASAL + KCRT)! Simon tied ROOTD to CRT like this, CRT is defined over entire area
!VERN    = 0.0
!VERND   = floor(VERNDI)                           ! Simon initialise count of cold days
!  if ((VERN==0).and.(VERND .ge. TVERND)) then ! copied from Vernalisation()
!	VERN = 1.0
! 	VERND = 0.0
!    DVERND = 0.0
!  end if
VERND   = VERNDI
VERN    = max(0.0, min(1.0, (VERND-TVERNDMN)/(TVERND-TVERNDMN))) ! FIXME does not include effect of new summer tillers
YIELD   = YIELDI
WAL     = 1000. * (ROOTDM - Fdepth) * WCFC        ! Simon set to WCFC
WALS    = min(WAL, 25.0)                          ! Simon added WALS rapid surface layer (see manual section 4.3)
WAPL    = WAPLI
WAPS    = WAPSI
WAS     = WASI
WETSTOR = WETSTORI

! Loop through days
do day = 1, NDAYS

  ! Calculate intermediate and rate variables (many variable and parameters are passed implicitly)
  !    SUBROUTINE      INPUTS                          OUTPUTS

  call set_weather_day(day,DRYSTOR,                    year,doy) ! set weather for the day, including DTR, PAR, which depend on DRYSTOR

  call Harvest        (CLV,CRES,CST,CSTUB,CLVD,year,doy,DAYS_HARVEST,LAI,PHEN,TILG2,TILG1,TILV, &
                                                       GSTUB,HARVLA,HARVLV,HARVLVD,HARVPH,HARVRE,HARVST, &
                                                       HARVTILG2,HARVFR,HARVFRIN,HARV,RDRHARV)
  LAI     = LAI     - HARVLA * (1 + RDRHARV)
  CLV     = CLV     - HARVLV * (1 + RDRHARV)
  CLVD    = CLVD    - HARVLVD     + (HARVLV + HARVRE) * RDRHARV
  CRES    = CRES    - HARVRE * (1 + RDRHARV)
  CST     = CST     - HARVST   - GSTUB
  CSTUB   = CSTUB              + GSTUB
  TILV    = TILV    - TILV * RDRHARV
  TILG1   = TILG1   - TILG1 * RDRHARV
  TILG2   = TILG2   - HARVTILG2
  if (TILG2 < 1.0) then
	TILG2 = 0.0                              ! Simon avoid roundoff error in TILG2
  end if
  TILTOT  = TILG1 + TILG2 + TILV
  PHEN    = PHEN    - HARVPH
  if (doy.eq.152) then                               ! Reset yield on 1 June
    YIELD = 0.0
  end if
  YIELD     = YIELD + ((HARVLV + HARVLVD + HARVST) / 0.45 + HARVRE / 0.40) * 10.0 / 1000.0 ! tDM ha-1 Simon cumulative harvest

  call SoilWaterContent(Fdepth,ROOTD,WAL,WALS)                   ! calculate WCL
  call Physics        (DAVTMP,Fdepth,ROOTD,Sdepth,WAS, Frate)    ! calculate Tsurf, Frate
  call MicroClimate   (doy,DRYSTOR,Fdepth,Frate,LAI,BASAL,Sdepth,Tsurf,WAPL,WAPS,WETSTOR, &
                                                       FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil, &
                                                       pSnow,reFreeze,SnowMelt,THAWPS,wRemain) ! calculate water, snow and ice
  call DDAYL          (doy)                                      ! calculate DAYL, DAYLMX
#ifdef weathergen
  call PEVAPINPUT     (LAI,BASAL)                                      ! calculate PEVAP, PTRAN, depend on LAY, RNINTC
#else
  call PENMAN         (LAI,BASAL)                                      ! calculate PEVAP, PTRAN, depend on LAY, RNINTC
#endif

  call Light          (DAYL,DTR,LAI,BASAL,PAR)                   ! calculate light interception DTRINT,PARINT,PARAV
  call EVAPTRTRF      (Fdepth,PEVAP,PTRAN,CRT,ROOTD,WAL,WCLM,WCL,EVAP,TRAN)! calculate EVAP,TRAN,TRANRF

  call FRDRUNIR       (EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS, &
                                                       DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS) ! calculate water movement etc DRAIN,FREEZEL,IRRIG,RUNOFF,THAWS
  call O2status       (O2,ROOTD)                                 ! calculate FO2

  call Vernalisation  (DAYL,PHEN,YDAYL,TMMN,TMMX,DAVTMP,Tsurf,VERN,VERND,DVERND) ! Simon calculate VERN,VERND,DVERND
  call Phenology      (DAYL,TILG2,PHEN,         DPHEN,GPHEN,HARVPH) ! calculate GPHEN, DPHEN, PHENRF, DAYLGE
  call Biomass        (AGE,CLV,CRES,CST,CSTUB)                   ! calculate RESNOR
  call CalcSLA                                                   ! calculate LERV,LERG,SLANEW
  call LUECO2TM       (PARAV,BASAL)                              ! calculate LUEMXQ
  call HardeningSink  (CLV,DAYL,doy,LT50,Tsurf)                  ! calculate RESPHARDSI
  call Growth         (CLV,CRES,CST,PARINT,TILG2,TILG1,TILV,TRANRF,AGE,LAI, GLV,GRES,GRT,GST) ! calculate assimilate partitioning
  call PlantRespiration(FO2,RESPHARD)                            ! calculate RplantAer
  call Senescence     (CLV,CRT,CSTUB,doy,LAI,PARBASE,BASAL,LT50,PERMgas,TRANRF,TANAER,TILV,Tsurf,AGE, &
                                                       DeHardRate,DLAI,DLV,DRT,DSTUB,dTANAER,DTILV,HardRate,RDRS,RDRW)
  call Decomposition  (CLVD,DAVTMP,WCLM,                DLVD,RDLVD)    ! Simon decomposition function
  call Tillering      (DAYL,GLV,LAI,BASAL,TILV,TILG1,TRANRF,Tsurf,VERN,AGE, &
                                                       GLAI,RGRTV,GTILV,TILVG1,TILG1G2)

  call ROOTDG         (Fdepth,ROOTD,WAL,WCL,CRT,GRT,DRT, EXPLOR,RROOTD)! calculate root depth increase rate RROOTD,EXPLOR
  call O2fluxes       (O2,PERMgas,ROOTD,RplantAer,     O2IN,O2OUT)


  !================
  ! Outputs
  !================

! structural carbohydrate is 45% carbon C6H12O5
! soluble carbohydrate is 40% carbon C6H12O6

  Time      = year + (doy-0.5)/366 ! "Time" = Decimal year (approximation)
  DM        = ((CLV+CST+CSTUB)/0.45 + CRES/0.40 + CLVD/0.45) * 10.0 ! "DM"  = Aboveground dry matter in kgDM ha-1 (Simon included CLVD, changed units)
  RES       = (CRES/0.40) / ((CLV+CST+CSTUB)/0.45 + CRES/0.40)      ! "RES" = Reserves in gDM gDM-1 aboveground green matter
  SLA       = LAI / CLV                          ! SLA     = m2 leaf area gC-1 dry matter vegetative tillers (Note units and RES not included)
  TSIZE     = (CLV+CST) / (TILG1+TILG2+TILV)     ! gC tillers-1 Average tiller size
  FRTILG    = (TILG1+TILG2) / (TILG1+TILG2+TILV) ! "FRTILG"  = Fraction of tillers that is generative
  FRTILG1   =  TILG1        / (TILG1+TILG2+TILV) ! "FRTILG1" = Fraction of tillers that is in TILG1
  FRTILG2   =        TILG2  / (TILG1+TILG2+TILV) ! "FRTILG2" = Fraction of tillers that is in TILG2
  LINT      = PARINT / PAR                       ! = Percentage light interception
  DEBUG     = LAI/BASAL                          ! Output any variable as "DEBUG" for debugging purposes

  ! a script checks that these variable names match what is expected in output_names.tsv (Simon)

  y(day, 1) = Time
  y(day, 2) = year
  y(day, 3) = doy
  y(day, 4) = DAVTMP

  y(day, 5) = CLV
  y(day, 6) = CLVD
  y(day, 7) = TRANRF * 100.0
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
  y(day,18) = LT50
  y(day,19) = DAYL         ! (Simon changed)
  y(day,20) = TILG2        ! (Simon changed)
  y(day,21) = TILG1        ! (Simon changed)
  y(day,22) = TILV
  y(day,23) = WAL          ! mm Soil water amount liquid
  y(day,24) = WCLM * 100.0 ! Soil moisture to ROOTDM (Simon changed)
  y(day,25) = DAYLGE       ! (Simon changed)
  y(day,26) = RDLVD        ! (Simon changed)
  y(day,27) = HARVFR * HARV! (Simon changed)

  ! Extra derived variables for calibration
  y(day,28) = DM
  y(day,29) = RES
  y(day,30) = LERG                               ! = m d-1 Leaf elongation rate per leaf for generative tillers
  y(day,31) = PHENRF                             ! Phenology effect
  y(day,32) = RLEAF                              ! = leaves tiller-1 d-1 Leaf appearance rate per tiller
  y(day,33) = SLA
  y(day,34) = TILTOT
  y(day,35) = RGRTV
  y(day,36) = RDRTIL
  y(day,37) = GRT
  y(day,38) = RDRL                               ! = d-1 Relative leaf death rate
  y(day,39) = VERN * 100.0                       ! = Vernalisation degree

  ! Simon added additional output variables
  y(day,40) = DRAIN
  y(day,41) = RUNOFF
  y(day,42) = EVAP
  y(day,43) = TRAN
  y(day,44) = LINT
  y(day,45) = DEBUG
  y(day,46) = ROOTD
  y(day,47) = TSIZE
  y(day,48) = LERV
  y(day,49) = WCL * 100.0
  y(day,50) = HARVFRIN * HARV
  y(day,51) = SLANEW
  y(day,52) = YIELD
  y(day,53) = BASAL * 100.0
  y(day,54) = GTILV
  y(day,55) = DTILV
  y(day,56) = FS

  ! Update state variables
  AGE     = AGE     + 1.0
  CLV     = CLV     + GLV   - DLV
  CLVD    = CLVD            + DLV             - DLVD       ! Simon included decomposition of dead material
  CRES    = CRES    + GRES  - RESMOB                       ! Simon modified harvest logic
  CRT     = CRT     + GRT   - DRT
  CST     = CST     + GST
  CSTUB   = CSTUB   - DSTUB
  DRYSTOR = DRYSTOR + reFreeze + Psnow - SnowMelt
  Fdepth  = Fdepth  + Frate
  LAI     = LAI     + GLAI - DLAI
  LT50    = LT50    + DeHardRate - HardRate
  O2      = O2      + O2IN - O2OUT
  PHEN    = PHEN    + GPHEN - DPHEN
  Sdepth  = Sdepth  + Psnow/RHOnewSnow - PackMelt
  TANAER  = TANAER  + dTANAER
  TILV    = TILV    + GTILV - TILVG1           - DTILV
  TILG1   = TILG1           + TILVG1 - TILG1G2
  TILG2   = TILG2                    + TILG1G2
  TILTOT  = TILG1 + TILG2 + TILV                           ! "TILTOT"  = Total tiller number in # m-2
!  BASAL   = BASAL * (1 - ABASAL) + TILTOT / (TILTOT + KBASAL) * ABASAL   ! Simon model grass basal area
!  BASAL   = BASAL * (1 - ABASAL) + min(1.0, TILTOT / KBASAL) * ABASAL   ! Simon model grass basal area
  BASAL   = BASAL * (1 - ABASAL) + min(1.0, LAI / KBASAL) * ABASAL   ! Simon model grass basal area
!  BASAL   = BASAL * (1 - ABASAL) + min(1.0, (1-exp(-KLAI*LAI)) / (1-exp(-KLAI*KBASAL)) ) * ABASAL   ! Simon model grass basal area
!  ROOTD   = ROOTD   + RROOTD                              ! Simon tied ROOTD to CRT
  ROOTD   = ROOTDM * CRT/BASAL / (CRT/BASAL + KCRT)                    ! Simon tied ROOTD to CRT like this
  VERND   = VERND   + DVERND
!  VERN    = VERN
  ! Simon treat VERN as a dynamic variable to capture effect of new summer tillers
  if (TILV>0) then
	VERN    = min(1.0, VERN + max(0.0, (VERND       -TVERNDMN)/(TVERND-TVERNDMN)) &
							- max(0.0, (VERND-DVERND-TVERNDMN)/(TVERND-TVERNDMN)) &
	                        - VERN * GTILV / TILV)
  else
	VERN    = 0.
  end if
  WAL     = WAL  + THAWS  - FREEZEL  + poolDrain + INFIL + EXPLOR + IRRIG - DRAIN - RUNOFF - EVAP - TRAN
  WALS    = max(0.0, min(25.0, WALS + THAWS - FREEZEL  + poolDrain + INFIL + IRRIG - DRAIN - RUNOFF - EVAP - TRAN)) ! Simon added WALS rapid surface pool
  WAPL    = WAPL + THAWPS - FREEZEPL + poolInfil - poolDrain
  WAPS    = WAPS - THAWPS + FREEZEPL
  WAS     = WAS  - THAWS  + FREEZEL
  WETSTOR = WETSTOR + Wremain - WETSTOR

enddo

end

end module basgramodule

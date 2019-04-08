module environment

use parameters_site
use parameters_plant

implicit none

! Environment variables
integer, parameter :: NMAXDAYS = 10000
real :: GR, TMMN, TMMX, VP, WN
real :: YEARI(NMAXDAYS), DOYI(NMAXDAYS) , RAINI(NMAXDAYS), GRI(NMAXDAYS)
real :: TMMNI(NMAXDAYS), TMMXI(NMAXDAYS), VPI(NMAXDAYS)  , WNI(NMAXDAYS)
#ifdef weathergen
real :: PETI(NMAXDAYS)
#endif
real :: DAVTMP,DAYL,YDAYL,DAYLMX,DTR,PAR,PERMgas,PEVAP,poolRUNOFF,PTRAN,pWater,RAIN,RNINTC
real :: runOn,StayWet,WmaxStore,Wsupply
#ifdef weathergen
real :: PET
#endif

contains

! Set all time and weather variables for day
#ifdef weathergen
  Subroutine set_weather_day(day,DRYSTOR, year,doy)
    integer :: day, doy, year
    real    :: DRYSTOR
    year   = YEARI(day) ! day of the year (d)
    doy    = DOYI(day)  ! day of the year (d)
    RAIN   = RAINI(day) ! precipitation (mm d-1)
    GR     = GRI(day)   ! irradiation (MJ m-2 d-1)
    TMMN   = TMMNI(day) ! minimum temperature (degrees Celsius)
    TMMX   = TMMXI(day) ! maximum temperature (degrees Celsius)
    DAVTMP = (TMMN + TMMX)/2.0         ! daily average temperature
    DTR    = GR * exp(-KSNOW*DRYSTOR)  ! MJ GR m-2 d-1 Daily global radiation on leaves
    PAR    = 0.5*4.56*DTR              ! mol PAR m-2 d-1 Daily photosynthetically active radiation
    PET    = PETI(day)                 ! mm d-1 Daily potential evapotranspiration
  end Subroutine set_weather_day
#else
  Subroutine set_weather_day(day,DRYSTOR, year,doy)
    integer :: day, doy, year
    real    :: DRYSTOR
    year   = YEARI(day) ! day of the year (d)
    doy    = DOYI(day)  ! day of the year (d)
    RAIN   = RAINI(day) ! precipitation (mm d-1)
    GR     = GRI(day)   ! irradiation (MJ m-2 d-1)
    TMMN   = TMMNI(day) ! minimum (or average) temperature (degrees Celsius)
    TMMX   = TMMXI(day) ! maximum (or average) temperature (degrees Celsius)
    VP     = VPI(day)   ! vapour pressure (kPa)
    WN     = WNI(day)   ! mean wind speed (m s-1)
    DAVTMP = (TMMN + TMMX)/2.0
    DTR    = GR * exp(-KSNOW*DRYSTOR)
    PAR    = 0.5*4.56*DTR
  end Subroutine set_weather_day
#endif

Subroutine MicroClimate(doy,DRYSTOR,Fdepth,Frate,LAI,BASAL,Sdepth,Tsurf,WAPL,WAPS,WETSTOR, &
          FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil,pSnow,reFreeze,SnowMelt,THAWPS,wRemain)
  integer :: doy
  real :: DRYSTOR,Fdepth,Frate,LAI,BASAL,Sdepth,Tsurf,WAPL,WAPS,WETSTOR
  real :: FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil,pSnow,reFreeze,SnowMelt,THAWPS,wRemain
  call RainSnowSurfacePool(doy,DRYSTOR,Fdepth,Frate,LAI,BASAL,Sdepth,Tsurf,WAPL,WAPS,WETSTOR, &
       FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil,pSnow,reFreeze,SnowMelt,THAWPS,Wremain)
  if (WAPS == 0.) then
    PERMgas = 1.                ! Permeable to gas if no pool ice
  else
    PERMgas = 0.
  end if
end Subroutine MicroClimate

   ! See equation in Marcel van Oijen and Peter Leffelaar Crop Ecology 2010
   Subroutine RainSnowSurfacePool(doy,DRYSTOR,Fdepth,Frate,LAI,BASAL,Sdepth,Tsurf,WAPL,WAPS,WETSTOR, &
       FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil,pSnow,reFreeze,SnowMelt,THAWPS,Wremain)
     integer :: doy
     real :: DRYSTOR,Fdepth,Frate,LAI,BASAL,Sdepth,Tsurf,WAPL,WAPS,WETSTOR
     real :: FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil,pSnow,reFreeze,SnowMelt,THAWPS,Wremain
     real :: PINFIL
     call precForm(Psnow)
     call WaterSnow(doy,DRYSTOR,Psnow,Sdepth,WETSTOR, PackMelt,reFreeze,SnowMelt,Wremain)
     RNINTC = min( Wsupply, 0.25*LAI/BASAL ) ! Leaf can intercept 0.25 mm of water (Eqn 12)
     PINFIL = Wsupply - RNINTC         ! Not-intercepted fraction
     call INFILrunOn(Fdepth,PINFIL, INFIL)
     call SurfacePool(Fdepth,Frate,Tsurf,WAPL,WAPS, &
                                            FREEZEPL,poolDrain,poolInfil,THAWPS)
   end Subroutine RainSnowSurfacePool

      ! Determine form of precipitation, based on average daily Temp.
      Subroutine precForm(Psnow)
        real :: Psnow
        if (DAVTMP > TrainSnow) then ! TrainSnow is a parameter ~ 0.01 deg C?
          Pwater = RAIN
          Psnow  = 0.
        else
          Pwater = 0.
          Psnow  = RAIN
        end if
      end Subroutine precForm

      !
      Subroutine WaterSnow(doy,DRYSTOR,Psnow,Sdepth,WETSTOR, &
                                             PackMelt,reFreeze,SnowMelt,Wremain)
        integer :: doy
        real :: DRYSTOR,Psnow,Sdepth,WETSTOR
        real :: PackMelt,reFreeze,SnowMelt,Wremain
        real :: DENSITY
        call SnowMeltWmaxStore      (doy,DRYSTOR,             SnowMelt)
        call WETSTORdynamics        (WETSTOR,                 reFreeze)
        call LiquidWaterDistribution(SnowMelt,                Wremain)
        call SnowDensity            (DRYSTOR,Sdepth,WETSTOR,  DENSITY)
        call SnowDepthDecrease      (DENSITY,Sdepth,SnowMelt, PackMelt)
      end Subroutine WaterSnow

         !
         Subroutine SnowMeltWmaxStore(doy,DRYSTOR, SnowMelt)
           integer :: doy
           real :: DRYSTOR
           real :: SnowMelt
           real :: Melt
           Melt = Bias + Ampl * DAYL
           if (DAVTMP > TmeltFreeze) then ! TmeltFreeze is a parameter ~ 0.0 deg C?
             SnowMelt = max( 0., min( DRYSTOR/DELT, Melt*(DAVTMP-TmeltFreeze) ))
           else
             SnowMelt = 0.
           end if
           WmaxStore = DRYSTOR * SWret ! SWret is liquid water retention capacity of snow (~0.1 mm mm-1 d-1)
         end Subroutine SnowMeltWmaxStore

         Subroutine WETSTORdynamics(WETSTOR, reFreeze)
           real :: WETSTOR
           real :: reFreeze
           real :: reFreezeMax
           reFreezeMax = SWrf * (TmeltFreeze-DAVTMP)
           if ((WETSTOR>0).and.(DAVTMP<TmeltFreeze)) then
             reFreeze = min(WETSTOR/DELT,reFreezeMax)
           else
             reFreeze = 0.
           end if
           StayWet = WETSTOR/DELT - reFreeze
         end Subroutine WETSTORdynamics

         Subroutine LiquidWaterDistribution(SnowMelt, Wremain)
           real :: SnowMelt
           real :: Wremain
           real :: Wavail
           Wavail  = StayWet + SnowMelt + Pwater
           Wremain = min(Wavail,WmaxStore)
           Wsupply = Wavail - Wremain
         end Subroutine LiquidWaterDistribution

         Subroutine SnowDensity(DRYSTOR,Sdepth,WETSTOR, DENSITY)
           real :: DRYSTOR,Sdepth,WETSTOR
           real :: DENSITY
           real :: SWE
           SWE = DRYSTOR + WETSTOR
           if (Sdepth > 0.) then
             DENSITY = min(480., SWE/Sdepth)
           else
             DENSITY = 0.
           end if
         end Subroutine SnowDensity

         Subroutine SnowDepthDecrease(DENSITY,Sdepth,SnowMelt, PackMelt)
           real :: DENSITY,Sdepth,SnowMelt
           real :: PackMelt
           if (Sdepth > 0.) then
             PackMelt = max(0.,min( Sdepth/DELT, Sdepth*RHOpack - SnowMelt/DENSITY ))
           else
             PackMelt = 0.
           end if
         end Subroutine SnowDepthDecrease

      Subroutine INFILrunOn(Fdepth,PINFIL, INFIL)
        real :: Fdepth,PINFIL
        real :: INFIL
        if (Fdepth <= poolInfilLimit) then
          INFIL = PINFIL
        else
          INFIL = 0.
        end if
        runOn = PINFIL - INFIL
      end Subroutine INFILrunOn

      Subroutine SurfacePool(Fdepth,Frate,Tsurf,WAPL,WAPS, &
                                            FREEZEPL,poolDrain,poolInfil,THAWPS)
        real :: Fdepth,Frate,Tsurf,WAPL,WAPS
        real :: FREEZEPL,poolDrain,poolInfil,THAWPS
        real :: eta,PIrate,poolVolRemain,poolWavail
        poolVolRemain = max(0., WpoolMax - WAPL - WAPS)
        poolInfil     = min(runOn,poolVolRemain)
        poolRUNOFF    = runOn - poolInfil
        poolWavail    = poolInfil + WAPL/DELT
        if (poolWavail == 0.) then
          poolDrain = 0.
        else if (Fdepth <= poolInfilLimit) then
          poolDrain = poolWavail
        else
          poolDrain = max(0.,min( -Frate*1000., poolWavail ))
        end if
        if ((Tsurf>0.).and.(WAPL==0).and.(WAPS==0.)) then
          PIrate    = 0.
        else
          eta       = LAMBDAice / ( RHOwater * LatentHeat )                                         ! [m2 C-1 day-1]
          PIrate    = (sqrt( max(0.,(0.001*WAPS)**2 - 2.*eta*Tsurf*DELT)))/DELT - (0.001*WAPS)/DELT ! [m day-1]
        end if
        if (PIrate < 0.) then
          FREEZEPL  = 0.
          THAWPS    = min( WAPS/DELT , -PIrate*1000. )
        else
          FREEZEPL  = max( 0.,min( poolInfil + WAPL/DELT - poolDrain*DELT, PIrate*1000. ))
          THAWPS    = 0.
        end if
      end Subroutine SurfacePool

Subroutine DDAYL(doy)
!=============================================================================
! Calculate day length (d d-1) from Julian day and latitude (LAT, degN)
! Author - Marcel van Oijen (CEH-Edinburgh)
!=============================================================================
  integer :: doy                                                      ! (d)
  real    :: DEC, DECC, RAD, DECLIM, DECCMN
  RAD  = pi / 180.                                                    ! (radians deg-1)
  DEC  = -asin (sin (23.45*RAD)*cos (2.*pi*(doy+10.)/365.))           ! (radians)
!  DECC = max(atan(-1./tan(RAD*LAT)),min( atan( 1./tan(RAD*LAT)),DEC)) ! (radians) (Old version)
  if (LAT == 0.) then
    DECLIM = pi/2.
  else
    DECLIM = abs(atan(1./tan(LAT*rad)))
  end if
  DECC = max(-DECLIM, min(DECLIM, DEC))                                ! Simon corrected for polar regions
  YDAYL = DAYL                                                         ! Simon recorded yesterday DAYL
  DAYL  = 0.5 * ( 1. + 2. * asin(tan(RAD*LAT)*tan(DECC)) / pi )        ! (d d-1)
  DECCMN  = max(-DECLIM, 23.45*RAD)
  DAYLMX = 0.5 * ( 1. + 2. * asin(tan(RAD*LAT)*tan(DECCMN)) / pi )     ! (d d-1) Maximum daylength at this latitude

end Subroutine DDAYL

! Calculate PEVAP and PTRAN = potential evaporation and transpiration rates
#ifdef weathergen
  Subroutine PEVAPINPUT(LAI,BASAL)
    real :: LAI,BASAL ! use BASAL to estimate whole sward
    PEVAP  =     exp(-0.5*LAI/BASAL)  * PET                      ! mm d-1 = Partitioning of PET into PEVAP (http://www.fao.org/docrep/x0490e/x0490e04.htm)
    PTRAN  = (1.-exp(-0.5*LAI/BASAL)) * PET                      ! mm d-1 = Partitioning of PET into PTRAN
    PTRAN  = max( 0., PTRAN-0.5*RNINTC )                   ! mm d-1 = Reduction in PTRAN due to wet leaves?
  end Subroutine PEVAPINPUT
#else
  Subroutine PENMAN(LAI,BASAL)
  !=============================================================================
  ! Calculate potential rates of evaporation and transpiration (mm d-1)
  ! Inputs: LAI (m2 m-2), DTR (MJ GR m-2 d-1), RNINTC (mm d-1)
  ! Inputs not in header: VP (kPa), WN (m s-1)
  ! Outputs: PEVAP & PTRAN (mm d-1)
  ! Author - Marcel van Oijen (CEH-Edinburgh)
  !=============================================================================
    real :: LAI,BASAL ! use BASAL to estimate whole sward
    real :: BBRAD, BOLTZM, DTRJM2, LHVAP, NRADC, NRADS
    real :: PENMD, PENMRC, PENMRS, PSYCH, RLWN, SLOPE, SVP, WDF
    DTRJM2 = DTR * 1.E6                                    ! (J GR m-2 d-1)
    BOLTZM = 5.668E-8                                      ! (J m-2 s-1 K-4)
    LHVAP  = 2.4E6                                         ! (J kg-1)
    PSYCH  = 0.067                                         ! (kPA degC-1))
    BBRAD  = BOLTZM * (DAVTMP+273.)**4 * 86400.            ! (J m-2 d-1)
    SVP    = 0.611 * exp(17.4 * DAVTMP / (DAVTMP + 239.))  ! (kPa)
    SLOPE  = 4158.6 * SVP / (DAVTMP + 239.)**2             ! (kPA degC-1)
    RLWN   = BBRAD * max(0.,0.55*(1.-VP/SVP))              ! (J m-2 d-1)
    NRADS  = DTRJM2 * (1.-0.15) - RLWN                     ! (J m-2 d-1)
    NRADC  = DTRJM2 * (1.-0.25) - RLWN                     ! (J m-2 d-1)
    PENMRS = NRADS * SLOPE/(SLOPE+PSYCH)                   ! (J m-2 d-1)
    PENMRC = NRADC * SLOPE/(SLOPE+PSYCH)                   ! (J m-2 d-1)
    WDF    = 2.63 * (1.0 + 0.54 * WN)                      ! (kg m-2 d-1 kPa-1)
    PENMD  = LHVAP * WDF * (SVP-VP) * PSYCH/(SLOPE+PSYCH)  ! (J m-2 d-1)
    PEVAP  =     exp(-0.5*LAI/BASAL)  * (PENMRS + PENMD) / LHVAP ! (mm d-1)
    PTRAN  = (1.-exp(-0.5*LAI/BASAL)) * (PENMRC + PENMD) / LHVAP ! (mm d-1)
    PTRAN  = max( 0., PTRAN-0.5*RNINTC )                   ! (mm d-1)
  end Subroutine PENMAN
#endif

end module environment






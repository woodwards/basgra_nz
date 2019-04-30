module parameters_site

implicit none

!%%%%% Settings for Saerheim 2000

! Simulation period and time step
  real, parameter       :: DELT   =   1.0 ! Model time step

! Geography
  real                  :: LAT            ! Latitude in degrees north

! Atmospheric conditions
  real, parameter       :: CO2A   = 350   ! CO2 concentration in atmosphere (ppm)

! Soil
  real, parameter       :: DRATE  =  50   ! mm d-1 Maximum soil drainage rate
  real                  :: WCI
  real                  :: FWCAD, FWCWP, FWCFC, FWCWET, WCST, BD
  real                  ::  WCAD,  WCWP,  WCFC,  WCWET

! Soil - WINTER PARAMETERS
  real                  :: FGAS, FO2MX, KTSNOW, KRTOTAER, KSNOW ! Simon renamed gamma as KTSNOW
  real, parameter       :: LAMBDAice      = 1.9354e+005  ! J m-1 K-1 d-1 Thermal conductivity of ice
  real                  :: LAMBDAsoil
  real, parameter       :: LatentHeat     = 335000.      ! J kg-1 Latent heat of water fusion
  real, parameter       :: poolInfilLimit =      0.2     ! m Soil frost depth limit for water infiltration
  real                  :: RHOnewSnow, RHOpack
  real, parameter       :: RHOwater       =   1000.      ! kg m-3	Density of water
  real                  :: SWret, SWrf, TmeltFreeze, TrainSnow
  real                  :: WpoolMax

! Soil initial values
  real, parameter       :: DRYSTORI = 0.
  real, parameter       :: FdepthI  = 0.
  real, parameter       :: SDEPTHI  = 0.
  real, parameter       :: TANAERI  = 0.
  real, parameter       :: WAPLI    = 0.
  real, parameter       :: WAPSI    = 0.
  real, parameter       :: WASI     = 0.
  real, parameter       :: WETSTORI = 0.

! Management: harvest dates and irrigation
  integer, dimension(3) :: doyHA
  real, parameter       :: IRRIGF = 0.                   ! Relative irrigation rate (not currently used)

! Mathematical constants
  real, parameter       :: pi   = 3.141592653589793
  real, parameter       :: Freq = 2.*pi / 365.
  real, parameter       :: Kmin = 4.           ! mm C-1 d-1 in SnowMeltWmaxStore()
  real, parameter       :: Ampl = 0.625        ! mm C-1 d-1 Intra-annual amplitude snow melt at 1 degree > 'TmeltFreeze' in SnowMeltWmaxStore()
  real, parameter       :: Bias = Kmin + Ampl  ! mm C-1 d-1 Average snow melting rate at 1 degree above 'TmeltFreeze' in SnowMeltWmaxStore()

end module parameters_site


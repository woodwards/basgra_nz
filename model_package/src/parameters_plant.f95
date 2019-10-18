module parameters_plant

implicit none

! Initial constants
  real :: LOG10CLVI, LOG10CRESI, LOG10CRTI, CSTI, LOG10LAII
  real ::      CLVI,      CRESI,      CRTI,            LAII
  real :: PHENI, TILTOTI, FRTILGI, FRTILGG1I, VERNDI

! Initial constants, continued
  real, parameter :: CLVDI  = 0.
  real, parameter :: YIELDI = 0.
  real, parameter :: CSTUBI = 0.
  real            :: LT50I

! Process parameters
  real :: CLAIV   , COCRESMX, CSTAVM, DAYLB   , DAYLG1G2, DAYLP  , DLMXGE, FSLAMIN
  real :: FSMAX   , HAGERE  , KLAI  , KLUETILG, LAICR   , LAIEFT , LAITIL, LFWIDG ! Simon renamed K to KLAI
  real :: LFWIDV  , NELLVM  , PHENCR, PHY     , RDRSCO  , RDRSMX , RDRTEM, RGENMX
  real :: RGRTG1G2, ROOTDM  , RRDMAX, RUBISC  , LSHAPE  , SIMAX1T, SLAMAX, SLAMIN ! Simon renamed SHAPE as LSHAPE
  real :: TBASE   , TCRES   , TOPTGE, TRANCO  , YG
  real :: RDRTMIN , TVERN   , TVERND, TVERNDMN, AGEH    , RDRROOT, DAYLA, DAYLRV
  real :: FCOCRESMN, KCRT   , RDRTILMIN, RDRHARVMAX, FGRESSI, EBIOMAX, HARVFRD, KBASAL, RDRWMAX, DAYLGEMN, BASALI, ABASAL
  real :: LERVA   , LERVB   , LERGA , LERGB, TRANRFCR, RDRSTUB, DELE, DELD

! Process parameters, continued
  real            :: Dparam, Hparam, KRDRANAER, KRESPHARD, KRSR3H
  real            :: LDT50A, LDT50B, LT50MN, LT50MX, RATEDMX
!  real, parameter :: RDRROOT      =  0.005 ! relative death rate of root CRT. Root currently doesn't do anything in the model. FIXME
!  real, parameter :: RDRSTUB      =  0.2   ! relative death rate of stubble CSTUB. Stubble currently doesn't do anything in the model. FIXME
  real            :: reHardRedDay
  real, parameter :: reHardRedEnd = 91.    ! end date of rehardening reduction period for Northern hemisphere. Adjusted for hemisphere in HardeningSink().
  real            :: THARDMX, TsurfDiff

end module parameters_plant

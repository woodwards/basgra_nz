module plant

use parameters_site
use parameters_plant
use environment

implicit none

! Plant variables
integer :: NOHARV
real :: CRESMX,DAYLGE,FRACTV,GLVSI,GSTSI,LERG,LERV,LUEMXQ,NELLVG,PHENRF,PHOT,RESMOB
real :: RDLVD, ALLOTOT,GRESSI,GSHSI,GLAISI,SOURCE,SINK1T,CSTAV,TGE
real :: RDRFROST,RDRT,RDRTOX,RESPGRT,RESPGSH,RESPHARD,RESPHARDSI,RESNOR,RLEAF,RplantAer,SLANEW
real :: RATEH,reHardPeriod,TV2TIL
real :: CRESMN
real :: ALLOSH, ALLORT, ALLOLV, ALLOST

contains

! Calculate Harvest GSTUB,HARVLA,HARVLV,HARVPH,HARVRE,HARVST,HARVTILG2,HARVFR
Subroutine Harvest(CLV,CRES,CST,year,doy,DAYS_HARVEST,LAI,PHEN,TILG2,TILG1,TILV, &
                             GSTUB,HARVLA,HARVLV,HARVPH,HARVRE,HARVST,HARVTILG2,HARVFR,HARV)
  integer :: doy,year
  integer, dimension(100,3) :: DAYS_HARVEST     ! Simon added third column (percent leaf removed)
  real    :: CLV, CRES, CST, LAI, PHEN, TILG2, TILG1, TILV
  real    :: GSTUB, HARVLV, HARVLA, HARVRE, HARVTILG2, HARVST, HARVPH
  real    :: CLAI, HARVFR, TV1
  integer :: HARV
  integer :: i

  HARV   = 0
  NOHARV = 1
  HARVFR = 0.0
  do i=1,100
    if ( (year==DAYS_HARVEST(i,1)) .and. (doy==DAYS_HARVEST(i,2)) ) then
      HARV   = 1
      NOHARV = 0
      HARVFR = DAYS_HARVEST(i,3) / 100.0  ! Simon read in fraction havested (however this is defined)
	end if
  end do

! CLAIV	m2 leaf m-2 Maximum LAI remaining after harvest, when no tillers elongate
! CLAI	m2 leaf m-2	Maximum LAI remaining after harvest
! (1-FRACTV) * CLAIV = LAI on elongating tillers, assumed to all be harvested
  FRACTV = (TILV + TILG1)/(TILG2 + TILG1 + TILV) ! = Fraction non-elongating tillers (Simon included TILG1)
!  CLAI   = FRACTV * CLAIV                       ! = Maximum LAI remaining after harvest
!  if (LAI <= CLAI) then
!    HARVFR = 0.0             ! Forgot to harvest (1-FRACTV)*LAI?
!  else
!    HARVFR = 1.0 - CLAI/LAI  ! Fraction of leaf that is harvested
!  end if
  HARVFR = max(0.0, 1.0-CLAIV/LAI ) * FRACTV + 1.0 * (1.0 - FRACTV)    ! Simon included harvest of TILG2 leaf in harvest logic
  TV1    = max(0.0, 1.0-CLAIV/LAI ) * FRACTV + HAGERE * (1.0 - FRACTV) ! Simon calculated proportion of CRES harvested
  HARVFR = HARVFR * HARV                                           ! Simon only return HARVFR on HARV days
! HARVFR = Fraction of leaf and reservse in non-elongating tillers that is harvested
! 1.0    = Fraction of leaf              in     elongating tillers that is harvested (I think, or is it HAGERE?)
! HAGERE = Fraction of stem and reserves in     elongating tillers that is harvested

  HARVLA    = (HARV   * LAI * HARVFR) / DELT
  HARVLV    = (HARV   * CLV * HARVFR) / DELT
  HARVPH    = (HARV   * PHEN        ) / DELT           ! PHEN zeroed after each harvest
!  TV1       = (HARVFR * FRACTV) + (HAGERE * (1-FRACTV))! Proportion of CRES harvested
  HARVRE    = (HARV   * CRES * TV1  ) / DELT
  HARVST    = (HARV   * CST * HAGERE) / DELT           ! CST zeroed after each harvest. Simon separated out GSTUB from HARVST
  GSTUB     = (HARV   * CST * (1-HAGERE) ) / DELT      ! Non harvested portion of CST becomes CSTUB
!  HARVTILV   = 0.                                     ! FIXME add effect of grazing on tiller death
!  HARVTILG1   = 0.                                    ! FIXME add effect of grazing on tiller death
  HARVTILG2 = (HARV   * TILG2       ) / DELT           ! TILG2 zeroed after each harvest
end Subroutine Harvest

! Calculate RESNOR (relative amount of CRES)
Subroutine Biomass(CLV,CRES,CST,CSTUB)
  real :: CLV, CRES, CST, CSTUB
  CRESMX = COCRESMX*(CLV + CRES + CST)      ! Maximum reserves in aboveground biomass (not stubble) in terms of C
  CRESMN = COCRESMN*(CLV + CRES + CST)      ! Minimum reserves in aboveground biomass (not stubble) in terms of C
  if (CRESMX>CRESMN) then
    RESNOR = max(0.0, min(1.0, (CRES-CRESMN)/(CRESMX-CRESMN) )) ! Simon revised normalisation of CRES
  else if (CRES.ge.CRESMX) then
	RESNOR = 1.0
  else
	RESNOR = 0.0
  end if
!  RESNOR = max(0.0, min(1.0, CRES/CRESMX )) ! CRES normalised as a proportion of maximum
end Subroutine Biomass

! Calculate phenological changes
Subroutine Phenology(DAYL,PHEN,AGE, DPHEN,GPHEN,HARVPH,FAGE)
  real :: DAYL,PHEN, AGE
  real :: DPHEN,GPHEN,HARVPH,FAGE
  FAGE  = min(1.0, exp( -KAGE * (AGE - AGEH) ))                               ! Simon added effect of sward age factor
  GPHEN = max(0., (DAVTMP-0.01)*0.000144*24. * (min(DAYLP,DAYL)-0.24) )       ! Basically degree days * day length
  DPHEN = 0.
  if (DAYL < DAYLB) then                     ! Simon adjusted resetting of PHEN whenever DAYL < DAYLB
    GPHEN  = 0.0
    DPHEN  = PHEN / DELT
    HARVPH = 0.0
  end if
  PHENRF = max(0.0, min(1.0, (1 - PHEN)/(1 - PHENCR) ))        ! Effect of phenological stage on leaf elongation and appearance on elongating tillers
  DAYLGE = max(0.0, min(1.0, (DAYL - DAYLB)/(DLMXGE - DAYLB) ))! Day length effect on allocation, tillering, leaf appearance, leaf elongation (very crude)
end Subroutine Phenology

! Simon added vernalisation function
! Calculate vernalisation VERN, which allows RGRTVG1 = relative growth rate of generative tillers
Subroutine Vernalisation(DAYL,PHEN,YDAYL,TMMN,TMMX,VERN,VERND, DVERND)
  real :: DAYL, PHEN, YDAYL, TMMN, TMMX
  integer :: VERN
  real :: VERND, DVERND
  real :: X, Y
  ! assume no change in vernalisation
  DVERND = 0.
  ! accumulate sum of cold temperatures
  if (VERN==0) then
    if (TVERN.le.TMMN) then
      DVERND = 0.0
    else if (TVERN.ge.TMMX) then
      DVERND = 1.0
	else
      DVERND = 1.0 ! FIXME. Also FIXME effect of vernalisation on heading date
!	  Y      = (TVERN-TMMN)/(TMMX-TMMN) * 2.0 - 1.0   ! TVERN relative to max and min
!	  X      = acos(Y)                                ! position on first half of cos wave
!      DVERND = 1.0 - X / pi                           ! proportion of day below TVERN
	end if
  end if
  ! does vernalisation occur?
  if ((VERN==0).and.(VERND .ge. TVERND)) then
	VERN = 1
 	VERND = 0.0
    DVERND = 0.0
  end if
  ! reset vernalisation
  if ((VERN==1).and.(DAYL<YDAYL).and.(DAYL<DAYLRV).and.(DAYL>DAYLRV-0.01)) then ! Reset vernalisation when daylength shortens after Solstice
	VERN = 0
	VERND = 0.0
    DVERND = 0.0
  end if
end Subroutine Vernalisation

! Simon renamed Foliage1() to CalcSLA()
! Calculate leaf elongation rates LERV, LERG and SLANEW of new leaves
Subroutine CalcSLA
  real :: EFFTMP, SLAMIN
  EFFTMP = max(TBASE, DAVTMP)
  ! Linear relationship based on Peacock 1976 (who did not include daylength effect)
  LERV   =          max(0., (-0.76 + 0.52*EFFTMP)/1000. ) ! m d-1 leaf elongation rate on vegetative tillers
  LERG   = DAYLGE * max(0., (-5.46 + 2.80*EFFTMP)/1000. ) ! Simon thinks this implies that DAYLGE should have a max of 1.0, so DLMXGE < maximum DAYL
  SLAMIN = SLAMAX * FSLAMIN
  SLANEW = SLAMAX - RESNOR*(SLAMAX-SLAMIN)                ! m2 leaf gC-1 SLA of new leaves (depends on CRES) Simon note unusual units!
end Subroutine CalcSLA

! Calculate light use efficiency LUEMXQ
Subroutine LUECO2TM(PARAV) ! also uses KLUETILG, FRACTV, KLAI
!=============================================================================
! Calculate LUEMXQ (mol CO2 mol-1 PAR quanta)
! Inputs : PARAV (micromol PAR quanta m-2 s-1)
! See equations in M. van Oijen et al. / Ecological Modelling 179 (2004) 39-60
!=============================================================================
  real :: PARAV
  real :: CO2I, EA, EAKMC, EAKMO, EAVCMX, EFF, GAMMAX, KC25, KMC, KMC25
  real :: KMO, KMO25, KOKC, O2, PMAX, R, RUBISCN, T, TMPFAC, VCMAX
  T      = DAVTMP                                            ! degC
  RUBISCN = RUBISC * (1.E6/550000.)                          ! mumol m-2 leaf Rubisco content of upper leaves
  EAVCMX =  68000                                            ! J mol-1 Activation energy for VCMAX
  EAKMC  =  65800                                            ! J mol-1 Activation energy for KMC
  EAKMO  =   1400                                            ! J mol-1 Activation energy for KMO
  KC25   =     20                                            ! mol CO2 mol-1 Rubisco s-1 Catalytic efficiency of Rubisco at 25 degC
  KMC25  =    460                                            ! ppm CO2 Km-value Rubisco for carboxylation at 25 degC
  KMO25  =     33                                            ! % O2	Km-value Rubisco for oxygenation at 25 degC
  KOKC   =      0.21                                         ! Catalytic efficiency ratio Rubisco oxygenation/carboxylation
  O2     =     21                                            ! % O2	Oxygen concentration in chloroplasts
  R      =      8.314                                        ! J K-1 mol-1 Universal gas constant
  CO2I   = 0.7 * CO2A                                        ! ppm CO2 concentration in chloroplasts (Eqn 8)
  VCMAX  = RUBISCN * KC25 * exp((1/298.-1/(T+273))*EAVCMX/R) ! micromol CO2 m-2 leaf s-1 Maximum carboxylation rate in upper leaves (Eqn 7a)
  KMC    =         KMC25 * exp((1/298.-1/(T+273))*EAKMC /R)  ! ppm CO2 Km-value Rubisco for carboxylation (Eqn 7b)
  KMO    =         KMO25 * exp((1/298.-1/(T+273))*EAKMO /R)  ! % O2	Km-value Rubisco for oxygenation (Eqn 7c)
  GAMMAX = 0.5 * KOKC * KMC * O2 / KMO                       ! ppm CO2)	CO2 compensation point at no mitochondrial respiration (Eqn 7d)
  PMAX   = VCMAX * (CO2I-GAMMAX) / (CO2I + KMC * (1+O2/KMO)) ! micromol CO2 m-2 s-1	Photosynthesis rate of upper leaves at light saturation (Eqn 6a)
  TMPFAC = max( 0., min( 1., (T+4.)/5. ) )                   ! Linear decrease of photosynthetic quantum yield at low temperature
  EFF    = TMPFAC * (1/2.1) * (CO2I-GAMMAX) / (4.5*CO2I+10.5*GAMMAX) ! mol CO2 mol-1 PAR quanta	Quantum yield of photosynthesis (Eqn 6b)
  LUEMXQ = EFF*PMAX*(1+KLUETILG*(1-FRACTV)) / (EFF*KLAI*PARAV + PMAX)   ! mol CO2 mol-1 PAR Light-use efficiency (Eqn 5)
end Subroutine LUECO2TM

! Calculate RESPHARDSI respiration for use in Growth()
Subroutine HardeningSink(CLV,DAYL,doy,LT50,Tsurf)
  integer :: doy
  real :: CLV,DAYL,LT50,Tsurf
  real :: doySinceStart, reHardRedStart
  if ( LAT > 0 ) then ! correct for hemisphere
    reHardRedStart = modulo( reHardRedEnd - reHardRedDay, 365. ) ! Rehardening reduction start
  else
    reHardRedStart = modulo( reHardRedEnd + 183 - reHardRedDay, 365. ) ! Rehardening reduction adjusted for hemisphere
  end if
  doySinceStart  = modulo( doy-reHardRedStart       , 365. )
  if ( doySinceStart < (reHardRedDay+0.5*(365.-reHardRedDay)) ) then
    reHardPeriod = max( 0., 1.-doySinceStart/reHardRedDay )
  else
    reHardPeriod = 1.
  end if
  if ( (Tsurf>THARDMX) .or. (LT50<LT50MN) ) then
    RATEH = 0.
  else
    RATEH = reHardPeriod * Hparam * (THARDMX-Tsurf) * (LT50-LT50MN)
  end if
  RESPHARDSI = RATEH * CLV * KRESPHARD * max(0.,min(1., RESNOR*5. )) ! gC m-2 d-1 Sink strength from carbohydrate demand of hardening
end Subroutine HardeningSink

! Calculate all the growth rates
Subroutine Growth(CLV,CRES,CST,PARINT,TILG2,TILG1,TILV,TRANRF, GLV,GRES,GRT,GST)
  real :: CLV,CRES,CST,PARINT,TILG2,TILG1,TILV,TRANRF
  real :: GLV,GRES,GRT,GST
  PHOT     = PARINT * TRANRF * 12. * LUEMXQ * NOHARV               ! gC m-2 d-1 Photosynthesis (12. = gC mol-1) FIXME remove NOHARV
  RESMOB   = (CRES * NOHARV / TCRES) * max(0.,min( 1.,DAVTMP/5. )) ! gC m-2 d-1	Mobilisation of reserves FIXME remove NOHARV
  SOURCE   = RESMOB + PHOT                                         ! gC m-2 d-1	Source strength from photsynthesis and reserve mobilisation
  RESPHARD = min(SOURCE,RESPHARDSI)                                ! gC m-2 d-1	Plant hardening respiration
  ALLOTOT  = SOURCE - RESPHARD                                     ! gC m-2 d-1	Allocation of carbohydrates to sinks other than hardening
  GRESSI   = 0.5 * (RESMOB + max(0.,(CRESMX-CRES)/DELT))           ! gC m-2 d-1	Sink strength of reserve pool
  if (TILG2 > 0.0) then
    CSTAV  = CST/TILG2                                             ! gC tiller-1 Average size of elongating tillers
  else
    CSTAV  = 0.
  end if
  SINK1T   = max(0., 1 - (CSTAV/CSTAVM)) * SIMAX1T                 ! gC tiller-1 d-1 Sink strength of individual elongating tillers
  NELLVG   = PHENRF * NELLVM
  GLAISI   = ((LERV*(TILV+TILG1)*NELLVM*LFWIDV) + (LERG*TILG2*NELLVG*LFWIDG)) * LSHAPE * TRANRF ! m2 leaf m-2 d-1 Potential growth rate of leaf area (Simon added TILG1)
!  GLAISI   = ((LERV*TILV*NELLVM*LFWIDV) + (LERG*TILG2*NELLVG*LFWIDG)) * LSHAPE * TRANRF ! m2 leaf m-2 d-1 Potential growth rate of leaf area
  GLVSI    = max(0.0, (GLAISI * NOHARV / SLANEW) / YG)              ! gC m-2 d-1 Potential growth rate of leaf mass FIXME remove NOHARV
  GSTSI    = max(0.0, (SINK1T * TILG2 * TRANRF * NOHARV) / YG)      ! gC m-2 d-1 Potential growth rate of stems FIXME remove NOHARV
  call Allocation(GRES,GRT,GLV,GST)
end Subroutine Growth

   ! Calculate allocation of CRES to GRES,GRT,GLV,GST
   Subroutine Allocation(GRES,GRT,GLV,GST)
     real :: GRES, GRT, GLV, GST
     GSHSI = GLVSI + GSTSI
!     if (DAYLGE >= 0.1) then   ! Simon thinks maybe this value should be a parameter
     if (DAYLGE >= DAYLGEA) then   ! Simon thinks maybe this value should be a parameter
     ! Situation 1: Growth has priority over storage (spring and growth period)
       ! Calculate amount of assimilates allocated to shoot
       ALLOSH = min( ALLOTOT, GSHSI )
       ! Calculate amount of assimilates allocated to reserves
       GRES   = min( ALLOTOT - ALLOSH, GRESSI)
     else
     ! Situation 2: Storage has priority over shoot (autumn)
       ! Calculate amount of assimilates allocated to reserves
       GRES   = min( ALLOTOT, GRESSI )
       ! Calculate amount of assimilates allocated to shoot
       ALLOSH = min( ALLOTOT - GRES, GSHSI )
     end if
     ! All surplus carbohydrate goes to roots
     ALLORT  = ALLOTOT - ALLOSH - GRES
     if (GSHSI == 0.) GSHSI = 1
     ALLOLV  = GLVSI * (ALLOSH / GSHSI)
     ALLOST  = GSTSI * (ALLOSH / GSHSI)
     GLV     = ALLOLV * YG
     GST     = ALLOST * YG
     GRT     = ALLORT * YG
     RESPGSH = (ALLOLV + ALLOST) * (1-YG) ! gC m-2 d-1 Respiration associated with shoot growth
     RESPGRT =  ALLORT           * (1-YG) ! gC m-2 d-1 Respiration associated with root growth
   end Subroutine Allocation

! Calculate RplantAer = gC m-2 d-1 Aerobic plant respiration
Subroutine PlantRespiration(FO2,RESPHARD)
  real :: FO2,RESPHARD
  real :: fAer
  fAer      = max(0.,min(1., FO2/FO2MX ))
  RplantAer = fAer * ( RESPGRT + RESPGSH + RESPHARD )
end Subroutine PlantRespiration

! Calculate death rates
Subroutine Senescence(CLV,CRT,CSTUB,doy,LAI,LT50,PERMgas,TANAER,TILV,Tsurf, &
                                 DeHardRate,DLAI,DLV,DRT,DSTUB,dTANAER,DTILV,HardRate)
  integer :: doy
  real :: CLV,CRT,CSTUB,DAYL,LAI,LT50,PERMgas,TANAER,TILV,Tsurf
  real :: DeHardRate,DLAI,DLV,DRT,DSTUB,dTANAER,DTILV,HardRate
  real :: RDRS, TV1, TV2
  call AnaerobicDamage(LT50,PERMgas,TANAER, dTANAER)
  call Hardening(CLV,LT50,Tsurf, DeHardRate,HardRate)
  if (LAI < LAICR) then
    TV1 = 0.0
  else
    TV1 = RDRSCO*(LAI-LAICR)/LAICR
  end if
  RDRS   = min(TV1, RDRSMX)
  RDRT   = max(RDRTMIN, RDRTEM * Tsurf)
  TV2    = NOHARV * max(RDRS,RDRT,RDRFROST,RDRTOX) ! d-1 Relative leaf death rate
  TV2TIL = NOHARV * max(RDRS,     RDRFROST,RDRTOX) ! d-1 Relative death rate of non-elongating tillers
  DLAI   = LAI    * TV2
  DLV    = CLV    * TV2
  DSTUB  = CSTUB  * RDRSTUB
  DTILV  = TILV   * TV2TIL
  DRT    = CRT    * RDRROOT

end Subroutine Senescence

   ! Calculate RDRTOX = d-1	Relative death rate of tillers due to anaerobic conditions
   Subroutine AnaerobicDamage(LT50,PERMgas,TANAER, dTANAER)
     real :: LT50,PERMgas,TANAER
     real :: dTANAER,LD50
     if (PERMgas==0.) then      ! d-1 Permeability of soil surface to gas exchange
       dTANAER = 1.             ! d d-1	Change in days since start anaerobic conditions
     else
       dTANAER = -TANAER / DELT ! d d-1	Change in days since start anaerobic conditions
     end if
     LD50 = LDT50A + LDT50B * LT50 ! d Duration of anaerobic conditions at which death rate is half the maximum
     if (TANAER > 0.) then      ! d	Time since start anaerobic conditions
       RDRTOX = KRDRANAER / (1.+exp(-KRDRANAER*(TANAER-LD50))) ! d-1 Relative death rate of tillers due to anaerobic conditions
     else
       RDRTOX = 0.
     end if
     end Subroutine AnaerobicDamage

   ! Calculate RDRFROST, DeHardRate, HardRate
   Subroutine Hardening(CLV,LT50,Tsurf, DeHardRate,HardRate)
     real :: CLV,LT50,Tsurf
     real :: DeHardRate,HardRate
     real :: RATED,RSR3H,RSRDAY
     RSR3H      = 1. / (1.+exp(-KRSR3H*(Tsurf-LT50))) ! d-1	Relative frost survival rate
     ! RDRFROST should be less than 1 to avoid numerical problems
     ! (loss of all biomass but keeping positive reserves). We cap it at 0.5.
     RSRDAY     = RSR3H ! In previous versions we had RSRDAY = RSR3H^8 which understimated survival
     RDRFROST   = min( 0.5, 1. - RSRDAY )             ! d-1 Relative death rate due to frost
     RATED      = min( Dparam*(LT50MX-LT50)*(Tsurf+TsurfDiff), (LT50MX-LT50)/DELT ) ! °C d-1 Potential rate of dehardening, if below limit set by RATEDMX
     DeHardRate = max(0.,min( RATEDMX, RATED ))
     if ( CLV > 0.0 ) then
       HardRate   = RESPHARD / (CLV * KRESPHARD)
     else
       HardRate   = 0.0
     end if
   end Subroutine Hardening

! Simon added decomposition function
! Calculate decompositon of dead leaf
Subroutine Decomposition(CLVD,DAVTMP,WCL, DLVD,RDLVD)
  real :: CLVD,DAVTMP,WCL
  real :: DLVD
  real :: PSIA,PSIB,SWCS,BD,PSIS,DELD,DELE
  real :: EBIOMASSMAX,EBIOMASS,CT,CP,WORMS
  real :: DTEMP,DWATER,DECOMP,RDLVD
  EBIOMASSMAX = 131.0
  PSIA    = 3.0e-3                 ! Te Kowhai silt loam FIXME link to WC params
  PSIB    = 7.75                   ! Te Kowhai silt loam FIXME link to WC params
  SWCS    = WCL                    ! Volumetric soil water content near surface
  BD      = 1.1                    ! Bulk density (Singleton pers comm) FIXME link to soil params
  PSIS    =  -PSIA * (SWCS ** PSIB) ! Soil water tension near surface
  DELD    = 0.0148
  DELE    = 0.0005
  ! Calculate number of worms and their grazing of dead matter
  ! Numbers at surface based on Baker et al., driven by GWCS
  ! Activity based on Daniels
  EBIOMASS= max(0.0, min(1.0, 5.0*SWCS/BD-1.0)) * EBIOMASSMAX
  if (DAVTMP > 20.0) then
    CT    = 0.0
  else
    CT    = 0.515 * (20.0-DAVTMP) ** 1.84 * exp(-0.297*(20.0-DAVTMP))/2.345 ! Daniels
  end if
  if (PSIS > -12.3) then
    CP    = 1.0
  else
    CP    = 0.549 * (-PSIS) ** 0.793 * exp(0.113 * PSIS) ! Daniels
  end if
  WORMS   = DELE * EBIOMASS * CT * CP
  ! Calculate decomposition, based on Andren paper
  if (DAVTMP > 0.0) then
    DTEMP = 2.0 ** ((DAVTMP - 20.0)/10.0)
  else
    DTEMP = 0.0
  end if
  DWATER  = max(0.0, min(1.0, log(-7580.0 / PSIS) / log(-7580.0 / (-10.0))))
  if (RAIN > 0.0) DWATER = 1.0       ! decomp on rain days even if dry soil, McCall 1984
  DECOMP  = DELD * DTEMP * DWATER  ! total relative decomposition rate
  ! Total relative dead matter disappearence rate
  RDLVD   = DECOMP + WORMS
  DLVD    = CLVD    * RDLVD
end Subroutine Decomposition

! Simon renamed Foliage2() to Tillering()
! Calculate GLAI,GTILV,TILVG1,TILG1G2
Subroutine Tillering(DAYL,GLV,LAI,TILV,TILG1,TRANRF,Tsurf,VERN,FAGE, GLAI,GTILV,TILVG1,TILG1G2)
  real    :: DAYL,GLV,LAI,TILV,TILG1,TRANRF,Tsurf,FAGE
  integer :: VERN
  real    :: GLAI,GTILV,TILVG1,TILG1G2
  real    :: RGRTV,RGRTVG1,TV1,TV2
  GLAI    = SLANEW * GLV                                                      ! Simon note SLANEW is in m2 leaf gC-1
  if (Tsurf < TBASE) then
    TV1   = 0.
  else
    TV1   = Tsurf/PHY                                                         ! d-1 Potential leaf appearence rate
  end if
  RLEAF   = TV1 * TRANRF * DAYLGE * ( FRACTV + PHENRF * (1-FRACTV) )          ! d-1 Leaf appearance rate (Simon moved NOHARV switch)
!  RLEAF   = TV1 * NOHARV * TRANRF * DAYLGE * ( FRACTV + PHENRF * (1-FRACTV) ) ! d-1 Leaf appearance rate
  TV2     = max( 0.0, min(FSMAX, LAITIL - LAIEFT*LAI ))                       ! d-1 Ratio of tiller and leaf apearance
  RGRTV   = max( 0.0       , TV2 * RESNOR * RLEAF )                           ! d-1 Relative rate of vegetative tiller appearence
  GTILV   = TILV  * RGRTV * NOHARV                                            ! Simon moved NOHARV switch to here
  TGE     = max( 0.0       , 1.0 - (abs(DAVTMP - TOPTGE))/(TOPTGE-TBASE))     ! Temperature effect on initiation of elongation in tillers
  RGRTVG1 = NOHARV * DAYLGE * TGE * RGENMX * VERN                             ! d-1 Relative rate of vegetative tiller conversion to generative
  RGRTVG1 = min( 1.0 - TV2TIL, RGRTVG1)                                       ! Make sure TILV doesn't all disappear in one day
  TILVG1  = TILV  * RGRTVG1
  if (DAYL > DAYLG1G2) then                                                   ! Generative tiller conversion controlled by DAYL
    TILG1G2 = TILG1 * RGRTG1G2                                                ! d-1 Rate of generative tiller conversion to elongating
  else if (YDAYL < DAYL) then
    TILG1G2 = 0.                                                              ! no conversion yet
  else
    TILG1G2 = TILG1                                                           ! Simon convert all remaining tillers once DAYL < DAYLG1G2
  end if
end Subroutine Tillering

end module plant

# test the WCST code from BASGRA (turns out FORTRAN almost runs in R!)

BD = 0.9
TRANCO = 6

SP = 30:90
WCST = SP/100*BD
WCFC  = (2.62 + 0.595 * SP)/100*BD
WCWP  = (-7.92 + 0.593 * SP)/100*BD
WCWET = 0.95 * WCST                 # Simon rough estimate
WCAD  = 0.3 * WCWP                  # Simon rough estimate
WCCR = WCWP + (WCFC - WCWP) * max(0.0, PTRAN/(PTRAN+TRANCO))  # = m3 m-3 Critical water content below which transpiration is reduced (Eqn 1)

plot(WCST, WCST, col="black", ylim=c(0.0,1.2))
points(WCST, WCWET, col="blue") # ok
points(WCST, WCFC, col="green") # too high
points(WCST, WCCR, col="red") # ok
points(WCST, WCWP, col="lightblue") # ok
points(WCST, WCAD, col="grey") # ok

WCST = 0.32

PEVAP = 1.5
PTRAN = 4
WCLM = 0.25
WCL = 0.25
ROOTDM = 0.55
Fdepth = 0
WAL = 1000. * WCL * (ROOTDM-Fdepth)
DELT = 1

EVAP = PEVAP * max(0., min(1., (WCLM-WCAD)/(WCFC-WCAD) ))      # = mm d-1 Reduction in evaporation due to soil water content
WCCR = WCWP + (WCFC - WCWP) * max(0.0, PTRAN/(PTRAN+TRANCO))  # = m3 m-3 Critical water content below which transpiration is reduced (Eqn 1)
if (WCL > WCCR) {                                          # Transpiraiton reduction factor (Fig 4)
FR = max(0., min(1., (WCST-WCL)/(WCST-WCWET) ))             # Transpiration reduction in wet conditions
}else{ if (WCCR > WCWP) {                                    # Gradual reduction due to dryness
FR = max(0., min(1., (WCL-WCWP)/(WCCR-WCWP)  ))             # Transpiration reduction in dry conditions
}else{
  FR = 0.0                                                     # Simon added this case explicitly instead of setting lower bound to WCCR
}}
TRAN = PTRAN * FR                                             # = mm d-1 Transpiration reduction due to WCL
WAAD = 1000. * WCAD * (ROOTDM-Fdepth)                         # = mm Water in non-frozen soil at air dryness, Simon modified to ROOTDM
if (EVAP+TRAN > 0.) {
AVAILF = min( 1., ((WAL-WAAD)/DELT) / (EVAP+TRAN) )         # = Prevent WAL falling below WAAD
}else{
  AVAILF = 0
}
EVAP = EVAP * AVAILF                                          # (mm d-1)
TRAN = TRAN * AVAILF                                          # (mm d-1)
if (PTRAN > 0.) {
TRANRF = TRAN / PTRAN                                       # (-) Water restriction on plant processes
}else{
  TRANRF = 1                                                  # (-) Assume no restriction when PTRAN==0
}


!======================== include file "fwa.h" =========================

!      parameters for adding anomalous fresh water pulses

!      isfwa1     = starting i index of section 1 for freshwater
!      iefwa1     = ending i index of section 1 for freshwater
!      isfwa2     = starting i index of section 2 for freshwater
!      iefwa2     = ending i index of section 2 for freshwater
!      jsfwa      = starting j index for freshwater
!      jefwa      = ending j index for freshwater
!      mrfwa      = regional mask region for freshwater
!      fwaflxi    = initial fresh water flux (Sv)
!      fwayri     = year to start fresh water flux
!      fwayrf     = year to end fresh water flux
!      fwarate    = rate of increase in flux (Sv year-1)
!      areafwa    = area of fresh water anomaly
!      areafwc    = area of fresh water compensation
!      compensate = flag to compensate for flux everywhere else
!      fwawt      = weighting for fresh water anomaly area
!      fwcwt      = weighting for fresh water compensation area

       integer isfwa1, iefwa1, isfwa2, iefwa2, jsfwa, jefwa, mrfwa
       common /fwa_i/ isfwa1, iefwa1, isfwa2, iefwa2, jsfwa, jefwa
       common /fwa_i/ mrfwa

       logical compensate
       common /fwa_l/ compensate

       real fwaflxi, fwayri, fwayrf, fwarate, areafwa, areafwc
       real fwawt, fwcwt
       common /fwa_r/ fwaflxi, fwayri, fwayrf, fwarate, areafwa, areafwc
       common /fwa_r/ fwawt(imt,jmt), fwcwt(imt,jmt)

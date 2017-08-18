!====================== include file "ctavg.h" =========================

!     common for tracer averages within horizontal regions

!     sumbt = volume weighted tracer sum for a given region
!     avgbt = volume weighted tracer average for a given region
!     sumbk = volume weighted tracer sum for a region at a k level
!     avgbk = volume weighted tracer average for a region at a k level
!     sumgt = total global volume weighted tracer sum
!     avggt = total global volume weighted tracer average
!     sumgk = total global volume weighted tracer sum at a k level
!     avggk = total global volume weighted tracer average at a k level
!     sumbf = area weighted tracer flux sum for a given region
!     avgbf = average area weighted tracer flux for a region
!     sumgf = global area weighted tracer flux sum
!     avggf = global average area weighted tracer flux

      real sumbk, sumbt, sumgk, sumgt, sumbf, sumgf, avgbk, avgbt
      real avggk, avggt, avgbf, avggf

      common /tavg_r/ sumbk(nhreg,km,nt), sumbt(nhreg,nt), sumgk(km,nt)
      common /tavg_r/ sumgt(nt), sumbf(nhreg,nt), sumgf(nt)
      common /tavg_r/ avgbk(nhreg,km,nt), avgbt(nhreg,nt), avggk(km,nt)
      common /tavg_r/ avggt(nt), avgbf(nhreg,nt), avggf(nt)

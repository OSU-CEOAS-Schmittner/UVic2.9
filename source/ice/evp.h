!======================== include file "evp.h" =========================

!     arrays for the evp ice dynamics code
#if defined O_ice_evp

      integer ndte
      parameter (ndte=100)

      real dtei, xyminevp, floor, waterx, watery, strairx, strairy
      real ecc2, ecc2m, ecc2p, fmass, diff1, diff2, umass, eyc
      real zetamin, zetan, zetas, zetae, zetaw, etan, etas, etae, etaw
      real sig11n, sig11e, sig11s, sig11w, sig12n, sig12e
      real sig12s, sig12w, sig22n, sig22e, sig22s, sig22w
      real HTN4, HTE4, dxt8, dyt8, prssn, prsss, prsse, prssw
      real a2na, a2sa, a2ea, a2wa, b2n, b2s, b2e, b2w, h2n
      real h2s, h2e, h2w, edy, edx, eHN, eHE, eHNm, eHEm

      common /dyn1/ dtei, xyminevp, floor
      common /dyn1/ waterx(imt,jmt), watery(imt,jmt)
      common /dyn1/ strairx(imt,jmt), strairy(imt,jmt)
      common /dyn1/ ecc2, ecc2m, ecc2p
      common /dyn1/ fmass(imt,jmt),diff1,diff2

      common /dyn2/ umass(imt,jmt), eyc
      common /dyn2/ zetamin, zetan(imt,jmt), zetas(imt,jmt)
      common /dyn2/ zetae(imt,jmt), zetaw(imt,jmt), etan(imt,jmt)
      common /dyn2/ etas(imt,jmt), etae(imt,jmt), etaw(imt,jmt)
      common /dyn2/ sig11n(imt,jmt), sig11e(imt,jmt), sig11s(imt,jmt)
      common /dyn2/ sig11w(imt,jmt), sig12n(imt,jmt), sig12e(imt,jmt)
      common /dyn2/ sig12s(imt,jmt), sig12w(imt,jmt), sig22n(imt,jmt)
      common /dyn2/ sig22e(imt,jmt), sig22s(imt,jmt), sig22w(imt,jmt)

      common /dyn3/ HTN4(imt,jmt), HTE4(imt,jmt)
      common /dyn3/ dxt8(imt,jmt), dyt8(imt,jmt), prssn(imt,jmt)
      common /dyn3/ prsss(imt,jmt), prsse(imt,jmt), prssw(imt,jmt)
      common /dyn3/ a2na(imt,jmt), a2sa(imt,jmt), a2ea(imt,jmt)
      common /dyn3/ a2wa(imt,jmt), b2n(imt,jmt), b2s(imt,jmt)
      common /dyn3/ b2e(imt,jmt), b2w(imt,jmt), h2n(imt,jmt)
      common /dyn3/ h2s(imt,jmt), h2e(imt,jmt), h2w(imt,jmt)
      common /dyn3/ edy(imt,jmt), edx(imt,jmt), eHN(imt,jmt)
      common /dyn3/ eHE(imt,jmt), eHNm(imt,jmt), eHEm(imt,jmt)

#endif

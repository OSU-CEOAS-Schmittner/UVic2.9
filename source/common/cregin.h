!====================== include file "cregin.h" ========================

!     variables used for computing regional tracer averages (see
!     "reg1st.F" & sub "diagt1" in "tracer.F") and for computing term
!     balances for tracer and momentum equations (see "clinic.F",
!     "tracer.F" and "diag2.F")

!     mskhr = mask field defining regions in the horizontal
!              (eg: mskhr(i,j) = n indicates point (i,j) is in the
!               "nth" horizontal region   where n=1..nhreg)
!               The "mskhr" masks are used in "diagt1" when
!               computing volume weighted tracer averages and in
!               "clinic.F", "tracer.F" and "diag2.F" when computing
!               term balances for tracers and momentum.
!     mskvr  =  mask field defining regions in the vertical
!              (eg: mskvr(k) = m indicates all points under a horizontal
!               mask at level "k" are in the "mth" vertical region
!               where m=1..nvreg)
!               The "mskvr" masks are used in "diag.F", but not
!               in "diagt1", where tracer averages are calculated
!               for each k-level.

!     hregnm = horizontal region name
!     vregnm = vertical region name
!     volbt  = total volume under a given horizontal region
!     volbk  = volume contained in a horizontal region at level "k"
!     volgt  = total ocean volume
!     volgk  = total ocean volume at level "k"
!     areab  = total ocean surface area for a given horizontal region
!     areag  = total ocean surface area

!     volt   = ocean volume within a particular horizontal & vertical
!              region (on the "t" grid) for tracer term balances
!              volt(0) represents the sum of all regions
!     rvolt  = 1/volt ( 0.0 if volt = 0.0)
!     areat  = horizontal ocean surface area corresponding to "volt"
!              areat(0) represents the sum of all regions
!     rareat = 1/areat ( 0.0 if areat = 0.0)
!     volu   = ocean volume within a particular horizontal & vertical
!              region (on the "u" grid) for momentum term balances
!              volu(0) represents the sum of all regions
!     rvolu  = 1/volu ( 0.0 if volu = 0.0)
!     areau  = horizontal ocean area corresponding to "volu"
!              areau(0) represents the sum of all regions
!     rareau = 1/areau ( 0.0 if areau = 0.0)
!     llvreg = level limits for defining vertical regions in term
!              balance calculations (not used in computing volume
!              weighted tracer averages)
!              (eg: llvreg(3,1) = 4... means that starting level for
!                  the third region in the vertical is 4. similarly,
!                  llvreg(3,2) = 6 means the ending level is 6 for that
!                  region. note regions should not overlap.)

      character(40) :: hregnm
      character(20) :: vregnm
      common /cregn_c/ hregnm(nhreg), vregnm(nvreg)

      integer mskhr, mskvr, llvreg
      common /cregn_i/ mskhr(imt,jmt), mskvr(km), llvreg(numreg,2)

      real volbk, volbt, volgk, volgt, areab, areag, volt, volu
      real areat, areau, rvolt, rvolu, rareat, rareau
      common /cregn_r/ volbk(nhreg,km), volbt(nhreg), volgk(km)
      common /cregn_r/ volgt, areab(nhreg), areag, volt(0:numreg)
      common /cregn_r/ volu(0:numreg), areat(0:numreg), areau(0:numreg)
      common /cregn_r/ rvolt(0:numreg), rvolu(0:numreg)
      common /cregn_r/ rareat(0:numreg), rareau(0:numreg)

!======================= include file "ctdbc.h" ========================

!           time dependent S.B.C. data for UVic_ESCM

!     ntdbc  = number of ocean S.B.C.
!               1 is for taux (dynes/cm**2)
!               2 is for tauy (dynes/cm**2)
!               3 is sst (converted to heat flux in "atmos.F")
!               4 is sss (converted to fresh water flux in "atmos.F")

!               for generality in data assimilation situations, each
!               S.B.C. is treated separately.  the default case assumes,
!               however, that all S.B.C. are monthly.

!     ntdrec = number of time dependent records on disk.
!               note: default case assumes 12 monthly S.B.C. records
!                     (defined at the centres of the months) from
!                     the MOM climatological database.

!     tdrec   = times at which the data records are defined.
!               times are assumed to define the centres of the records.
!               (ie: when using  monthly climatological data and the
!                first record is January, "tdrec(1)" is defined in days
!                as day 15.5)

!     aprec   = period lengths for the data records. (number of days
!               per month)

!     isbcstart= index to time at beginning of first S.B.C. data record

!     method   = interpolation method = (0,1,2,3) = (none, linear,
!                equal linear, equal area) see timeinterp for details.

!     dstamp   = time stamps marking the end of each S.B.C. data record

!     period  = (true,false)  if the data records are to be treated as
!               (periodic, not periodic). if periodic, then the model
!               time is always mapped into the data records. if not,
!               then record 1 is used for all model time before the
!               beginning of the dataset and record "ntdrec" is used
!               for all model time after the end of the dataset.

!     daysbc  = model time in days after start of S.B.C.

!     iprevd  = index for pointing to the nearest S.B.C. data on disk
!               which the current model time has passed.
!     inextd  = index for pointing to the nearest S.B.C. data on disk
!               which the current model time has not yet reached.
!     iprevm  = index for pointing to the nearest S.B.C. data in memory
!               which the current model time has passed.
!     inextm  = index for pointing to the nearest S.B.C. data in memory
!               which the current model time has not yet reached.
!     wprev   = interpolation factor for defining data at the current
!               time step (the weight attributed to the data from the
!               "iprevd" time interval)
!               schematically the interpolation is defined at the time
!               step by:
!               data = wprev*data(iprevd) + (1-wprev)*data(inextd)
!     obc     = ocean S.B.C.

!     iotdbc  = units for time dependent monthly boundary conditions
!               based on hellerman stress and levitus temperature and
!               salt from the MOM database.

      integer maxbcr, ntdbc
      parameter (maxbcr=12)
      parameter (ntdbc = 8
#if defined O_restorst || defined O_replacst || defined O_read_my_stf
     &                 + 2
#endif
#if defined O_embm_awind
     &                 + 1
#endif
#if defined O_carbon_co2_2d
     &                 + 2
#endif
#if defined O_shortwave
     &                 + 1
#endif
#if defined O_mtlm
     &                 + 1
#endif
#if defined O_mobi_iron
     &                 + 1
#endif
#if defined O_mobi_silicon
     &                 + 1
#endif
     &                    )

      character(32) :: rstamp1, rstamp2, dstamp
      common /ctdbc_c/ rstamp1, rstamp2, dstamp(maxbcr,ntdbc)

      integer ntdrec, iprevd, inextd
      integer iprevm, inextm, iotdbc
      integer isbcstart, method
      common /ctdbc_i/ ntdrec(ntdbc), iprevd(ntdbc), inextd(ntdbc)
      common /ctdbc_i/ iprevm(ntdbc), inextm(ntdbc), iotdbc(ntdbc)
      common /ctdbc_i/ isbcstart(ntdbc), method

      logical period, rdtdbc
      common /ctdbc_l/ period(ntdbc), rdtdbc(ntdbc)

      real wprev, aprec, tdrec, avgper, daysbc, obc
      common /ctdbc_r/ wprev(ntdbc), aprec(maxbcr,ntdbc)
      common /ctdbc_r/ tdrec(maxbcr,ntdbc), avgper(maxbcr)
      common /ctdbc_r/ daysbc(ntdbc)
#if defined O_sbc_in_memory
!     this is set up only for monthly data.
      common /ctdbc_r/ obc(imt,jmt,ntdbc,12)
#else
      common /ctdbc_r/ obc(imt,jmt,ntdbc,2)
#endif

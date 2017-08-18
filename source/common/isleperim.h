!===================== include file "isleperim.h" ======================

!                      topography changes

!     kmt_opt     = a list of possible modifications to kmt to fix
!                   one flagged problem.
!     kmt_changes = list of collected changes to kmt
!     nchanges    = number of changes to kmt in kmt_changes
!     auto_kmt_changes = tells whether any changes have been made to the
!                        kmt field as a result define options or
!                        interactive user actions.
!     max_opt     = max number of options for kmt_opt
!     len_opt     = max number of kmt points changed per option
!     io_del_kmt  = io unit for writing delta.kmt.###x###x###.h file
!     n_del_kmt   = number of kmt changes in delta.kmt.###x###x###.h

      integer max_opt, len_opt, max_change
      parameter (max_opt=3, len_opt=10, max_change=100)

      integer kmt_opt, kmt_changes, nchanges, io_del_kmt, n_del_kmt
      common /kmtchg_i/ kmt_opt(max_opt, len_opt, 4)
      common /kmtchg_i/ kmt_changes(max_change, 4)
      common /kmtchg_i/ nchanges, io_del_kmt, n_del_kmt

      logical auto_kmt_changes
      common /kmtchg_l/ auto_kmt_changes

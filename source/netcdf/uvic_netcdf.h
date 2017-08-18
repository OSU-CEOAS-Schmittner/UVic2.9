!==================== include file "uvic_netcdf.h" =====================

!     arrays for netcdf interface

!     max_num_files = maximum number of files opened simultaneosly
!     list_ncid     = list of open netcdf file ids
!     list_names    = list of open netdcf file names
!     num           = number of open netcdf files

      integer list_ncid, max_num_files, num
      character(120) :: list_names

      parameter (max_num_files=200)

      common /uvic_netcdf_c/ list_names(max_num_files)
      common /uvic_netcdf_i/ list_ncid(max_num_files), num

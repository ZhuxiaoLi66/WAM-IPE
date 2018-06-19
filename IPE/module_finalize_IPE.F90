!nm20140919: subroutine to finalize IPE
! DATE: 19 September, 2014
!********************************************
!***      Copyright 2014 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Whole Atmosphere Model(WAM)-Ionosphere Plasmasphere Electrodynamics(IPE) model Interface
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!-------------------------------------------- 
!
module module_finalize_IPE
   implicit none	

   private

   public :: finalize_IPE
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine finalize_IPE (model,rc )
  use ESMF
  use nnt_types_module
  use module_decomp
  USE module_close_files,ONLY: close_files
  implicit none
  type(ESMF_GridComp)  :: model
  integer :: rc
!---
  call ESMF_LogWrite("sub-finalize_IPE start:", ESMF_LOGMSG_INFO, rc=rc)

!--- copied from driver_ipe_sms.f90

! DEallocate arrays
      CALL ESMF_LogWrite("sub-finalize_IPE: allocate_arrays1", ESMF_LOGMSG_INFO, rc=rc)
      CALL allocate_arrays ( 1 )

! close all open files
      CALL ESMF_LogWrite("sub-finalize_IPE: close_files", ESMF_LOGMSG_INFO, rc=rc)
      CALL close_files ( )




      CALL ESMF_LogWrite("sub-finalize_IPE: stop", ESMF_LOGMSG_INFO, rc=rc)

      call stop

      call nnt_stop('module_finalize_IPE.F90:L60',0,PPP_EXIT)
!---end copy
  call ESMF_LogWrite("sub-finalize_IPE finished:", ESMF_LOGMSG_INFO, rc=rc)

  rc = ESMF_SUCCESS

  end subroutine finalize_IPE
end module module_finalize_IPE

!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! Clay Blankenship, 22 Feb 2024 based on GLDAS 2 reader


! 
! !MODULE: HydroSCSobsMod
! \label(HydroSCSobsMod)
!
! !INTERFACE:
module HydroSCSobsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015   Sujay Kumar  Initial Specification
!  22 Feb 2024  Clay Blankenship changes for HydroSCS
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: HydroSCSobsinit !Initializes structures for reading HydroSCS data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: HydroSCSobs !Object to hold HydroSCS observation attributes
!EOP

  type, public :: HydroSCSdec
     character*100           :: odir
     character*20            :: model_name
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: yr
     integer                 :: mo
     logical                 :: startFlag
     real                    :: datares
  end type HydroSCSdec
     
  type(HydroSCSdec), allocatable :: HydroSCSObs(:)

contains
  
!BOP
! 
! !ROUTINE: HydroSCSobsInit
! \label{HydroSCSobsInit}
!
! !INTERFACE: 
  subroutine HydroSCSobsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the GIMMS NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    if(.not.allocated(HydroSCSobs)) then 
       allocate(HydroSCSobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, HydroSCSobs(i)%odir, &
         label='HydroSCS data directory:',rc=status)
    call LVT_verify(status, 'HydroSCS data directory: not defined')

    !CBB Not used?
    !call ESMF_ConfigGetAttribute(LVT_Config, HydroSCSobs(i)%model_name, &
    !     label='HydroSCS data model name:',rc=status)
    !call LVT_verify(status, 'HydroSCS data model name: not defined')

    call LVT_update_timestep(LVT_rc, 2592000)

    HydroSCSobs(i)%gridDesc = 0
        
    HydroSCSobs(i)%nc = 7200
    HydroSCSobs(i)%nr = 3000

    !filling the items needed by the interpolation library
    HydroSCSobs(i)%gridDesc(1) = 0                  !CBB 0=lat/lon
    HydroSCSobs(i)%gridDesc(2) = HydroSCSobs(i)%nc
    HydroSCSobs(i)%gridDesc(3) = HydroSCSobs(i)%nr
    HydroSCSobs(i)%gridDesc(4) = -59.975
    HydroSCSobs(i)%gridDesc(5) = -179.975
    HydroSCSobs(i)%gridDesc(6) = 128         !fixed acc. to ref. manual
    HydroSCSobs(i)%gridDesc(7) = 89.975
    HydroSCSobs(i)%gridDesc(8) = 179.975
    HydroSCSobs(i)%gridDesc(9) = 0.05
    HydroSCSobs(i)%gridDesc(10) = 0.05
    HydroSCSobs(i)%gridDesc(11) = 64         !fixed acc. to ref. manual
    HydroSCSobs(i)%gridDesc(20) = 64         !means N-S ordering acc. to LIS72_referencemanual
                                         !copied from GLDAS2 reader (CBB)

    HydroSCSobs(i)%datares  = 0.05

    if(LVT_isAtAfinerResolution(HydroSCSobs(i)%datares)) then
       
       allocate(HydroSCSobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(HydroSCSobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(HydroSCSobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(HydroSCSobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            HydroSCSobs(i)%rlat, &
            HydroSCSobs(i)%rlon, &
            HydroSCSobs(i)%n11)
    else
       allocate(HydroSCSobs(i)%n11(HydroSCSobs(i)%nc*HydroSCSobs(i)%nr))
       call upscaleByAveraging_input(HydroSCSobs(i)%gridDesc,&
            LVT_rc%gridDesc,HydroSCSobs(i)%nc*HydroSCSobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,HydroSCSobs(i)%n11)
    endif

    HydroSCSobs(i)%yr = -1
    HydroSCSobs(i)%mo = LVT_rc%mo
    HydroSCSobs(i)%startFlag = .false. 

  end subroutine HydroSCSobsinit


end module HydroSCSobsMod

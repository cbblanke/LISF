!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Clay Blankenship, 22 Feb 2024 from GLDAS2 reader
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readHydroSCSObs
! \label{readHydroSCSObs}
!
! !INTERFACE: 
subroutine readHydroSCSObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use HydroSCSobsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This plugin processes the HydroSCS downscaled forcing produced by the SERVIR-S2S team 
!   
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
!  12 Feb 2024: Clay Blankenship, update for HydroSCS
!EOP

  integer                :: c,r, k,nc,nr,tindex
  integer                :: nDims, nVars, nAtts, unlimId
  character(len=100)     :: var_name
  character(len=10)      :: var_suffix
  integer                :: flag
  integer                :: ftn
  character*100          :: fname
  logical                :: file_exists
  integer                :: lwdownid, swdownid, rainfid
  integer                :: qairid, windid, tairid, psurfid
  real                   :: tair(HydroSCSobs(source)%nc, HydroSCSobs(source)%nr)
  real                   :: qair(HydroSCSobs(source)%nc, HydroSCSobs(source)%nr)
  real                   :: swdown(HydroSCSobs(source)%nc, HydroSCSobs(source)%nr)
  real                   :: lwdown(HydroSCSobs(source)%nc, HydroSCSobs(source)%nr)
  real                   :: wind(HydroSCSobs(source)%nc, HydroSCSobs(source)%nr)
  real                   :: psurf(HydroSCSobs(source)%nc, HydroSCSobs(source)%nr)
  real                   :: rainf(HydroSCSobs(source)%nc, HydroSCSobs(source)%nr)


  real                   :: lwdown_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: swdown_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: rainf_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: qair_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: wind_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: tair_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: psurf_ip(LVT_rc%lnc,LVT_rc%lnr)

  integer                :: iret

  lwdown_ip   =  LVT_rc%udef
  swdown_ip   =  LVT_rc%udef
  rainf_ip    =  LVT_rc%udef
  qair_ip     =  LVT_rc%udef
  wind_ip     =  LVT_rc%udef
  tair_ip     =  LVT_rc%udef
  psurf_ip    =  LVT_rc%udef

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  nc = HydroSCSobs(source)%nc
  nr = HydroSCSobs(source)%nr

  if((HydroSCSobs(source)%mo.ne.LVT_rc%d_nmo(source)).or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 

     if(HydroSCSobs(source)%startFlag) then 
        HydroSCSobs(source)%startFlag = .false. 
     endif

     HydroSCSobs(source)%yr = LVT_rc%d_nyr(source)
     HydroSCSobs(source)%mo = LVT_rc%d_nmo(source)

     call create_HydroSCS_filename(HydroSCSobs(source)%odir,&
          LVT_rc%dyr(source),&
          LVT_rc%dmo(source),&
          LVT_rc%dda(source),& 
          LVT_rc%dhr(source),&
          fname)
     
     inquire(file=trim(fname),exist=file_exists) 

     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading HydroSCS file ',trim(fname)
        
        iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
             ncid = ftn)
        if(iret.eq.0) then 

           iret = nf90_inquire(ftn, nDims, nVars, nAtts, unlimId)
           iret = nf90_inquire_variable(ftn, nVars, var_name)
           !var_suffix = var_name(index(var_name, 'ave')+3:len(var_name))
           
           !call LVT_verify(nf90_inq_varid(ftn,"LWdown_GDS0_SFC_ave"//trim(var_suffix),lwdownid),&
           !      'nf90_inq_varid failed for LWdown_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Tair",tairid),'nf90_inq_varid failed for Tair')
           call LVT_verify(nf90_inq_varid(ftn,"Qair",qairid),'nf90_inq_varid failed for Qair')
           call LVT_verify(nf90_inq_varid(ftn,"SWdown",swdownid),'nf90_inq_varid failed for SWdown')
           call LVT_verify(nf90_inq_varid(ftn,"LWdown",lwdownid),'nf90_inq_varid failed for LWdown')
           call LVT_verify(nf90_inq_varid(ftn,"Wind",windid),'nf90_inq_varid failed for Wind')
           call LVT_verify(nf90_inq_varid(ftn,"Psurf",psurfid),'nf90_inq_varid failed for Psurf')
           call LVT_verify(nf90_inq_varid(ftn,"Rainf",rainfid),'nf90_inq_varid failed for Rainf')

           call LVT_verify(nf90_get_var(ftn,tairid,tair,&
                start=(/1,1/),&
                count=(/HydroSCSobs(source)%nc,HydroSCSobs(source)%nr/)),&
                'Error in nf90_get_var Tair')
           call LVT_verify(nf90_get_var(ftn,qairid,qair,&
                start=(/1,1/),&
                count=(/HydroSCSobs(source)%nc,HydroSCSobs(source)%nr/)),&
                'Error in nf90_get_var Qair')
           call LVT_verify(nf90_get_var(ftn,swdownid,swdown,&
                start=(/1,1/),&
                count=(/HydroSCSobs(source)%nc,HydroSCSobs(source)%nr/)),&
                'Error in nf90_get_var SWdown')
           call LVT_verify(nf90_get_var(ftn,lwdownid,lwdown,&
                start=(/1,1/),&
                count=(/HydroSCSobs(source)%nc,HydroSCSobs(source)%nr/)),&
                'Error in nf90_get_var LWdown')
           call LVT_verify(nf90_get_var(ftn,windid,wind,&
                start=(/1,1/),&
                count=(/HydroSCSobs(source)%nc,HydroSCSobs(source)%nr/)),&
                'Error in nf90_get_var Wind')
           call LVT_verify(nf90_get_var(ftn,psurfid,psurf,&
                start=(/1,1/),&
                count=(/HydroSCSobs(source)%nc,HydroSCSobs(source)%nr/)),&
                'Error in nf90_get_var Psurf')
           call LVT_verify(nf90_get_var(ftn,rainfid,rainf,&
                start=(/1,1/),&
                count=(/HydroSCSobs(source)%nc,HydroSCSobs(source)%nr/)),&
                'Error in nf90_get_var Rainf')

           call LVT_verify(nf90_close(ftn))

           call interp_HydroSCSvar2d(source,lwdown,lwdown_ip)
           call interp_HydroSCSvar2d(source,swdown,swdown_ip)
           call interp_HydroSCSvar2d(source,rainf,rainf_ip)
           call interp_HydroSCSvar2d(source,qair,qair_ip)
           call interp_HydroSCSvar2d(source,wind,wind_ip)
           call interp_HydroSCSvar2d(source,tair,tair_ip)
           call interp_HydroSCSvar2d(source,psurf,psurf_ip)

        endif
        
     endif
     
  endif
#endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_LWDOWNFORC,source,lwdown_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_SWDOWNFORC,source,swdown_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_RAINFFORC,source,rainf_ip,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_QAIRFORC,source,qair_ip,&
       vlevel=1,units="kg/kg")
  call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,source,tair_ip,&
       vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_PSURFFORC,source,psurf_ip,&
       vlevel=1,units="Pa")
  call LVT_logSingleDataStreamVar(LVT_MOC_WINDFORC,source,wind_ip,&
       vlevel=1,units="m/s")   !CBB added
end subroutine readHydroSCSObs


!BOP
!
! !ROUTINE: interp_HydroSCSvar2d
! \label{interp_HydroSCSvar2d}
!
! !INTERFACE: 
subroutine interp_HydroSCSvar2d(source, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use HydroSCSobsMod
! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(HydroSCSobs(source)%nc,HydroSCSobs(source)%nr)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the HydroSCS fields to the 
!  target LVT domain
!
!EOP

  real              :: var_inp_1d(HydroSCSobs(source)%nc*HydroSCSobs(source)%nr)
  logical*1         :: input_bitmap(HydroSCSobs(source)%nc*HydroSCSobs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r
  integer           :: iret

  nc = HydroSCSobs(source)%nc
  nr = HydroSCSobs(source)%nr
  
  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        if(var_inp(c,r).ne.1e20) then 
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
           input_bitmap(c+(r-1)*nc) = .true. 
        else
           var_inp(c,r) = LVT_rc%udef
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
        endif
     enddo
  enddo
  
  if(LVT_isAtAfinerResolution(HydroSCSobs(source)%datares)) then
     call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d, &
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          HydroSCSobs(source)%rlat, & 
          HydroSCSobs(source)%rlon, &
          HydroSCSobs(source)%n11, &
          LVT_rc%udef, iret)
     
  else
     call upscaleByAveraging(&
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          HydroSCSobs(source)%n11, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d)
     
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
           var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
        endif
     enddo
  enddo

end subroutine interp_HydroSCSvar2d

!BOP
!
! !ROUTINE: interp_HydroSCSvar3d
! \label{interp_HydroSCSvar3d}
!
! !INTERFACE: 
subroutine interp_HydroSCSvar3d(source, nlevs, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use HydroSCSobsMod
! !ARGUMENTS: 
  integer           :: source
  integer           :: nlevs
  real              :: var_inp(HydroSCSobs(source)%nc,HydroSCSobs(source)%nr,nlevs)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr,nlevs)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the HydroSCS fields to the 
!  target LVT domain
!
!EOP
  integer           :: k
  real              :: var_inp_1d(HydroSCSobs(source)%nc*HydroSCSobs(source)%nr)
  logical*1         :: input_bitmap(HydroSCSobs(source)%nc*HydroSCSobs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r
  integer           :: iret

  do k=1,nlevs

     nc = HydroSCSobs(source)%nc
     nr = HydroSCSobs(source)%nr
     
     input_bitmap = .false. 
     do r=1,nr
        do c=1,nc
           if(var_inp(c,r,k).ne.1e20) then 
              var_inp_1d(c+(r-1)*nc) = var_inp(c,r,k)
              input_bitmap(c+(r-1)*nc) = .true. 
           else
              var_inp(c,r,k) = LVT_rc%udef
              var_inp_1d(c+(r-1)*nc) = var_inp(c,r,k)
           endif
        enddo
     enddo
     
     if(LVT_isAtAfinerResolution(HydroSCSobs(source)%datares)) then
        call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
             var_inp_1d, output_bitmap, var_out_1d, &
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, &
             HydroSCSobs(source)%rlat, & 
             HydroSCSobs(source)%rlon, &
             HydroSCSobs(source)%n11, &
             LVT_rc%udef, iret)
        
     else
        call upscaleByAveraging(&
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
             HydroSCSobs(source)%n11, input_bitmap, &
             var_inp_1d, output_bitmap, var_out_1d)
        
     endif
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
              var_out(c,r,k) = var_out_1d(c+(r-1)*LVT_rc%lnc)
           endif
        enddo
     enddo
  enddo
     
end subroutine interp_HydroSCSvar3d

!BOP
! 
! !ROUTINE: create_HydroSCS_filename
! \label{create_HydroSCS_filename}
!
! !INTERFACE: 
subroutine create_HydroSCS_filename(odir,yr,mo,da,hr,filename)
! 
! !USES:   
  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr, mo, da, hr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the HydroSCS data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            HydroSCS base directory
!   \item[model\_name]     name of the model used in the HydroSCS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the HydroSCS file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*2             :: fmo, fda, fhr

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
 
  !CBB removing trims 
  filename = trim(odir)//'/MERRA2_CERES_IMERG_'//fyr//fmo//fda//fhr//'00.d01.nc'
  
end subroutine create_HydroSCS_filename



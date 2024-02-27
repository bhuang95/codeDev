  SUBROUTINE write_iodav3_viirsaod(outfile)

    USE, INTRINSIC :: iso_fortran_env

    USE module_viirs_vars, ONLY: inst,sat,retrieval_type,&
         &nchans,channels,viirs_wavelength,viirstimestr,&
         &viirs_aod_output,nobs_in,nobs_out,viirstdiff
   
    USE module_viirs2aeronet_args, ONLY: infile,&
        &validtimestr,validtime,viirs_errors
  
    USE module_constants

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: outfile

    !CHARACTER(len=10) :: validtimestr
    INTEGER :: ncid ! netCDF file ID
    INTEGER :: nlocsid, nchansid, ndtimeid ! dimension IDs
    CHARACTER(len=max_varname_length) :: attname,attvalue
    
    INTEGER :: varids(50),grpids(50)

    INTEGER :: i,j,nlocs,ndtime
    INTEGER :: validtimeint,yyyy,mm,dd,hh,min,sec

    TYPE(datetime_type) :: datatime
    TYPE(timedelta_type) dt
    
    REAL, DIMENSION(nchans) :: freqs

    REAL(real32), ALLOCATABLE :: zeros_r32(:)
    REAL(real64), ALLOCATABLE :: zeros_r64(:)
    INTEGER(int32), ALLOCATABLE :: zeros_i32(:),series(:)
    INTEGER(int64), ALLOCATABLE :: zeros_i64(:),secs(:)

    REAL(real32)    :: r32_missing
    REAL(real64)    :: r64_missing
    INTEGER(int32)  :: i32_missing
    INTEGER(int64)  :: i64_missing
    
    i32_missing = -9999
    i64_missing = -9999
    r32_missing = -9999.
    r64_missing = -9999.

    CALL check_nc(nf90_create(path=outfile,cmode=nf90_netcdf4,&
         &ncid=ncid))

    CALL check_nc(nf90_put_att(ncid,NF90_GLOBAL,"_ioda_layout","ObsGroup"))
    CALL check_nc(nf90_put_att(ncid,NF90_GLOBAL,"_ioda_layout_version",0))
    CALL check_nc(nf90_put_att(ncid,NF90_GLOBAL,"observation_type","Aod"))

    READ(validtimestr,"(i)") validtimeint
    CALL check_nc(nf90_put_att(ncid,NF90_GLOBAL,"date_time",validtimeint))
    CALL check_nc(nf90_put_att(ncid,NF90_GLOBAL,"satellite",sat))
    CALL check_nc(nf90_put_att(ncid,NF90_GLOBAL,"sensor",inst))
    CALL check_nc(nf90_put_att(ncid,NF90_GLOBAL,"retrieval_type",&
         &TRIM(retrieval_type)))

    nlocs=nobs_out

    ALLOCATE(zeros_r32(nlocs),zeros_r64(nlocs),&
         &zeros_i32(nlocs),zeros_i64(nlocs),secs(nlocs))

    zeros_r32=0.
    zeros_r64=0.
    zeros_i32=0
    zeros_i64=0
    series=(/(i,i=1,nlocs)/)
    secs=int(viirstdiff)

    CALL check_nc(nf90_def_dim(ncid,'Channel',nchans,nchansid)) 
    CALL check_nc(nf90_def_dim(ncid,'Location',nlocs,nlocsid))

    CALL check_nc(nf90_def_var(ncid,'Channel',nf90_int,nchansid,&
         &varids(1)))
    CALL check_nc(nf90_def_var(ncid,'Location',nf90_int,nlocsid,&
         &varids(2)))

    CALL check_nc(nf90_def_grp(ncid,'MetaData',grpids(1)))

    CALL check_nc(nf90_def_var(grpids(1),'dateTime',nf90_int64,&
         &(/nlocsid/),varids(3)))
    attvalue="seconds since "//validtime%isoformat()
    attname="units"
    CALL check_nc(nf90_put_att(grpids(1), varids(3), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(1),varids(3),0,i64_missing))

    CALL check_nc(nf90_def_var(grpids(1),'latitude',nf90_float,&
         &(/nlocsid/),varids(4)))
    attvalue="degrees_north"
    attname="units"
    CALL check_nc(nf90_put_att(grpids(1), varids(4), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(1),varids(4),0,r32_missing))

    CALL check_nc(nf90_def_var(grpids(1),'longitude',nf90_float,&
         &(/nlocsid/),varids(5)))
    attvalue="degrees_east"
    attname="units"
    CALL check_nc(nf90_put_att(grpids(1), varids(5), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(1),varids(5),0,r32_missing))

    CALL check_nc(nf90_def_var(grpids(1),'sensorCentralFrequency',&
         &nf90_float, (/nchansid/),varids(6)))
    attvalue="Hz"
    attname="units"
    CALL check_nc(nf90_put_att(grpids(1), varids(6), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(1),varids(6),0,r32_missing))

    CALL check_nc(nf90_def_var(grpids(1),'sensorCentralWavelength',&
         &nf90_float, (/nchansid/),varids(7)))
    attvalue="microns"
    attname="units"
    CALL check_nc(nf90_put_att(grpids(1), varids(7), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(1),varids(7),0,r32_missing))

    CALL check_nc(nf90_def_var(grpids(1),'sensorChannelNumber',&
         &nf90_int, (/nchansid/),varids(8)))
    attvalue=""
    attname="units"
    CALL check_nc(nf90_put_att(grpids(1), varids(8), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(1),varids(8),0,i32_missing))

    CALL check_nc(nf90_def_grp(ncid,'ObsError',grpids(2)))

    CALL check_nc(nf90_def_var(grpids(2),'aerosolOpticalDepth',&
         &nf90_float, (/nchansid,nlocsid/),varids(9)))
    attvalue="1"
    attname="units"
    CALL check_nc(nf90_put_att(grpids(2), varids(9), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(2),varids(9),0,r32_missing))

    CALL check_nc(nf90_def_grp(ncid,'ObsValue',grpids(3)))

    CALL check_nc(nf90_def_var(grpids(3),'aerosolOpticalDepth',&
         &nf90_float, (/nchansid,nlocsid/),varids(10)))
    attvalue="1"
    attname="units"
    CALL check_nc(nf90_put_att(grpids(3), varids(10), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(3),varids(10),0,r32_missing))

    CALL check_nc(nf90_def_grp(ncid,'PreQC',grpids(4)))

    CALL check_nc(nf90_def_var(grpids(4),'aerosolOpticalDepth',&
         &nf90_int, (/nchansid,nlocsid/),varids(11)))
    attvalue=""
    attname="units"
    CALL check_nc(nf90_put_att(grpids(4), varids(11), &
         &TRIM(attname),attvalue))
    CALL check_nc(nf90_def_var_fill(grpids(4),varids(11),0,i32_missing))

    CALL check_nc(nf90_enddef(ncid))

    freqs = speed_of_light/viirs_wavelength

    CALL check_nc(nf90_put_var(ncid,varids(1),channels))
    CALL check_nc(nf90_put_var(ncid,varids(2),zeros_i32))

    CALL check_nc(nf90_put_var(grpids(1),varids(3),secs))
    CALL check_nc(nf90_put_var(grpids(1),varids(4),&
         &viirs_aod_output(:)%lat))
    CALL check_nc(nf90_put_var(grpids(1),varids(5),&
         &viirs_aod_output(:)%lon))
    CALL check_nc(nf90_put_var(grpids(1),varids(6),&
         &freqs))
    CALL check_nc(nf90_put_var(grpids(1),varids(7),&
         &(/viirs_wavelength*1.e6/)))
    CALL check_nc(nf90_put_var(grpids(1),varids(8),&
         &channels))

    CALL check_nc(nf90_put_var(grpids(2),varids(9),&
         &RESHAPE(viirs_aod_output(:)%uncertainty,(/1,nlocs/))))

    CALL check_nc(nf90_put_var(grpids(3),varids(10),&
         &RESHAPE(viirs_aod_output(:)%value550,(/1,nlocs/))))

    CALL check_nc(nf90_put_var(grpids(4),varids(11),&
         &RESHAPE(viirs_aod_output(:)%qcall,(/1,nlocs/))))

    CALL check_nc(nf90_close(ncid))

    PRINT *, 'Wrote output to ', TRIM(outfile)

  END SUBROUTINE write_iodav3_viirsaod

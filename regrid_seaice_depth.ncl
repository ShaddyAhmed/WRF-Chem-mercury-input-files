;----------------------------------------------------------------------
; Script to regrid sea ice depth data from the Arctic Ocean Physics
; Reanalysis to the WRF grid. Data for sea ice depth can be retrieved from
; the Copernicus Marine Service at: 
; https://data.marine.copernicus.eu/product/ARCTIC_MULTIYEAR_PHY_002_003/description
; Sea ice depth data is on a curvilinear grid and is regridded to the
; WRF grid using the ESMF regrid function in NCL.
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin
  
    setfileoption("nc","Format","NetCDF4")

;--Data file containing source grid (Sea ice data)
    src_file = "/data/sahmed/OBS_DATA/ArcticSIT/cmems_mod_arc_phy_my_topaz4_P1D-m_1651157543957.nc"
    sfile    = addfile(src_file,"r")

    src_lat = sfile->latitude
    src_lon = sfile->longitude

;---Get variable to regrid
    varname = "sithick"                 
    var     = sfile->$varname$         ; Assumption is that "var"
                                       ; contains coordinate arrays
    printVarSummary(var)               ; Make sure it has coord arrays.

    var@lat2d = src_lat     ; Needed for plotting.
    var@lon2d = src_lon     ; var_regrid already has these

;---Data file containing destination grid (example met_em file with WRF grid)
    dst_file = "/scratchu/sahmed/WRFChem_MOSAIC/met_em_WRFChem_MOSAIC_lower_domain_2020/met_em.d02.2020-04-28_06:00:00.nc"
    dfile    = addfile(dst_file,"r")
    dst_lat  = dfile->XLAT_M(0,:,:)                    
    dst_lon  = dfile->XLONG_M(0,:,:)                   
    seaice = dfile->ICEDEPTH
;---Set up regridding options
    Opt                   = True

;---"bilinear" is the default. "patch" and "conserve" are other options.
    Opt@InterpMethod      = "bilinear"       

;---Set file type to netcdf4 for large files
    Opt@SrcNetCDFType     = "netcdf4"
    Opt@DstNetCDFType     = "netcdf4"
    Opt@WgtNetCDFType     = "netcdf4"
    Opt@WgtFileName       = "curv_to_curv.nc"

    Opt@SrcRegional       = True              
    Opt@SrcInputFileName  = src_file         

    Opt@DstGridLat        = dst_lat           ; destination grid
    Opt@DstGridLon        = dst_lon
    Opt@DstRegional       = True             

    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True
   
    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding

    printVarSummary(var_regrid)           ; Check that everything
    printMinMax(var_regrid,0)             ; looks okay.

;---Modify variable coordinate information
    sit_regridded = var_regrid
    sit_regridded@description="Sea ice thickness"
    sit_regridded@units="m"
    sit_regridded@_FillValue = 0
    delete(sit_regridded@missing_value)
    delete(sit_regridded@remap)
    delete(sit_regridded@lon2d)
    delete(sit_regridded@lat2d)
    delete(sit_regridded@_FillValue)
    delete(sit_regridded@time)
    delete(sit_regridded@cell_methods)
    delete(sit_regridded@grid_mapping)
    delete(sit_regridded@standard_name)
    delete(sit_regridded@_ChunkSizes)

    printVarSummary(sit_regridded)

;---Create new sea ice thickness variable (ICEDEPTH)    
    icedepth = new(dimsizes(seaice),typeof(seaice))
    delete(icedepth@_FillValue)
    copy_VarMeta(seaice,icedepth)
    printVarSummary(icedepth)

;---Set indices to loop through regridded sea ice variable
    i=0    ; Indice for time value of regridded ice depth
    j=0    ; Counter to loop to the next time of the regridded ice depth

;---Save ICEDEPTH to all met_em wrf files ; NOTE: This is only for domain 2 
    met_em_dir = "/scratchu/sahmed/WRFChem_MOSAIC/met_em_WRFChem_MOSAIC_lower_domain_2020/"
    met_em_files = systemfunc ("csh -c 'cd " + met_em_dir + ";ls *d02*:00:00.nc'")
    nfiles = dimsizes(met_em_files)

;---Loop over met_em files and save ICEDEPTH
    do it = 0,nfiles-1,1
      a = addfile(met_em_dir+met_em_files(it),"w")
      print("Writing ICEDEPTH variable to " + met_em_files(it))
      if (j .eq. 3) then
        icedepth = (/ sit_regridded(i,:,:) /)
        a->ICEDEPTH = icedepth
        j=0
        i=i+1
      else
        icedepth = (/ sit_regridded(i,:,:) /)
        a->ICEDEPTH = icedepth
        j=j+1
      end if
      delete(a)
    end do

;----------------------------------------------------------------------
; Plotting section
;
; This section creates filled contour plots of both the original
; data and the regridded data, and panels them.
;----------------------------------------------------------------------
;---Do plotting routine
    PLOT_GRIDS = True

 if (PLOT_GRIDS) then
    wks = gsn_open_wks("png","curv_to_curv")

    res                       = True

    res@gsnMaximize           = True

    res@gsnDraw               = False
    res@gsnFrame              = False

    res@cnFillOn              = True
    res@cnLinesOn             = False
    res@cnLineLabelsOn        = False
    res@cnFillMode            = "RasterFill"

    res@lbLabelBarOn          = False    ; Turn on later in panel

    res@mpMinLatF             = min(src_lat)
    res@mpMaxLatF             = max(src_lat)
    res@mpMinLonF             = min(src_lon)
    res@mpMaxLonF             = max(src_lon)

;;--Change (maybe)
    mnmxint = nice_mnmxintvl( min(var), max(var), 18, False)
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF       = mnmxint(0)
    res@cnMaxLevelValF       = mnmxint(1)
    res@cnLevelSpacingF      = mnmxint(2)

;---Resources for plotting regridded data
    res@gsnAddCyclic  = False            ;;---Change (maybe)
    res@tiMainString  = "Curvilinear grid (" + Opt@InterpMethod + ")"

    plot_regrid = gsn_csm_contour_map(wks,var_regrid(0,:,:),res)

;---Resources for plotting original data
    res@gsnAddCyclic = False            ;;---Change (maybe)
    res@tiMainString = "Original curvilinear grid"

    plot_orig = gsn_csm_contour_map(wks,var(0,:,:),res)  

;---Compare the plots in a panel
    pres                   = True
    pres@gsnMaximize       = True
    pres@gsnPanelLabelBar  = True

    gsn_panel(wks,(/plot_orig,plot_regrid/),(/2,1/),pres)
 end if
end

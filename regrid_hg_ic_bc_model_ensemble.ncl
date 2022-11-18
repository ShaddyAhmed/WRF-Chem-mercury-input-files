;----------------------------------------------------------------------
; Script to regrid Hg(0) concentrations from the global model ensemble of 
; Dastoor et al., Nature Geoscience, 2022, to CAM-Chem model output. This is 
; regridded to the WRF grid by the mozbc preprocessor. Hg concentrations 
; are on a curvilinear grid and are regridded the CAM-Chem grid using the 
; ESMF regrid function in NCL.
; 
; Shaddy Ahmed 18/11/2022
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin
  
    setfileoption("nc","Format","NetCDF4")

;--Data file containing source grid
    src_file = "/data/sahmed/mosaic_Ensemble_concHg0_ModelLevel_2020_Monthly.nc"              ;;---Change (likely)
    sfile    = addfile(src_file,"r")

    src_lat = sfile->lat
    src_lon = sfile->lon
;---Get variable to regrid
    varname = "concHg0"                      ;;---Change (likely)
    var     = sfile->$varname$         ; Assumption is that "var"
                                       ; contains coordinate arrays.
    printVarSummary(var)               ; Make sure it has coord arrays.

;---Convert units from ng m-3 to mol/mol
    hg_molar_mass = 200.59  ;g/mol
    temp = 273.15  ; K
    press = 1e5 ; Pa
    r_gas = 8.314 
    air_vol = (r_gas*temp)/press 
    convf = air_vol / (hg_molar_mass * 1e9) ;from ng/m3 to mol/mol
    var = var*convf ; mol/mol

;---Data file containing destination grid
    dst_file = "/scratchu/sahmed/CAMCHEM/MOSAIC/camchem_0320_0520.nc"                      ;;---Change (likely)
    dfile    = addfile(dst_file,"r")
    dst_lat  = dfile->lat(:)                     ;;---Change (likely)
    dst_lon  = dfile->lon(:)                     ;;---Change (likely)
    dst_lev  = dfile->lev(:)
    o3 = dfile->O3
    printVarSummary(o3)
;---Set up regridding options
    Opt                   = True

;---"bilinear" is the default. "patch" and "conserve" are other options.
    Opt@InterpMethod      = "bilinear"        ;;---Change (maybe)

;---Set file type to netcdf4 for large files
    Opt@SrcNetCDFType     = "netcdf4"
    Opt@DstNetCDFType     = "netcdf4"
    Opt@WgtNetCDFType     = "netcdf4"
    Opt@WgtFileName       = "rect_to_curv.nc"

    Opt@SrcRegional       = False              ;;--Change (maybe)
    Opt@SrcInputFileName  = src_file          ; optional, but good idea

    Opt@DstGridLat        = dst_lat           ; destination grid
    Opt@DstGridLon        = dst_lon
    Opt@DstRegional       = True              ;;--Change (maybe)

    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True
   
    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding

    printVarSummary(var_regrid)           ; Check that everything
    printMinMax(var_regrid,0)             ; looks okay.

;---Create new distance to coast variable (D2CST)    
    hg = new(dimsizes(o3),typeof(o3))
    copy_VarMeta(o3,hg)
    printVarSummary(hg)
    
    ; Set data for all times in March to the March average and April to the April average for all vertical levels
    do itim = 0, 123, 1
        do ilev = 0, 55, 1
            hg(itim,ilev,:,:) = var_regrid(2,:,:)
        end do
    end do
    
    do itim = 124, 245, 1
        do ilev = 0, 55, 1
            hg(itim,ilev,:,:) = var_regrid(3,:,:)
        end do
    end do
    printVarSummary(hg)

;---Copy variable meta data and set attribute information
    hg@long_name="Hg0 concentration"
    hg@units="mol/mol"
    delete(hg@missing_value)
    delete(hg@_FillValue)
    delete(hg@standard_name)
    delete(hg@remap)
    delete(hg@lon2d)
    delete(hg@lat2d)
    printVarSummary(hg)

   ;a = addfile("/home/sahmed/hgspecfile.nc", "c")
   ;a->hg = hg

;---Save concentration to CAMCHEM file 
    camchem_file = "/scratchu/sahmed/CAMCHEM/MOSAIC/camchem_0320_0520.nc"
    a = addfile(camchem_file,"w")
    a->HG0 = hg

end

;----------------------------------------------------------------------
; Script to regrid Hg emissions (Hg0, Hg2, HgP) from the GMA 2015 inventory 
; to the WRF grid. Hg emissions are on a curvilinear grid and are regridded 
; to the WRF grid using the ESMF regrid function in NCL.
; 
; Shaddy Ahmed 18/11/2022
;----------------------------------------------------------------------


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin

; Set up summed variables
hg0_emiss = new((/720,1440/), "float")
hg0_emiss=0
hg0_emiss!0 = "lat"
hg0_emiss!1 = "lon"

hg2_emiss = new((/720,1440/), "float")
hg2_emiss=0
hg2_emiss!0 = "lat"
hg2_emiss!1 = "lon"

hgp_emiss = new((/720,1440/), "float")
hgp_emiss=0
hgp_emiss!0 = "lat"
hgp_emiss!1 = "lon"

; Loop over netcdf files and extract variable
hg_emission_dir = "/data/sahmed/HG_EMISSIONS_INVENTORY/GMA_2018_anthropogenic/"
hg0_emission_files = systemfunc ("csh -c 'cd " + hg_emission_dir + ";ls *_Hg0*.nc'")
hg2_emission_files = systemfunc ("csh -c 'cd " + hg_emission_dir + ";ls *_Hg2*.nc'")
hgp_emission_files = systemfunc ("csh -c 'cd " + hg_emission_dir + ";ls *_HgP*.nc'")
nfiles = dimsizes(hg0_emission_files)

;---Loop over Hg0 emission files
do it = 0,nfiles-1,1
  a = addfile(hg_emission_dir+hg0_emission_files(it),"r")
  print("Reading " + hg0_emission_files(it))
  hg0 = a->emi_hg_0
  hg0_emiss = hg0_emiss + hg0
  delete(a)
end do

;---Loop over Hg2 emission files
do it = 0,nfiles-1,1
  a = addfile(hg_emission_dir+hg2_emission_files(it),"r")
  print("Reading " + hg2_emission_files(it))
  hg2 = a->emi_hg_2
  hg2_emiss = hg2_emiss + hg2
  delete(a)
end do

;---Loop over HgP emission files
do it = 0,nfiles-1,1
  a = addfile(hg_emission_dir+hgp_emission_files(it),"r")
  print("Reading " + hgp_emission_files(it))
  hgp = a->emi_hg_p
  hgp_emiss = hgp_emiss + hgp
  delete(a)
end do

; Convert from kg/m2/s to mol/km2/h
sec_to_hr = 3600
m2_to_km2 = 1e6
kg_to_g = 1e3
Hg_mol_mass = 200.59
Hg2_mol_mass = 297.1

; This is in mol/km2/h
hg0_emiss = hg0_emiss * kg_to_g * m2_to_km2 * sec_to_hr / Hg_mol_mass
hg2_emiss = hg2_emiss * kg_to_g * m2_to_km2 * sec_to_hr / Hg2_mol_mass
hgp_emiss = hgp_emiss * kg_to_g * m2_to_km2 * sec_to_hr / Hg_mol_mass


;;;;;;;;;; REGRID DATA TO WRF GRID ;;;;;;;;;;;
setfileoption("nc","Format","NetCDF4")

src_file = addfile("/data/sahmed/HG_EMISSIONS_INVENTORY/GMA_2018_anthropogenic/GMA_emissions_ASGM_Hg0.0.25x0.25.2015.nc", "r")
;--Data file containing source grid
src_lat = src_file->lat
src_lon = src_file->lon

;--Assign coordinate variables
hg0_emiss&lat = src_lat
hg0_emiss&lon = src_lon
hg2_emiss&lat = src_lat
hg2_emiss&lon = src_lon
hgp_emiss&lat = src_lat
hgp_emiss&lon = src_lon

;--Assign units to variables
hg0_emiss&lat@units = "degrees_north"
hg0_emiss&lon@units = "degrees_east"
hg2_emiss&lat@units = "degrees_north"
hg2_emiss&lon@units = "degrees_east"
hgp_emiss&lat@units = "degrees_north"
hgp_emiss&lon@units = "degrees_east"

;---Get one variable to regrid
var     = hg0_emiss
printVarSummary(var)               ; Make sure it has coord arrays.

;---Data file containing destination grid
dst_file = "/scratchu/sahmed/WRFChem_MOSAIC/met_em_WRFChem_MOSAIC_lower_domain_2020/met_em.d02.2020-04-29_18:00:00.nc"
dfile    = addfile(dst_file,"r")
dst_lat  = dfile->XLAT_M(0,:,:)                   
dst_lon  = dfile->XLONG_M(0,:,:)
icedepth = dfile->ICEDEPTH

;---Set up regridding options
Opt                   = True

;---"bilinear" is the default. "patch" and "conserve" are other options.
Opt@InterpMethod      = "bilinear"   

;---Set file type to netcdf4 for large files
Opt@SrcNetCDFType     = "netcdf4"
Opt@DstNetCDFType     = "netcdf4"
Opt@WgtNetCDFType     = "netcdf4"
Opt@WgtFileName       = "rect_to_curv.nc"

Opt@SrcRegional       = False            

Opt@DstGridLat        = dst_lat           ; destination grid
Opt@DstGridLon        = dst_lon
Opt@DstRegional       = True            

Opt@ForceOverwrite    = True
Opt@PrintTimings      = True
Opt@Debug             = True

hg0_emiss_regrid = ESMF_regrid(var,Opt)     ; Do the regridding

printVarSummary(hg0_emiss_regrid)           ; Check that everything
printMinMax(hg0_emiss_regrid,0)             ; looks okay.

; Do regridding of other variables
; Hg2
hg2_emiss_regrid = ESMF_regrid_with_weights(hg2_emiss,"rect_to_curv.nc",False)
; HgP
hgp_emiss_regrid = ESMF_regrid_with_weights(hgp_emiss,"rect_to_curv.nc",False)


;--Save to separate netcdf file to check everything looks okay
a = addfile("/scratchu/sahmed/hgemissnetcdf.nc", "c")
a->hg0_emiss = hg0_emiss_regrid
a->hg2_emiss = hg2_emiss_regrid
a->hgp_emiss = hgp_emiss_regrid


;---Copy variable meta data and set attribute information from wrfchemi files
emis_file = "/scratchu/sahmed/WRFChem_MOSAIC/snow_test/wrfchemi_d01_2020-03-28_00:00:00/"
emisfile    = addfile(emis_file,"r")
e_co = emisfile->E_CO

;---Create new emission variable (E_HG)    
e_hg = new(dimsizes(e_co),typeof(e_co))
e_hg2 = new(dimsizes(e_co),typeof(e_co))
e_hgp = new(dimsizes(e_co),typeof(e_co))

;--Copy metadata
e_hg = hg0_emiss_regrid
copy_VarMeta(e_co,e_hg)
e_hg2 = hg2_emiss_regrid
copy_VarMeta(e_co,e_hg2)
e_hgp = hgp_emiss_regrid
copy_VarMeta(e_co,e_hgp)

;--Assign coordinate variables
e_hg@description="EMISSIONS"
e_hg@units="mole km-2 hr-1"
delete(e_hg@missing_value)
delete(e_hg@remap)
delete(e_hg@lon2d)
delete(e_hg@lat2d)
delete(e_hg@_FillValue)
printVarSummary(e_hg)

e_hg2@description="EMISSIONS"
e_hg2@units="mole km-2 hr-1"
delete(e_hg2@missing_value)
delete(e_hg2@remap)
delete(e_hg2@lon2d)
delete(e_hg2@lat2d)
delete(e_hg2@_FillValue)

e_hgp@description="EMISSIONS"
e_hgp@units="mole km-2 hr-1"
delete(e_hgp@missing_value)
delete(e_hgp@remap)
delete(e_hgp@lon2d)
delete(e_hgp@lat2d)
delete(e_hgp@_FillValue)

;---Save Hg emission variables to all wrfchemi files 
real_dir = "/scratchu/sahmed/WRFChem_MOSAIC/snow_test/"
chemi_files = systemfunc ("csh -c 'cd " + real_dir + ";ls wrfchemi_d01*:00:00'")
nfiles = dimsizes(chemi_files)

;---Loop over met_em files and save E_HG0, E_HG2, and E_HGP
do it = 0,nfiles-1,1
   a = addfile(real_dir+chemi_files(it),"w")
   print("Writing HG emissions to " + chemi_files(it))
   a->E_HG0 = e_hg
   a->E_HG2 = e_hg2
   a->E_HGP = e_hgp
   delete(a)
end do

;---Do plotting routine
PLOT_GRIDS = True

if (PLOT_GRIDS) then

; MAKE PLOT OF REGRIDDING
   wks = gsn_open_wks("X11","HG0_EMISSIONS")

   res                       = True
   res@gsnMaximize           = True

   res@gsnDraw               = False
   res@gsnFrame              = False

   res@cnFillOn              = True
   res@cnLinesOn             = False
   res@cnLineLabelsOn        = False
   res@cnFillMode            = "RasterFill"

   res@lbLabelBarOn          = False    ; Turn on later in panel

   res@mpMinLatF             = min(dst_lat)
   res@mpMaxLatF             = max(dst_lat)
   res@mpMinLonF             = min(dst_lon)
   res@mpMaxLonF             = max(dst_lon)

   res@cnLevelSelectionMode         = "ExplicitLevels"
   clevels                         = (/1e-15,5e-15,1e-14,3e-14,5e-14,7e-14,1e-13,3e-13,5e-13,7e-13,1e-12,3e-12,5e-12,7e-12,1e-11/)
;clevels                         = (/0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000,5000/)
   res@cnLevels                    = clevels

   ;---Resources for plotting regridded data
   res@gsnAddCyclic  = False            ;;---Change (maybe)
   res@tiMainString  = "Regridded Hg emission inventory"

   plot_regrid = gsn_csm_contour_map(wks,hg0_emiss_regrid(:,:),res)

   ;---Resources for plotting original data
   res@gsnAddCyclic = False            ;;---Change (maybe)
   res@tiMainString = "Original Hg emission inventory"

   plot_orig = gsn_csm_contour_map(wks,hg0_emiss(:,:),res)

   ;---Compare the plots in a panel
   pres                   = True
   pres@gsnMaximize       = True
   pres@gsnPanelLabelBar  = True

   gsn_panel(wks,(/plot_orig,plot_regrid/),(/2,1/),pres)

end if

end

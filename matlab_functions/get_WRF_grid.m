function [moad_cen_lat, truelat1, truelat2, stdlon, imax, jmax, kmax, dx, dy, ref_lat, ref_lon, map_proj, hemi, ref_x, ref_y] = get_WRF_grid(WRFOUT_FILENAME)

% This routine retrieves the WRF grid properties/parameters from a wrfout file

ncid = netcdf.open(WRFOUT_FILENAME, 'NC_NOWRITE');

moad_cen_lat = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'MOAD_CEN_LAT'));
truelat1 = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'TRUELAT1'));
truelat2 = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'TRUELAT2'));
stdlon = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'STAND_LON'));
imax = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'WEST-EAST_GRID_DIMENSION')-1;
jmax = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'SOUTH-NORTH_GRID_DIMENSION')-1;
kmax = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'BOTTOM-TOP_GRID_DIMENSION')-1;
dx = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'DX')) / 1000; % km
dy = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'DY')) / 1000; % km
ref_lat = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'CEN_LAT'));
ref_lon = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'CEN_LON'));
map_proj = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'MAP_PROJ');
if(truelat1 > 0)
    hemi = 1;
else
    hemi = -1;
end

ref_x = double((imax + 1) / 2);
ref_y = double((jmax + 1) / 2);

netcdf.close(ncid);

end

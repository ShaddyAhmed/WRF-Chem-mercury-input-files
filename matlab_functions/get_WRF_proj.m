function [proj] = get_WRF_proj(WRF_FILENAME)
%
% 2017/11/22, Louis Marelle
%
%-------- This routine retrieves the WRF grid properties/parameters from a WRF netCDF file (e.g. wrfout* file, wrfinput*, met_em*...) --------
%
% This should be used instead of get_WRF_grid, which should be deprecated at some point, since this is a lot more flexible for adding new output fields
%
% Inputs:
%   WRF_FILE, full path to a WRF netCDF file
% Outputs:
%   proj, a matlab structure containing WRF projection details (e.g., proj.imax for the west-east grid dimension)

% Open WRF file WRF_FILENAME
ncid = netcdf.open(WRF_FILENAME, 'NC_NOWRITE');
proj.map_proj = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'MAP_PROJ');
proj.imax = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'WEST-EAST_GRID_DIMENSION') - 1);
proj.jmax = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'SOUTH-NORTH_GRID_DIMENSION') - 1);
proj.kmax = double(netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'BOTTOM-TOP_GRID_DIMENSION') - 1);
proj.dx = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'DX') / 1000; % km
proj.dy = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'DY') / 1000; % km
proj.moad_cen_lat = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'MOAD_CEN_LAT');
proj.truelat1 = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'TRUELAT1');
proj.truelat2 = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'TRUELAT2');
proj.stdlon = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'STAND_LON');
proj.ref_lat = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'CEN_LAT');
proj.ref_lon = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'CEN_LON');
proj.pole_lat = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'POLE_LAT');
proj.pole_lon = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'POLE_LON');
proj.mminlu = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'MMINLU');
proj.iswater = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'ISWATER');
proj.islake = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'ISLAKE');
proj.isice = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'ISICE');
proj.isurban = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'ISURBAN');
proj.isoilwater = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'ISOILWATER');
% Get info about the parent nest, if applicable
try
  proj.parent_id  = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'PARENT_ID');
  proj.i_parent_start = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'I_PARENT_START');
  proj.j_parent_start = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'J_PARENT_START');
  proj.parent_grid_ratio = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'PARENT_GRID_RATIO');
catch
  proj.parent_id  = NaN;
  proj.i_parent_start = NaN;
  proj.j_parent_start = NaN;
  proj.parent_grid_ratio = NaN;
end

netcdf.close(ncid);

% Calculate hemisphere
% Not sure if this holds for all projections
if(proj.truelat1 > 0)
  proj.hemi = 1;
else
  proj.hemi = -1;
end

% Calulate the reference (= known; = center) WRF x,y (corresponding to ref_lat and ref_lon)
proj.ref_x = double((proj.imax + 1) / 2);
proj.ref_y = double((proj.jmax + 1) / 2);

end


function []=create_netcdffile_saprc_mosaic_dms(wrfchemi_filename, specnames, timestring, zdim, dim_we, dim_sn, dx, dy, cenlat, cenlon, truelat1, truelat2, moad_cen_lat, mapproj, mminlu, COMMENT_EMISSIONS)
% ncfiles[i],timestring[i],wrflevs,wrf_ix,wrf_jx,dx*1e3,dx*1e3,lat1,lon1,truelat1,truelat2,moad_cen_lat,2,nei_path


try
ncid_wrfchemi = netcdf.create(wrfchemi_filename,'NC_NOCLOBBER');
catch exception
system(['rm -f ' wrfchemi_filename]);
ncid_wrfchemi = netcdf.create(wrfchemi_filename,'NC_NOCLOBBER');
end

% create dimensions
dimid_time = netcdf.defDim(ncid_wrfchemi,'Time',1);
dimid_datestrlen = netcdf.defDim(ncid_wrfchemi, 'DateStrLen',19); 
dimid_we = netcdf.defDim(ncid_wrfchemi,'west_east',dim_we);
dimid_sn = netcdf.defDim(ncid_wrfchemi,'south_north',dim_sn);
dimid_bt = netcdf.defDim(ncid_wrfchemi,'bottom_top',zdim); % !!!
dimid_stag = netcdf.defDim(ncid_wrfchemi,'emissions_zdim_stag',zdim);

%netcdf
% write global attributes
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'WEST-EAST_GRID_DIMENSION',dim_we+1)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'SOUTH-NORTH_GRID_DIMENSION',dim_sn+1)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'BOTTOM-TOP_GRID_DIMENSION',int32(zdim)) % !!!
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'TITLE', COMMENT_EMISSIONS)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'DX',dx*1000)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'DY',dy*1000)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'CEN_LAT',cenlat)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'CEN_LON',cenlon)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'TRUELAT1',truelat1)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'TRUELAT2',truelat2)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'MOAD_CEN_LAT',moad_cen_lat)
netcdf.putAtt(ncid_wrfchemi,netcdf.getConstant('NC_GLOBAL'),'MAP_PROJ',mapproj)
netcdf.putAtt(ncid_wrfchemi, netcdf.getConstant('NC_GLOBAL'), 'MMINLU', mminlu)

% specnames = {'CO','NH3','SO2','NO','C2H5OH','CH3OH','HC5','ETH','TOL','OL2','KET','XYL','HCHO','OLT','ALD','OLI','CSL','ORGJ','ECJ'};
names=strcat('E_',specnames);
nnames=length(names);


% define varaiable fields
tvarid = netcdf.defVar(ncid_wrfchemi,'Times','char',[dimid_datestrlen, dimid_time]);
for n=1:nnames
    if (strcmp(specnames{n},'ECJ') || strcmp(specnames{n},'ORGJ') || strcmp(specnames{n},'SO4J')|| strcmp(specnames{n},'PM25J'))
        varid = netcdf.defVar(ncid_wrfchemi,names{n},'float',[dimid_we,dimid_sn,dimid_stag,dimid_time ]);
        netcdf.putAtt(ncid_wrfchemi,varid,'FieldType',int32(104))
        netcdf.putAtt(ncid_wrfchemi,varid,'MemoryOrder','XYZ')
        netcdf.putAtt(ncid_wrfchemi,varid,'description','EMISSIONS')
        netcdf.putAtt(ncid_wrfchemi,varid,'units','ug m-2 s-1')
        netcdf.putAtt(ncid_wrfchemi,varid,'stagger','Z')
    elseif(strcmp(specnames{n},'DMS_OC'))
        varid = netcdf.defVar(ncid_wrfchemi,names{n},'float',[dimid_we,dimid_sn,dimid_stag,dimid_time ]);
        netcdf.putAtt(ncid_wrfchemi,varid,'FieldType',int32(104))
        netcdf.putAtt(ncid_wrfchemi,varid,'MemoryOrder','XYZ')
        netcdf.putAtt(ncid_wrfchemi,varid,'description','OCEANIC DMS CONTENT')
        netcdf.putAtt(ncid_wrfchemi,varid,'units','mole m^-3')
        netcdf.putAtt(ncid_wrfchemi,varid,'stagger','Z')
    else
        varid = netcdf.defVar(ncid_wrfchemi,names{n},'float',[dimid_we,dimid_sn,dimid_stag,dimid_time ]);
        netcdf.putAtt(ncid_wrfchemi,varid,'FieldType',int32(104))
        netcdf.putAtt(ncid_wrfchemi,varid,'MemoryOrder','XYZ')
        netcdf.putAtt(ncid_wrfchemi,varid,'description','EMISSIONS')
        netcdf.putAtt(ncid_wrfchemi,varid,'units','mole km-2 hr-1')
        netcdf.putAtt(ncid_wrfchemi,varid,'stagger','Z')
    end
end

% Leave define mode and enter data mode
netcdf.endDef(ncid_wrfchemi)
% write fields
netcdf.putVar(ncid_wrfchemi,tvarid,timestring);

%Close netcdf file
netcdf.close(ncid_wrfchemi);


end

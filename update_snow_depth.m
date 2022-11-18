% Shaddy Ahmed, Louis Marelle, 2022/11/18
%
clear
close all

%-------- Input --------
WPS_DIR = '/scratchu/sahmed/WRFChem_MOSAIC/snow_test/'


%-------- Parameters --------
SNOWDEPTH_SEAICE = 0.3; % meters, N-ICE campaign, before the melt season begins
ICEDEPTH_SEAICE = 1.5; % meters, if USE_ASR2_SEAICE = false
ALBEDO_SEAICE = 0.82; % N-ICE, before the melt season begins
SNOWDENSITY_SEAICE = 200.0; % kg/m3

%-------- Initialize --------

for idomain = 1:2
  domain = ['d0', num2str(idomain)];
  disp(domain)

  % Get met_em dates from the met_em_<domain>* filenames in WRFCHEMI_PATH/
  path_files = dir([WPS_DIR, '/met_em.', domain, '*.nc']);
  met_em_filenames = { path_files.name };
  met_em_filenames = cell2mat(met_em_filenames');
  met_em_dates = met_em_filenames(:, end-21:end-3);
  met_em_dates = datenum(met_em_dates, 'yyyy-mm-dd_HH:MM:SS');
  years_list = year(met_em_dates)';
  months_list = month(met_em_dates)';
  days_list = day(met_em_dates)';
  hours_list = hour(met_em_dates)';

  % Loop on days
  ndates = length(met_em_dates);
  for idate = 1:ndates
    met_em_date = met_em_dates(idate);
    disp(['  ', num2str(idate), ' ', datestr(met_em_date)])

    % Open WRF met_em SEAICE for the current date
    met_em_filename = [WPS_DIR, '/met_em.', domain, '.', datestr(met_em_date, 'yyyy-mm-dd_HH:MM:SS'), '.nc'];
    ncid = netcdf.open(met_em_filename, 'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid, 'SEAICE');
    wrf_seaice = netcdf.getVar(ncid, varid);
    varid = netcdf.inqVarID(ncid, 'SNOWH');
    wrf_snowh = netcdf.getVar(ncid, varid);
    varid = netcdf.inqVarID(ncid, 'SNOW');
    wrf_snow = netcdf.getVar(ncid, varid);
    netcdf.close(ncid);

    % Set ICEDEPTH to the fixed value of ICEDEPTH_SEAICE
    wrf_icedepth = zeros(size(wrf_seaice));
    wrf_icedepth(wrf_seaice > 0.0) = ICEDEPTH_SEAICE;

    % Set SNOWH on seaice to the value of SNOWDEPTH_SEAICE
    wrf_snowh(wrf_seaice > 0.0) = SNOWDEPTH_SEAICE;
    wrf_snow(wrf_seaice > 0.0) = wrf_snowh(wrf_seaice > 0.0) * SNOWDENSITY_SEAICE;

    % Set ALBSI to the default value of ALBEDO_SEAICE
    wrf_albsi = zeros(size(wrf_seaice));
    wrf_albsi(wrf_seaice > 0.0) = ALBEDO_SEAICE;

    % Create and write ICEDEPTH & ALBSI to all met_em files and write
    % FLAG_ICEDEPTH,FLAG_ALBSI
    ncid = netcdf.open(met_em_filename, 'NC_WRITE');
    % If the ICEDEPTH variable already exists do not create it
    try
      varid = netcdf.inqVarID(ncid, 'ICEDEPTH');
    catch
      netcdf.reDef(ncid)
      dimid_we = netcdf.inqDimID(ncid, 'west_east');
      dimid_sn = netcdf.inqDimID(ncid, 'south_north');
      dimid_time = netcdf.inqDimID(ncid, 'Time');
      varid = netcdf.defVar(ncid, 'ICEDEPTH', 'float', [dimid_we, dimid_sn, dimid_time]);
      varid = netcdf.inqVarID(ncid, 'ICEDEPTH');
      netcdf.putAtt(ncid, varid, 'FieldType', int32(104))
      netcdf.putAtt(ncid, varid, 'MemoryOrder', 'XY')
      netcdf.putAtt(ncid, varid, 'description', 'ICE DEPTH')
      netcdf.putAtt(ncid, varid, 'units', 'm')
      netcdf.putAtt(ncid, varid, 'stagger', '')
      netcdf.putAtt(ncid, varid, 'coordinates', 'XLONG XLAT')
      netcdf.endDef(ncid)
    end
    % If the ALBSI variable exists do not create it
    try
      varid = netcdf.inqVarID(ncid, 'ALBSI');
    catch
      netcdf.reDef(ncid)
      dimid_we = netcdf.inqDimID(ncid, 'west_east');
      dimid_sn = netcdf.inqDimID(ncid, 'south_north');
      dimid_time = netcdf.inqDimID(ncid, 'Time');
      varid = netcdf.defVar(ncid, 'ALBSI', 'float', [dimid_we, dimid_sn, dimid_time]);
      varid = netcdf.inqVarID(ncid, 'ALBSI');
      netcdf.putAtt(ncid, varid, 'FieldType', int32(104))
      netcdf.putAtt(ncid, varid, 'MemoryOrder', 'XY')
      netcdf.putAtt(ncid, varid, 'description', 'SEA ICE ALBEDO')
      netcdf.putAtt(ncid, varid, 'units', ' ')
      netcdf.putAtt(ncid, varid, 'stagger', '')
      netcdf.putAtt(ncid, varid, 'coordinates', 'XLONG XLAT')
      netcdf.endDef(ncid)
    end
    netcdf.reDef(ncid)
    % Write variables and attributes to met_em files
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid, varid, 'FLAG_ICEDEPTH', int32(1));
    netcdf.putAtt(ncid, varid, 'FLAG_ALBSI', int32(1));
    netcdf.endDef(ncid)
    varid = netcdf.inqVarID(ncid, 'ICEDEPTH');
    netcdf.putVar(ncid, varid, wrf_icedepth);
    varid = netcdf.inqVarID(ncid, 'ALBSI');
    netcdf.putVar(ncid, varid, wrf_albsi);
    varid = netcdf.inqVarID(ncid, 'SNOWH');
    netcdf.putVar(ncid, varid, wrf_snowh);
    varid = netcdf.inqVarID(ncid, 'SNOW');
    netcdf.putVar(ncid, varid, wrf_snow);
    netcdf.close(ncid);

  end % for idate = 1:ndates
end % for idomain


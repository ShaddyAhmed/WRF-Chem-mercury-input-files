function [] = prep_anthro_emissions_saprcmosaic_eclipsev6_dms(startdate, enddate)

  %-------- Create WRF-Chem emissions (wrfchemi_d0*) for the SAPRCMOSAIC mechanism ----------%
  %
  % Louis Marelle, 2022/04/20
  %
  % Purpose:
  %   This MATLAB routine (function) creates daily wrfchemi_dXX input emission
  %   files for WRFChem using:
  %   - Anthropogenic emissions from the ECLIPSEv6b inventory
  %   - This routine also includes E_DMS_OC in the wrfchemi file, which is the
  %     oceanic DMS content from Lana et al. 2011 interpolated on the WRF grid.
  %     E_DMS_OC used to compute oceanic dms emissions in
  %     module_nightingale_dmsemiss, a custom online dms emission routine that I
  %     developped for the CBMZ-MOSAIC and SAPRC-MOSAIC mechanisms including some
  %     DMS chemistry. In the future, this should probably be read in a different
  %     input file.
  %
  % Inputs: 
  %  startdate: first date at which the emissions must be created, in MATLAB datenum format
  %  enddate: last date at which the emissions must be created, in MATLAB datenum format
  %
  % Outputs:
  %   Creates daily wrfchemi files in the path where this routine was run
  %
  % This routine can be modified easily to create hourly files, by replacing the 
  % "for ihour = 0:0" loop by "for ihour = 0:23", and commenting the line:
  % "hourly_factors_now = 1". There is no daily/hourly variation for shipping 
  % emisssions.
  %
  % The routine writes wrfchemi files in the path where the routine is located/run. 
  % These files can take a lot of disk space for big runs
  %
  % Always check the emissions before running WRFChem to make sure that the routine
  % produced reasonable results (I usually check one VOC, one trace gas and one 
  % aerosol). In the future I might print some diagnosis at the end of this routine
  % to make sure that emission mass is conserved and the regridding or netcdf writing
  % goes well.
  %------------------------------------------------------------------------------------------%


  %-------- Input --------
  % This needs to contain the wrfinput file, can be different from the WRF run
  RUN_DIRECTORY = './';
  max_domains = 1;


  %-------- Parameters --------
  % Emission datasets directories
  ECLIPSE_DIRECTORY = '/data/onishi/ECLIPSE_V6b';
  LANA_DMS_DIRECTORY = '/data/marelle/EMISSIONS/DMS_LANA';
  
  % Years present in the ECLIPSE inventory
  years_avail_ECLIPSE = [1990, 1995, 2000, 2005, 2008, 2009, 2010, 2014, 2015, 2016, 2020, 2025, 2030, 2040, 2050];
  
  %---- Species
  % Mechanism species, written in wrfchemi file
  SPECNAMES_SAPRC = {'CO', 'NH3', 'SO2', 'NO', 'NO2', 'C2H6', 'C3H8', 'C2H2', 'ALK3', 'ALK4',...
                     'ALK5', 'ETHENE', 'C3H6', 'OLE2', 'ARO1', 'ARO2', 'HCHO', 'CCHO', 'ACET', 'MEK',...
                     'TERP', 'MEOH', 'PROD2', 'ORGJ', 'ECJ', 'DMS_OC'};
  % Emission datasets species
  SPECNAMES_ECLIPSE = {'CO', 'NH3', 'SO2', 'NOx', 'NOx', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC',...
                       'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC',...
                       'VOC', 'VOC', 'VOC', 'OC', 'BC', 0};
  ECLIPSE_SECTORS = {'ene', 'ind', 'dom', 'tra', 'agr', 'wst', 'slv', 'shp', 'flr'};

  % Molarweights for mechanism species, NOx = NO2 in inventory
  MOLARWEIGHTS_ECLIPSE = [28, 17, 64, 46, 46, NaN, NaN, NaN, NaN, NaN,...
                          NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,...
                          NaN, NaN, NaN, NaN, NaN, NaN];

  % NOx proportion as NO
  % Proportion of NOx stored as NO for ship emissions %EPA, 2000, von glasow et
  % al for full ref
  nox_as_no_ships = 0.94;
  % Proportion of NOx stored as NO for anthropogenic emissions % ~Kaynak et al
  % 2009, Finlayson-Pitts and Pitts, 1999, google "nox "90% no" "10% no2" for
  % more references
  nox_as_no_anthro = 0.90;

  % Number of vertical levels in the wrfchemi file
  ZDIM = 1;

  % Conversion factor from OC to POM
  % 1 g oc = 1.25 g POM (value from saprc_mosaic, M. Shrivastava)
  oc_pom_factor=1.25;

  %---- Time variation factors
  % Daily cycle - hourly factors
  % Factors to apply a daily emission cycle: 1h to 24h, from TNO-MACC (Denier van
  % der Gon et al., full ref in Marelle et al., ACP, 2016)
  hourly_factors = [0.79, 0.72, 0.72, 0.71, 0.74, 0.80, 0.92, 1.08, 1.19, 1.22, 1.21, 1.21,...
                    1.17, 1.15, 1.14, 1.13, 1.10, 1.07, 1.04, 1.02, 1.02, 1.01, 0.96, 0.88;... ENERGY
                    0.75, 0.75, 0.78, 0.82, 0.88, 0.95, 1.02, 1.09, 1.16, 1.22, 1.28, 1.30,...
                    1.22, 1.24, 1.25, 1.16, 1.08, 1.01, 0.95, 0.90, 0.85, 0.81, 0.78, 0.75;... INDUSTRY
                    0.40, 0.40, 0.40, 0.40, 0.40, 0.50, 1.20, 1.50, 1.60, 1.60, 1.40, 1.20,...
                    1.10, 1.10, 1.00, 1.00, 1.00, 1.10, 1.40, 1.50, 1.40, 1.40, 1.00, 0.40;... RESIDENTIAL
                    0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84, 1.86, 1.41, 1.24, 1.20,...
                    1.32, 1.44, 1.45, 1.59, 2.03, 2.08, 1.51, 1.06, 0.74, 0.62, 0.61, 0.44;... TRANSPORT
                    0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75, 0.90, 1.10, 1.25, 1.45, 1.60,...
                    1.80, 1.75, 1.70, 1.55, 1.35, 1.10, 0.90, 0.75, 0.65, 0.60, 0.60, 0.60;... AGRICULTURE
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,...
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... WASTE = 1.00
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,...
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... SOLVENTS = 1.00
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,...
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... SHIPS = 1.00
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,...
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]; % FLARING = 1.00
  % Weekly cycle - daily factors
  % Daily factors to apply a weekly emission cycle, Sunday to Saturday from TNO-MACC
  daily_factors = [0.85, 1.06, 1.06, 1.06, 1.06, 1.06, 0.85;... ENERGY
                   0.80, 1.08, 1.08, 1.08, 1.08, 1.08, 0.80;... INDUSTRY
                   0.80, 1.08, 1.08, 1.08, 1.08, 1.08, 0.80;... RESIDENTIAL
                   0.79, 1.02, 1.06, 1.08, 1.10, 1.14, 0.81;... TRANSPORT
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... AGRICULTURE
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... WASTE = 1.00
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... SOLVENTS = 1.00
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... SHIPS = 1.00
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]; % FLARING = 1.00


  %-------- Initialization --------
  routine_filename = mfilename;
  disp(' ')
  disp(['-------- ', routine_filename, ' - produce wrfchemi files --------'])
  disp(' ')
  disp('Initializing')

  % Include my MATLAB toolbox (ijll and llij projection routines,
  % create_netcdf and speciate_voc functions)
  addpath('./matlab_functions/');

  % Cleanup existing wrfchemi files
  system(['rm -f wrfchemi_*']);

  % Year for the ECLIPSE inventory
  year_ECLIPSE = year(startdate);

  %---- Initialization: WRF emission dates
  wrf_dates = startdate:enddate;
  ndates = length(wrf_dates);
  years_list = year(wrf_dates);
  months_list = month(wrf_dates);
  days_list = day(wrf_dates);
  % List of individual simulation years
  years_list_y = years_list([true, diff(years_list)~=0]);
  % List of individual simulation months (e.g., for a run from 2012/12/1 to
  % 2013/3/1, months_list_m = [12 1 2 3])
  months_list_m = months_list([true, diff(months_list)~=0]);
  % List of years corresponding to individual simulation months (e.g., same
  % example than above, years_list_m = [2012 2013 2013 2013]
  years_list_m = years_list([true, diff(months_list)~=0]);
  % List of individual simulation days
  days_list_d = days_list([true, diff(days_list)~=0]);
  % List of months corresponding to individual simulation days
  months_list_d = months_list([true, diff(days_list)~=0]);
  % List of years corresponding to individual simulation days
  years_list_d = years_list([true, diff(days_list)~=0]);

  %---- Initialization: emission inventory years
  % Find the closest years to year_ECLIPSE in the ECLIPSE inventory
  if(year_ECLIPSE < min(years_avail_ECLIPSE) || year_ECLIPSE >= max(years_avail_ECLIPSE))
    error('Error, year_ECLIPSE is out of range')
  end
  % Find the closest year in the list. If the requested year is after the
  % closest year, pick the closest year and the next year in the list for
  % interpolation. If the requested year is before the closest year, pick the
  % closest year and the previous year in the list for interpolation
  [min_val_ECLIPSE, indx_closest_ECLIPSE_year] = min(abs(year_ECLIPSE - years_avail_ECLIPSE));
  if(min_val_ECLIPSE > 0)
    % Closest year is before
    years_closest_ECLIPSE = years_avail_ECLIPSE(indx_closest_ECLIPSE_year:indx_closest_ECLIPSE_year + 1);
  elseif(min_val_ECLIPSE < 0)
    % Closest year is after
    years_closest_ECLIPSE = years_avail_ECLIPSE(indx_closest_ECLIPSE_year - 1:indx_closest_ECLIPSE_year);
  else
    % Closest year is in the list
    years_closest_ECLIPSE = [year_ECLIPSE year_ECLIPSE];
  end
  % Compute the factors for the linear interpolation of emissions to year_ECLIPSE
  if(any(~(years_avail_ECLIPSE - year_ECLIPSE)))
    % If the requested year is available in the list, no need to interpolate
    interp_factor_before_ECLIPSE = 1;
    interp_factor_after_ECLIPSE = 0; 
  else
    % If the year is not in the list, need to interpolate between closest years
    interp_factor_before_ECLIPSE = ((year_ECLIPSE) - (years_closest_ECLIPSE(1))) ...
                          / ((years_closest_ECLIPSE(2)) - (years_closest_ECLIPSE(1)));
    interp_factor_after_ECLIPSE = ((years_closest_ECLIPSE(2)) - (year_ECLIPSE)) ... 
                          / ((years_closest_ECLIPSE(2)) - (years_closest_ECLIPSE(1)));
  end
  disp(['  Year is ', num2str(year_ECLIPSE), ', interpolating between closest years for ECLIPSE, ' num2str(years_closest_ECLIPSE)])

  %---- Initialization: get emission grid dimensions
  %-- Get ECLIPSE grid dimensions
  disp('  Reading lat & lon from ECLIPSE emission files')
  eclipse_filename = [ ECLIPSE_DIRECTORY '/ETP_base_CLE_V6_CO_' num2str(years_closest_ECLIPSE(1)) '.nc'];
  ncid = netcdf.open(eclipse_filename, 'NC_NOWRITE');
  varid = netcdf.inqVarID(ncid, 'lat');
  emissions_lat_ECLIPSE = netcdf.getVar(ncid, varid);
  varid = netcdf.inqVarID(ncid, 'lon');
  emissions_lon_ECLIPSE = netcdf.getVar(ncid, varid);
  netcdf.close(ncid);
  % Dimensions of ECLIPSE grid
  DLAT_ECLIPSE = emissions_lat_ECLIPSE(2) - emissions_lat_ECLIPSE(1);
  DLON_ECLIPSE = emissions_lon_ECLIPSE(2) - emissions_lon_ECLIPSE(1);
  imax_emiss_ECLIPSE = length(emissions_lon_ECLIPSE);
  jmax_emiss_ECLIPSE = length(emissions_lat_ECLIPSE);
  %-- Get DMS climatology grid dimensions
  disp('  Reading lat & lon for the DMS climatology')
  emissions_lat_DMS = -89.5:1:89.5;
  emissions_lon_DMS = -179.5:1:179.5;

  %---- Initialization: Check dimensions
  nspecies = length(SPECNAMES_SAPRC);
  nsectors_ECLIPSE = length(ECLIPSE_SECTORS);
  if(nspecies ~= length(SPECNAMES_ECLIPSE) | nspecies ~= length(MOLARWEIGHTS_ECLIPSE))
    error(['Error, nspecies should be the same for all inventories and for the mechanism, ',...
            num2str([nspecies, length(SPECNAMES_ECLIPSE), length(MOLARWEIGHTS_ECLIPSE)]) ])
  end


  %-------- Main prepemis_anthro program  --------
  % WRF domains loop
  for idomain = 1:max_domains
    % Domain name
    domain = ['d0', num2str(idomain)];
    disp(' ')
    disp(['WRF Domain ', domain])
    disp(['  Initializing for WRF domain ', domain])

    %---- Initializing grid for domain idomain
    % Read projection & grid info from wrfinput file
    disp('  Reading WRF domain info')
    wrfinput_file = [RUN_DIRECTORY, '/wrfinput_', domain];
    [moad_cen_lat, truelat1, truelat2, stdlon, imax, jmax, kmax, dx, dy, ref_lat,...
      ref_lon, map_proj, hemi, ref_x, ref_y] = get_WRF_grid(wrfinput_file);
    wrf_proj = get_WRF_proj(wrfinput_file);
    % Find min & max lat & lon from the WRF domain
    disp('  Opening wrfinput file to get WRF lat & lon boundaries')
    ncid = netcdf.open(wrfinput_file), 'NC_NOWRITE';
    data_lat_wrf = ncread(ncid, 'XLAT');
    data_lon_wrf = ncread(ncid, 'XLONG');
    % Adding -1 & + 1 degrees to avoid missing emissions at the domain
    % boundaries
    lat_min = min(min(data_lat_wrf)) - 1;
    lat_max = max(max(data_lat_wrf)) + 1;
    lon_min = min(min(data_lon_wrf)) - 1;
    lon_max = max(max(data_lon_wrf)) + 1;
    % Read the map factors from the wrfinput file
    data_mapfacx_wrf = ncread(ncid, 'MAPFAC_M');
    data_mapfacy_wrf = ncread(ncid, 'MAPFAC_M');
    data_landmask_wrf = ncread(ncid, 'LANDMASK');
    netcdf.close(ncid);
    % Choose ncells, number of minicells per grid length, used for regridding emissions 
    % This was determined by trial and error, and might not work in all cases, check the 
    % emissions for spurious "lines" or "waves" through the data, which indicate a too low
    % minicell number was chosen
    if(min(dx, dy) >= 100)
      ncells_ECLIPSE = 8;
    elseif(min(dx, dy) >= 50)
      ncells_ECLIPSE = 16;
    elseif ( min(dx, dy) >= 25)
      ncells_ECLIPSE = 32;
    elseif ( min(dx, dy) >= 10)
      ncells_ECLIPSE = 64;
    else
      ncells_ECLIPSE = 128;
    end % if min(dx, dy)
    disp(['  Minicell number for domain ', domain, ', ncells_ECLIPSE = ' num2str(ncells_ECLIPSE)])

    
    % -------- Emission mapping on wrf grid --------
    % The objective of this step is to create a mapping between the grid of
    % each inventory and the WRF grid
    % This mapping is later used to distribute emissions from each inventory on the WRF grid 
    disp(['  Emission mapping from emission inventory to WRF grid for domain ', domain])
  
    % Emission mapping - Creating the mapping  - ECLIPSE
    tic
    disp('    Creating the mapping, ECLIPSE')
    [mapping_array_ECLIPSE] = map_emissions(imax, jmax, imax_emiss_ECLIPSE, jmax_emiss_ECLIPSE,...
                                            emissions_lon_ECLIPSE, emissions_lat_ECLIPSE,...
                                            DLON_ECLIPSE, DLAT_ECLIPSE, lon_min, lon_max,...
                                            lat_min, lat_max, ncells_ECLIPSE, truelat1,...
                                            truelat2, hemi, stdlon, ref_lat, ref_lon,...
                                            ref_x, ref_y, dx, map_proj);
    t_elapsed = toc;
    disp(['      Elapsed time ' num2str(t_elapsed) 's']);
    

    % -------- Loop over months and interpolate monthly emission data -------- 
    for imonth = 1:length(months_list_m) 
      % months_list_m contains all the individual months (e.g., for a run from
      % 2012/12/1 to 2013/3/1, months_list_m = [12 1 2 3])
      month_i = months_list_m(imonth);
      % Years corresponding to this month (month_i)
      year_i = years_list_m(imonth);
      disp(' ')
      disp([ 'Regridding monthly emission data to WRF domain ', domain, ' for ' datestr(datenum(year_i, month_i, 1), 'yyyy/mm')]);
      
      %---- Monthly scaling factors for ECLIPSE
      % Open ECLIPSE monthly factors for current month
      ECLIPSE_MONTHS_FILENAME = [ECLIPSE_DIRECTORY '/ECLIPSE_V6a_monthly_pattern.nc'];
      ncid = netcdf.open(ECLIPSE_MONTHS_FILENAME, 'NC_NOWRITE');
      eclipse_mfac_dom = ncread(ncid, 'dom');
      eclipse_mfac_dom = eclipse_mfac_dom(:, :, month_i);
      eclipse_mfac_ene = ncread(ncid, 'ene');
      eclipse_mfac_ene = eclipse_mfac_ene(:, :, month_i);
      eclipse_mfac_ind = ncread(ncid, 'ind');
      eclipse_mfac_ind = eclipse_mfac_ind(:, :, month_i);
      eclipse_mfac_tra = ncread(ncid, 'tra');
      eclipse_mfac_tra = eclipse_mfac_tra(:, :, month_i);
      eclipse_mfac_wst = ncread(ncid, 'wst');
      eclipse_mfac_wst = eclipse_mfac_wst(:, :, month_i);
      eclipse_mfac_slv = ncread(ncid, 'slv');
      eclipse_mfac_slv = eclipse_mfac_slv(:, :, month_i);
      eclipse_mfac_shp = ncread(ncid, 'shp');
      eclipse_mfac_shp = eclipse_mfac_shp(:, :, month_i);
      eclipse_mfac_flr = ncread(ncid, 'flr');
      eclipse_mfac_flr = eclipse_mfac_flr(:, :, month_i);
      eclipse_mfac_agr = ncread(ncid, 'agr');
      eclipse_mfac_agr = eclipse_mfac_agr(:, :, month_i);
      eclipse_mfac_agr_NH3 = ncread(ncid, 'agr_NH3');
      eclipse_mfac_agr_NH3 = eclipse_mfac_agr_NH3(:, :, month_i);
      netcdf.close(ncid);
      
      %---- Regrid the monthly DMS oceanic values from Lana et al
      % Oceanic DMS: This is only done once per month because the data is monthly
      DMS_FILE = [LANA_DMS_DIRECTORY '/DMSclim_' upper(datestr(datenum(1, month_i, 1), 'mmm')) '.csv'];
      disp(['  Opening ' DMS_FILE])
      fid = fopen(DMS_FILE, 'r');
      DMS_string_format = repmat('%f', 1, 360);
      C = textscan(fid, DMS_string_format, 'delimiter', ',');
      data_DMS = flipud((cell2mat(C)))'*1E-6; %mol.m^-3
      fclose(fid);
      clear C
      data_wrf_DMS = zeros(imax, jmax);
      for ii = 1:imax
        for jj = 1:jmax
          if(~data_landmask_wrf(ii, jj))
            data_wrf_DMS(ii, jj) = interp2(emissions_lat_DMS, emissions_lon_DMS,...
                                     data_DMS, data_lat_wrf(ii, jj), data_lon_wrf(ii, jj), 'nearest');
            if(isnan(data_wrf_DMS(ii, jj)))
                
            end
          end
        end
      end
      % Monthly oceanic DMS : Some values of oceanic content at sea (coastal)
      % are still NaN because of the coarse resolution of the climatology, if
      % coastal values are NaN replace them with nearby values until all ocean
      % gridpoints are filled.  This is not a very elegant way to do this, but
      % you only need to regrid the DMS oceanic content once per month so it
      % does not matter if this is a bit slow
      is_filling = 1;
      while(is_filling)
        data_wrf_DMS_2 = data_wrf_DMS;
        data_wrf_DMS_2(data_wrf_DMS==0) = NaN;
        is_filling = 0;
        for ii = 1:imax
          for jj = 1:jmax
            if(isnan(data_wrf_DMS(ii, jj)))
              i1_average = ii - 1;
              i2_average = ii+1;
              j1_average = jj - 1;
              j2_average = jj+1;
              if(ii==1)
                i1_average = 1;
              end
              if(ii==imax)
                i2_average = imax;
              end
              if(jj==1)
                j1_average = 1;
              end
              if(jj==jmax)
                j2_average = jmax;
              end
              if(any(any(~isnan(data_wrf_DMS_2(i1_average:i2_average, j1_average:j2_average)))))
                data_wrf_DMS(ii, jj) = nan_mean(nan_mean((data_wrf_DMS_2(i1_average:i2_average, j1_average:j2_average))));
                is_filling = 1;
              end
            end
          end
        end
      end
      data_wrf_DMS(isnan(data_wrf_DMS)) = 0;
      

      % -------- Regrid monthly VOCs --------
      disp(['  Regridding monthly NMVOC'])

      %---- Monthly ECLIPSE VOCs
      disp(['    Regridding ECLIPSE NMVOC'])
      % This is done only once per month (ECLIPSE emissions are monthly), for bulk VOCs
      % The speciation from bulk VOCs to specific WRFChem VOCs is done below and is very fast compared to this step

      %---- Monthly ECLIPSE VOCs: regrid for the closest year before year_ECLIPSE
      ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V6_VOC_' num2str(years_closest_ECLIPSE(1)) '.nc'];
      disp(['    Opening ' ECLIPSE_FILENAME])
      ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
      % ECLIPSE VOCs year_before: Loop over sectors, regrid apply monthly factors
      for isector = 1:nsectors_ECLIPSE
        % Check if emissions exist for this sector
        errorID = 1;
        try
          data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
          errorID = 0;
        catch exception
          errorID = 1;
        end
        if(~errorID)
          data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
          eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
          data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
          data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
          eval(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} '_before = data_wrf_ECLIPSE;']);
        else
          %disp(['      ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
        end % if error
      end % for isector (ECLIPSE emission sectors)
      netcdf.close(ncid); % Close ECLIPSE emission file
      %---- Monthly ECLIPSE VOCs: regrid for the closest year after year_ECLIPSE
      ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V6_VOC_' num2str(years_closest_ECLIPSE(2)) '.nc'];
      disp(['    Opening ' ECLIPSE_FILENAME])
      ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
      % ECLIPSE VOCs year_after: Loop over sectors, regrid apply monthly factors
      for isector = 1:nsectors_ECLIPSE
        % Check if emissions exist for this sector
        errorID = 1;
        try
          data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
          errorID = 0;
        catch exception
          errorID = 1;
        end
        if(~errorID)
          data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
          eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
          data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
          data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
          eval(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} '_after = data_wrf_ECLIPSE;']);
        else
          %disp(['      ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
        end % if error
      end % for isector (ECLIPSE emission sectors)
      netcdf.close(ncid); % Close ECLIPSE emission file
      %----- Monthly ECLIPSE VOCs: linear interpolation emissions to
      % year_ECLIPSE between year_before and year_after
      for isector = 1:nsectors_ECLIPSE
        if(exist(['data_wrf_ECLIPSE_VOC_', ECLIPSE_SECTORS{isector}, '_before'], 'var')...
            & exist(['data_wrf_ECLIPSE_VOC_', ECLIPSE_SECTORS{isector}, '_after'], 'var'))
          eval(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ' = data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ...
                '_after * interp_factor_after_ECLIPSE  + data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ...
                '_before * interp_factor_before_ECLIPSE;']);
        end % if exist (ECLIPSE regridded variables exist)
      end % for isector (ECLIPSE emission sectors)

      
      %-------- Regrid remaining monthly (non-VOC) emissions --------
      disp(['  Regridding remaining ECLIPSE species'])
      for ispecies = 1:nspecies

        %-------- Monthly ECLIPSE non-VOC emissions regridding
        specname_ECLIPSE = SPECNAMES_ECLIPSE{ispecies};
        if(specname_ECLIPSE & ~strcmp(specname_ECLIPSE,'VOC'))
          disp(['    Regridding ECLIPSE emissions for ' specname_ECLIPSE ' : '])
          %---- Monthly ECLIPSE non-VOC: regrid species for the closest year before year_ECLIPSE
          ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V6_' specname_ECLIPSE '_' num2str(years_closest_ECLIPSE(1)) '.nc'];
          disp(['      Opening ' ECLIPSE_FILENAME])
          ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
          % Loop over ECLIPSE sectors
          for isector = 1:nsectors_ECLIPSE
            errorID = 1;
            try
              data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
              errorID = 0;
            catch exception
              errorID = 1;
            end
            if(~errorID)
              data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
              eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
              if(strcmp(specname_ECLIPSE, 'NH3') && strcmp(ECLIPSE_SECTORS{isector}, 'agr'))
                eclipse_mfac = eclipse_mfac_agr_NH3;
              end
              data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
              data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
              eval(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} '_before = data_wrf_ECLIPSE;']);
            else
              % disp(['        ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
            end % if error
          end % for isector (ECLIPSE sectors)
          netcdf.close(ncid); % Close ECLIPSE emission file
          %---- Monthly ECLIPSE non-VOC: regrid species for the closest year after year_ECLIPSE
          ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V6_' specname_ECLIPSE '_' num2str(years_closest_ECLIPSE(2)) '.nc'];
          disp(['      Opening ' ECLIPSE_FILENAME])
          ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
          for isector = 1:nsectors_ECLIPSE
            errorID = 1;
            try
              data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
              errorID = 0;
            catch exception
              errorID = 1;
            end
            if(~errorID)
              data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
              eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
              if(strcmp(specname_ECLIPSE, 'NH3') && strcmp(ECLIPSE_SECTORS{isector}, 'agr'))
                  eclipse_mfac = eclipse_mfac_agr_NH3;
              end
              data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
              data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
              eval(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} '_after = data_wrf_ECLIPSE;']);
            else
              %disp(['        ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
            end % if error
          end % for isector (ECLIPSE emission sectors)
          netcdf.close(ncid); % Close ECLIPSE emission file
          %----- Monthly ECLIPSE non-VOC: Linear time interpolation of species ispecies in ECLIPSE to year_ECLIPSE
          for isector = 1:nsectors_ECLIPSE
            if(exist(['data_wrf_ECLIPSE_' specname_ECLIPSE '_', ECLIPSE_SECTORS{isector}, '_before'], 'var')...
                & exist(['data_wrf_ECLIPSE_' specname_ECLIPSE '_', ECLIPSE_SECTORS{isector}, '_after'], 'var'))
              % Linear interpolation of ECLIPSE emission data to year_ECLIPSE
              eval(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector}...
                    ' = data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} ...
                    '_after * interp_factor_after_ECLIPSE  + data_wrf_ECLIPSE_' ...
                    specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} ...
                    '_before * interp_factor_before_ECLIPSE;']);
            end % if exist (ECLIPSE regridded variables exist)
          end % for isector (ECLIPSE emission sectors)
        end % if(specname_ECLIPSE & ~strcmp(specname_ECLIPSE,'VOC'))
      end % for ispecies
      
      
      %-------- Create wrfchemi, assemble emissions, write to file --------
      days_in_month = days_list_d(months_list_d == month_i & years_list_d == year_i); 
      for iday = 1:length(days_in_month)
        day_i = days_in_month(iday);

        % Hourly loop, 1 emission file per day. Change to ihour = 0:23 to get
        % hourly files, also comment hourly_factors_now = 1.
        for ihour = 0:0
          %---- Set exact time
          % Date in matlab format
          wrfchemi_date = datenum(year_i, month_i, day_i, ihour, 0, 0);
          % Calculate day of week for applying daily emission variation factors
          dayofweek_now = weekday(wrfchemi_date);
          dayofweek_now(dayofweek_now == 0) = 1;
          % Calculate the local "solar" time in order to add hourly variations if needed
          local_time = round(mod(ihour + data_lon_wrf / 180 * 12, 24));
          local_time(local_time == 0) = 24;
          % Date in WRF format
          date_wrf = datestr(wrfchemi_date, 'yyyy-mm-dd_HH:MM:SS');

          %---- Create wrfchemi and open it
          disp(['  Create wrfchemi_', domain, ' for date ', datestr(wrfchemi_date, 'yyyy/mm/dd-HH:MM:SS')])
          wrfchemi_filename = ['wrfchemi_', domain, '_', date_wrf];
          create_netcdffile_saprc_mosaic_dms(wrfchemi_filename, SPECNAMES_SAPRC,...
                                             date_wrf, ZDIM, imax, jmax, ...
                                             dx, dy, ref_lat, ref_lon, truelat1,...
                                             truelat2, moad_cen_lat, map_proj, wrf_proj.mminlu,...
                                             ['Created by Louis Marelle for WRF V4.0 ',...
                                             'using ECLIPSE_v6b emission files']);
          ncid_wrfchemi = netcdf.open(wrfchemi_filename, 'NC_WRITE');


          %-------- Final emission assembly --------- 
          % Emissions speciation, time variation and unit conversion ; writing
          % to wrfchemi
          for ispecies = 1:nspecies
            specname_mechanism = SPECNAMES_SAPRC{ispecies};
            disp(['    Retrieve data and write to wrfchemi, ', specname_mechanism])

            %---- Assemble emissions: ECLIPSE emissions
            specname_ECLIPSE = SPECNAMES_ECLIPSE{ispecies};
            data_wrf_ECLIPSE = zeros(imax, jmax);
            if(specname_ECLIPSE)
              if(strcmp(specname_ECLIPSE, 'VOC'))
                %-- Assemble emissions: ECLIPSE VOCs
                % If species is a VOC, retrieve ECLIPSE bulk VOCs for
                % isector and speciate to WRFChem VOCs
                for isector = 1:nsectors_ECLIPSE
                  if(exist(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector}], 'var'))
                    eval(['data_wrf_ECLIPSE_VOC = data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ';']);
                    % Hour of day emission variation
                    hourly_factors_now = reshape(hourly_factors(isector, local_time), size(local_time));
                    %turn off the hourly factors for daily emissions
                    hourly_factors_now = 1;
                    % Day of week emission variation
                    daily_factor_now = daily_factors(isector, dayofweek_now);
                    % Convert from g(tot VOC)/cell/month to moles(spec VOC)/cell/month
                    data_wrf_ECLIPSE_VOC = voc_speciation_saprc(data_wrf_ECLIPSE_VOC * 1E9,...
                                             specname_mechanism, ECLIPSE_SECTORS{isector});
                    data_wrf_ECLIPSE = data_wrf_ECLIPSE + data_wrf_ECLIPSE_VOC...
                                       .* hourly_factors_now * daily_factor_now;
                  end
                end % for isector
              else 
                %-- Assemble emissions: ECLIPSE non-VOCs
                % If species is not a VOC, just retrieve the already regridded
                % emissions for this species and apply the time variation
                % factors
                for isector = 1:nsectors_ECLIPSE
                  if(exist(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector}], 'var'))
                    eval(['data_wrf_ECLIPSE_SPECIES = data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} ';']);
                    % Hour of day emission variation
                    hourly_factors_now = reshape(hourly_factors(isector, local_time), size(local_time));
                    % Turn off the hourly factors for daily emissions
                    hourly_factors_now = 1;
                    % Day of week emission variation
                    daily_factor_now = daily_factors(isector, dayofweek_now);
                    data_wrf_ECLIPSE = data_wrf_ECLIPSE + data_wrf_ECLIPSE_SPECIES .* hourly_factors_now * daily_factor_now;
                  end
                end % for isector
              end % if(strcmp(specname_ECLIPSE, 'VOC'))
              %-- Assemble emissions: ECLIPSE speciation and unit conversion
              cell_area_km2 = dx ./ data_mapfacx_wrf .* dy ./ data_mapfacy_wrf * 1.0E6;
              cell_area_m2 = cell_area_km2 * 1.0E6;
              hours_in_month = (datenum(0, month_i + 1, 0) - datenum(0, month_i, 0)) * 24;
              seconds_in_month = hours_in_month * 3600;
              if(strcmp(specname_ECLIPSE, 'BC'))
                % kT/month to ug/m2/s
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * 1E15 ./ cell_area_m2 / seconds_in_month;
              elseif(strcmp(specname_ECLIPSE, 'OC'))
                % kT/month OC to ug/m2/s POM
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * 1E15 ./ cell_area_m2 / seconds_in_month * oc_pom_factor;
              elseif(strcmp(specname_ECLIPSE, 'VOC'))
                % moles/month to moles/km2/hour
                data_wrf_ECLIPSE = data_wrf_ECLIPSE ./ cell_area_km2 / hours_in_month;
              else
                % kT/month to moles/km2/hour
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * 1E9 ./ cell_area_km2 / hours_in_month / MOLARWEIGHTS_ECLIPSE(ispecies);
              end
              %-- Assemble emissions: ECLIPSE NOx speciation to NO and NO2
              if (strcmp(specname_mechanism, 'NO'))
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * nox_as_no_anthro;
              elseif (strcmp(specname_mechanism, 'NO2'))
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * (1 - nox_as_no_anthro);
              end
            end % if(specname_ECLIPSE)

            %---- Assemble emissions: Write emissions to wfchemi file
            % Sum all the regridded emissions in data_wrf
            data_wrf = data_wrf_ECLIPSE;
            % Special case for DMS oceanic content data 
            if(strcmp(specname_mechanism,'DMS_OC'))
                data_wrf = data_wrf_DMS;
            end
            % Variable name in the wrfchemi file
            wrfchemi_varname = ['E_' specname_mechanism];
            % Write regridded emissions (data_wrf) to the wrfchemi netcdf
            % variable (E_...)
            ncwrite(ncid_wrfchemi, data_wrf, wrfchemi_varname);

          end % for ispecies

          % Close wrfchemi NetCDF
          netcdf.close(ncid_wrfchemi);

        end % for ihour
      end % for iday
    end % for imonth
  end % for idomain

  disp(' ')
  disp(['-------- ', routine_filename, ' - done --------'])
  disp(' ')

end % function



% -------- Additional sub-functions --------


function [] = ncwrite(ncid, data_var, VARNAME)
  % Write data_var to netcdf variable VARNAME
  varid = netcdf.inqVarID(ncid, VARNAME);
  netcdf.putVar(ncid, varid, data_var);
end


function [data_var] = ncread(ncid, VARNAME)
  % Read netcdf variable VARNAME
  varid = netcdf.inqVarID(ncid, VARNAME);
  data_var = netcdf.getVar(ncid, varid);
end


function [mapping_array] = map_emissions(imax, jmax, imax_emiss, jmax_emiss,...
                             emissions_lon, emissions_lat, DLON_emiss, DLAT_emiss,...
                             lon_min, lon_max, lat_min, lat_max, ncells, truelat1,...
                             truelat2, hemi, stdlon, ref_lat, ref_lon, ref_x,...
                             ref_y, dx, map_proj)
  % Create a "mapping" linking each point of the WRF grid to points of the grid
  % of the emission inventory
  %
  % To regrid emissions, we cut emission inventory cells into minicells
  % We then build for each inventory a table of correspondance linking each
  % inventory minicell to a WRF cell. This table (or mapping) has
  % the horizontal dimensions of the destination WRF grid (xmax, ymax). each
  % value mapping_array{x, y} contains a list
  % (i_cell_1, j_cell_1, i_cell_2, j_cell_2...) where i_cell_n, j_cell_n are
  % the coordinates of the emission inventory cells where the minicell n has
  % been cut

  mapping_array = cell(imax, jmax);
  for i = 1:imax_emiss
    for j = 1:jmax_emiss
      if(emissions_lat(j) > lat_min && emissions_lat(j) < lat_max && emissions_lon(i) > lon_min && emissions_lon(i) < lon_max)
        lat_cells = repmat(linspace(emissions_lat(j) - (DLAT_emiss - DLAT_emiss / ncells) / 2,...
          emissions_lat(j) + (DLAT_emiss - DLAT_emiss / ncells) / 2, ncells), 1, ncells);
        lon_cells = reshape(repmat(linspace(emissions_lon(i) - (DLON_emiss - DLON_emiss / ncells) / 2,...
          emissions_lon(i) + (DLON_emiss - DLON_emiss / ncells) / 2, ncells), ncells, 1), 1, ncells * ncells);
        [x_cells, y_cells] = llij(lat_cells, lon_cells, truelat1, truelat2, hemi, stdlon, ref_lat, ref_lon, ref_x, ref_y, dx, map_proj);
        x_cells = round(x_cells);
        y_cells = round(y_cells);
        for c = 1:ncells * ncells
          if(x_cells(c) > 0 && x_cells(c) < imax + 1 && y_cells(c) > 0 && y_cells(c) < jmax + 1)
            mapping_array{x_cells(c), y_cells(c)} = [mapping_array{x_cells(c), y_cells(c)}, i, j];
          end
        end
      end
    end
  end
end


function [data_wrf_emissions] = regrid_area(data_emissions, mapping_array)
  % Regrid emissions for an inventory with emissions in mass or moles /surface/time

  [imax, jmax] = size(mapping_array);
  data_wrf_emissions = zeros(imax, jmax);
  if(~isempty(data_emissions))
    data_wrf_n_emissions = zeros(imax, jmax);
    for i = 1:imax
      for j = 1:jmax
        if(~isempty(mapping_array{i, j}))
          a = mapping_array{i, j};
          for w = 1:2:length(a)
            % Surface emission data, summed (this sum has no physical meaning without the following average)
            data_wrf_emissions(i, j) = data_wrf_emissions(i, j) + data_emissions(a(w), a(w + 1));
            % Counter for each cell, used to calculate average surface emissions
            data_wrf_n_emissions(i, j) = data_wrf_n_emissions(i, j) + 1;
          end
        end
      end
    end
    % Divide emissions by the counter to calculate average emissions for each wrf grid
    data_wrf_n_emissions(data_wrf_n_emissions == 0) = 1;
    data_wrf_emissions = data_wrf_emissions ./ data_wrf_n_emissions;
  end
end


function [data_wrf_emissions] = regrid(data_emissions, mapping_array, ncells)
  % Regrid emissions for the inventoires with emissions in mass/time

  [imax, jmax] = size(mapping_array);
  data_wrf_emissions = zeros(imax, jmax);
  if(~isempty(data_emissions))
    for i = 1:imax
      for j = 1:jmax
        if(~isempty(mapping_array{i, j}))
          a = mapping_array{i, j};
          for w = 1:2:length(a)
            data_wrf_emissions(i, j) = data_wrf_emissions(i, j) + data_emissions(a(w), a(w + 1)) / (ncells * ncells);
          end
        end
      end
    end
  end
end


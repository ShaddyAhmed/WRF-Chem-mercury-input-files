function [i,j]=llij(lat,lon,truelat1, truelat2, hemi, stdlon,lat1, lon1, knowni, knownj, dx, map_proj)

%-------- Convert latitude and longitude to WRF (i,j) coordinates --------

rad_per_deg =  pi/180;
deg_per_rad =  180/pi;
% Earth radius in kilometers divided by dx
rebydx = 6370./dx;

if(map_proj == 1)
%-------- Lambert conformal projection --------
% Subroutines set_lc, lc_cone and llij_lc from WPS: geogrid/src/module_map_utils.f90

  %---- Subroutine to compute the cone factor of a Lambert Conformal projection

  % Input Args
  %REAL, INTENT(IN)             :: truelat1  % (-90 -> 90 degrees N)
  %REAL, INTENT(IN)             :: truelat2  %   "   "  "   "     "
  % Output Args
  %REAL, INTENT(OUT)            :: cone

  % First, see if this is a secant or tangent projection.  For tangent
  % projections, truelat1 = truelat2 and the cone is tangent to the 
  % Earth's surface at this latitude.  For secant projections, the cone
  % intersects the Earth's surface at each of the distinctly different
  % latitudes
  if(abs(truelat1 - truelat2) > 0.1)
     cone = log10(cos(truelat1 * rad_per_deg)) - ...
            log10(cos(truelat2 * rad_per_deg));
     cone = cone / (log10(tan((45.0 - abs(truelat1) / 2.0) * rad_per_deg)) - ...
            log10(tan((45.0 - abs(truelat2) / 2.0) * rad_per_deg)));
  else
     cone = sin(abs(truelat1)*rad_per_deg );
  end

  %---- Initialize the remaining items in the proj structure for a
  %---- lambert conformal grid.

  % Compute longitude differences and ensure we stay out of the
  % forbidden "cut zone"
  deltalon1 = lon1 - stdlon;
  if(deltalon1 > +180.) 
    deltalon1 = deltalon1 - 360.;
  end
  if(deltalon1 < -180.) 
    deltalon1 = deltalon1 + 360.;
  end

  % Convert truelat1 to radian and compute COS for later use
  tl1r = truelat1 * rad_per_deg;
  ctl1r = cos(tl1r);

  % compute the radius to our known lower-left (sw) corner
  rsw = rebydx * ctl1r/cone * ...
         (tan((90.*hemi-lat1)*rad_per_deg/2.) / ...
          tan((90.*hemi-truelat1)*rad_per_deg/2.)) ^ cone;

  % find pole point
  arg = cone*(deltalon1*rad_per_deg);
  polei = hemi*knowni - hemi * rsw * sin(arg);
  polej = hemi*knownj + rsw * cos(arg);


  %---- Subroutine to compute the geographical latitude and longitude values
  %---- to the cartesian x/y on a Lambert Conformal projection.

  % Input Args
  %REAL, INTENT(IN)              :: lat      % Latitude (-90->90 deg N)
  %REAL, INTENT(IN)              :: lon      % Longitude (-180->180 E)

  % Output Args                 
  %REAL, INTENT(OUT)             :: i        % Cartesian X coordinate
  %REAL, INTENT(OUT)             :: j        % Cartesian Y coordinate

  %---- Compute deltalon between known longitude and standard lon and ensure
  % it is not in the cut zone
  deltalon = lon - stdlon;
  deltalon(deltalon > 180.) = deltalon(deltalon > 180.) - 360.;
  deltalon(deltalon < -180.) = deltalon(deltalon < -180.) + 360.;

  % Convert truelat1 to radian and compute COS for later use
  tl1r = truelat1 * rad_per_deg;
  ctl1r = cos(tl1r);

  % Radius to desired point
  rm = rebydx * ctl1r/cone * ...
       (tan((90.*hemi-lat)*rad_per_deg/2.) / ...
        tan((90.*hemi-truelat1)*rad_per_deg/2.)) .^ cone;

  arg = cone*(deltalon*rad_per_deg);
  i = polei + hemi * rm .* sin(arg);
  j = polej - rm .* cos(arg);

  % Finally, if we are in the southern hemisphere, flip the i/j
  % values to a coordinate system where (1,1) is the SW corner
  % (what we assume) which is different than the original NCEP
  % algorithms which used the NE corner as the origin in the 
  % southern hemisphere (left-hand vs. right-hand coordinate?)
  i = hemi * i;
  j = hemi * j;


elseif(map_proj == 2)
%-------- Polar stereographic projection --------

  %deltalon1 = xlon1 - stdlon
  %if (deltalon1 >180.) 
  %deltalon1 = deltalon1 - 360.;
  %if(deltalon1 < -180.) 
  %deltalon1 = deltalon1 + 360.;

  %Compute the reference longitude by rotating 90 degrees to the east
  %to find the longitude line parallel to the positive x-axis.
  reflon = stdlon + 90;
  %Compute numerator term of map scale factor
  scale_top = 1. + hemi .* sin(truelat1 .* rad_per_deg);
  %Compute radius to lower-left (SW) corner
  ala1 = lat1 .* rad_per_deg;
  rsw = rebydx.*cos(ala1).*scale_top./(1.+hemi.*sin(ala1));
  %Find the pole point
  alo1 = (lon1 - reflon) .* rad_per_deg;
  polei = knowni - rsw .* cos(alo1);
  polej = knownj - hemi .* rsw .* sin(alo1);
  %Find radius to desired point
  ala = lat .* rad_per_deg;
  rm = rebydx .* cos(ala) .* scale_top./(1. + hemi .*sin(ala));
  alo = (lon - reflon) .* rad_per_deg;
  i = polei + rm .* cos(alo);
  j = polej + hemi .* rm .* sin(alo);


elseif(map_proj == 6)
%-------- Rotated lat-lon --------


else
%-------- Other projections are not yet supported --------
  error(['WRF map projection map_proj not recognized'])
end

end

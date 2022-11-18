function [data_VOC_speciated] = voc_speciation_saprc(data_VOC, SPECIES_NAME, SECTOR_NAME)
  % Do the speciation of NMVOC according to UK Emissions of Air Pollutants 1970 to 2007, UK Emissions Inventory Team, AEA 2010 & Zaveri et al, 1999 JGR
  % mapping to SAPRC based on SpecDB.xls from "Development of an Improved Chemical Speciation Database for Processing Emissions of Volatile Organic Compounds for Air Quality Models" By William P. L. Carter 
  % data_VOC in units of mass, output data_VOC_speciated in units of mol
  
  molar_weights=[46.1 , 58.1 , 30.1 , 32 , 44.1 , 92.1 , 28.1 , 58.1 , 72.2, ...
      72.2 , 106.2 , 86.2 , 78.1 , 30 , 131.4 , 58.12 , 72.1 , 84.9 , 142.3, ...
      116.2 , 42.1 , 120.2 , 106.2 , 60.1 , 88.1 , 100.2 , 100.2 , 114.2 , 106.2, ...
      106.2 , 165.8 , 128.3 , 156.3 , 74.1 , 56.1 , 26 , 44.1 , 60.1 , 118.2 , 86.2, ...
      136.2 , 120.2 , 74.1 , 90.1 , 120.2 , 120.2 , 156.3 , 54.1 , 86.2 , 132.2];
  
  % Energy production (tonnes),Commercial/residential combustion (tonnes),Industrial combustion (tonnes),Production processes (tonnes),Extraction/ distribution fossil fuels (tonnes),Solvents (tonnes),Road transport (tonnes),Other transport (tonnes),Waste treatment and disposal (tonnes),NAME,CBMZ Name,Molar mass (g/mol)
  table=[0,6552 , 54 , 55044 , 0 , 41527 , 0 , 0 , 630; ...  , ethanol , ALK3
  186 , 1314 , 291 , 2936 , 41243 , 19314 , 4667 , 228 , 53; ...  , butane , ALK3
  242 , 3778 , 95 , 1195 , 29762 , 0 , 1282 , 221 , 5956; ...  , ethane , ALK1, C2H6
  0 , 0 , 0 , 1486 , 0 , 28719 , 0 , 0 , 163; ...  , methanol , MEOH
  140 , 1530 , 139 , 1594 , 20396 , 3777 , 430 , 152 , 5642; ... , propane , ALK2, C3H8
  94 , 689 , 89 , 3392 , 165 , 10766 , 5431 , 1792 , 364; ...  , toluene , ARO1/ARO1NBZ
  49 , 7186 , 107 , 3570 , 24 , 0 , 6107 , 3734 , 1013; ...  , ethylene , ETHE
  18 , 6 , 15 , 1560 , 0 , 18078 , 578 , 48 , 3; ... , acetone , ACET
  125 , 713 , 320 , 1593 , 15446 , 434 , 2866 , 135 , 46; ...  , pentane , ALK4
  57 , 1241 , 112 , 929 , 7613 , 48 , 6095 , 318 , 34; ...  , 2-methylbutane , ALK4
  460 , 102 , 41 , 1916 , 58 , 12319 , 1658 , 540 , 167; ...  , m-xylene , ARO2 
  124 , 105 , 51 , 3226 , 8004 , 2612 , 2613 , 84 , 240; ...  , hexane , ALK4
  88 , 8213 , 339 , 1584 , 541 , 0.043 , 2218 , 2853 , 985; ...  , benzene , 0.295 ARO1/ BENZENE
  2002 , 673 , 527 , 292 , 39 , 24 , 3815 , 2511 , 3761; ...  , formaldehyde , HCHO
  0 , 0 , 0 , 575 , 0 , 12335 , 0 , 0 , 132; ...  , trichloroethene, ALK3
  17 , 447 , 12 , 210 , 9084 , 973 , 1998 , 120 , 17; ...  , 2-methylpropane , ALK3
  0 , 0 , 0 , 636 , 0 , 11677 , 176 , 8 , 30; ...  , 2-butanone, MEK
  0 , 0 , 0 , 2086 , 0 , 9983 , 0 , 0 , 148; ...  , dichloromethane, ALK1 
  0.107 , 16 , 0 , 681 , 20 , 8284 , 562 , 1345 , 0; ...  , decane , ALK5
  0 , 0 , 0 , 193 , 0 , 10132 , 0 , 0 , 47; ...  , butyl acetate, ALK4
  59 , 1322 , 29 , 3690 , 14 , 0.009 , 2547 , 1092 , 58; ...  , propylene , OLE1, C3H6
  0.107 , 0.003 , 0 , 464 , 4 , 5408 , 1906 , 453 , 0; ...  , 1 , 2 , 4-trimethylbenzene, ARO2
  134 , 35 , 25 , 1538 , 17 , 4717 , 1278 , 348 , 270; ...  , ethylbenzene , ARO1/ARO1NBZ
  0 , 4 , 0 , 561 , 0 , 7717 , 0 , 0 , 36; ...  , 2-propanol, ALK4
  0 , 0 , 0 , 1166 , 0 , 7010 , 0 , 0 , 50; ...  , ethyl acetate, ALK2
  18 , 281 , 1 , 284 , 7857 , 1472 , 620 , 102 , 0; ...  , heptane ,  ALK4
  0 , 0 , 0 , 673 , 0 , 5806 , 0 , 0 , 0; ...  , 4-methyl-2-pentanone , PROD2
  0.214 , 27 , 0 , 182 , 6907 , 1277 , 274 , 37 , 0; ...  , octane , ALK5
  3 , 79 , 18 , 819 , 12 , 3314 , 1282 , 417 , 130; ...  , p-xylene , ARO2
  102 , 58 , 13 , 660 , 27 , 3091 , 1432 , 468 , 95; ...  , o-xylene , ARO2
  0 , 0 , 0 , 119 , 0 , 5661 , 0 , 0 , 272; ...  , tetrachloroethene, ALK1
  0.107 , 23 , 0 , 425 , 51 , 4994 , 140 , 338 , 0; ...  , nonane , ALK5
  0.107 , 0.003 , 0 , 354 , 0 , 4317 , 0 , 639 , 0; ... , undecane , ALK5
  0 , 0 , 0 , 213 , 0 , 4163 , 0 , 0 , 15; ... , 1-butanol, ALK5
  0 , 59 , 0 , 575 , 178 , 0 , 1718 , 1344 , 11; ... , 2-methylpropene , OLE2
  15 , 7 , 57 , 625 , 11 , 0.177 , 2109 , 876 , 0; ... , acetylene , ALK2, C2H2
  0.321 , 0.009 , 0 , 681 , 0 , 0 , 1940 , 1281 , 0; ... , acetaldehyde , CCHO 
  0 , 0 , 0 , 60 , 0 , 3535 , 0 , 0 , 86; ... , 1-propanol, ALK4
  0 , 0 , 0 , 96 , 0 , 3360 , 0 , 0 , 0; ... , 2-butoxyethanol, ALK5
  3 , 5 , 7 , 876 , 1389 , 1236 , 0 , 6 , 114; ... , 2-methylpentane , ALK4
  0 , 0 , 0 , 13 , 0 , 3364 , 0 , 0 , 0; ... , dipentene , TRP1/TERP
  0.214 , 0.006 , 0 , 164 , 0 , 1873 , 721 , 255 , 0; ... , 1,3,5-trimethylbenzene, ARO2
  0 , 0 , 0 , 3061 , 0 , 0 , 0 , 0 , 0; ... , methyl acetate, ALK2
  0 , 0 , 0 , 93 , 0 , 2879 , 0 , 0 , 0; ... , 1-methoxy-2-propanol, ALK5
  0 , 0 , 0 , 225 , 0 , 2657 , 0 , 0 , 0; ... , methylethylbenzene , ARO1NBZ/ARO1
  0.107 , 0.003 , 0 , 155 , 0 , 1874 , 438 , 215 , 0; ... , 1 , 2 , 3-trimethylbenzene , ARO2
  0 , 0 , 0 , 206 , 0 , 2514 , 0 , 0 , 0; ... , 4-methyldecane , ALK5
  1 , 0 , 0 , 402 , 5 , 0 , 1267 , 693 , 16; ... , 1 , 3-butadiene , OLE2
  2 , 3 , 5 , 599 , 755 , 964 , 0 , 0 , 77; ... , 3-methylpentane , ALK4
  0 , 0 , 0 , 48 , 0 , 2216 , 0 , 0 , 0]; % , 1-methoxy-2-propyl acetate, ALK5
  
  
  % Mass ratio for each species
  if(strcmp(SECTOR_NAME,'ene') || strcmp(SECTOR_NAME,'fef') || strcmp(SECTOR_NAME,'energy') || strcmp(SECTOR_NAME,'flr') )
      table_mr=(table(:,1)+table(:,5))/(sum(table(:,1)+table(:,5))); %g/gtotal
  elseif(strcmp(SECTOR_NAME,'ind') || strcmp(SECTOR_NAME,'industry'))
      table_mr=(table(:,3)+table(:,4))/(sum(table(:,3)+table(:,4))); %g/gtotal
  elseif(strcmp(SECTOR_NAME,'dom') || strcmp(SECTOR_NAME,'res') || strcmp(SECTOR_NAME,'residental'))
      table_mr=(table(:,2))/(sum(table(:,2))); %g/gtotal
  elseif(strcmp(SECTOR_NAME,'tra') || strcmp(SECTOR_NAME,'tro') || strcmp(SECTOR_NAME,'tnr') || strcmp(SECTOR_NAME,'transportation'))
      table_mr=(table(:,7))/(sum(table(:,7))); %g/gtotal
  elseif(strcmp(SECTOR_NAME,'agr') || strcmp(SECTOR_NAME,'ags') || strcmp(SECTOR_NAME,'agl') || strcmp(SECTOR_NAME,'agriculture'))
      table_mr=table(:,1)*0; %g/gtotal, no voc emissions from agr
  elseif(strcmp(SECTOR_NAME,'wst') || strcmp(SECTOR_NAME,'swd') || strcmp(SECTOR_NAME,'waste'))
      table_mr=(table(:,9))/(sum(table(:,9))); %g/gtotal
  elseif(strcmp(SECTOR_NAME,'slv') || strcmp(SECTOR_NAME,'solvents'))
      table_mr=(table(:,6))/(sum(table(:,6))); %g/gtotal
  elseif(strcmp(SECTOR_NAME,'shp') || strcmp(SECTOR_NAME,'shipping'))
      table_mr=(table(:,8))/(sum(table(:,8))); %g/gtotal
  else
      disp(['SECTOR_NAME unknown: ',SECTOR_NAME])
  end
  table_factors=table_mr./molar_weights'; %mol/gtotal
  
  % Sections for each SAPRC species
  C2H6_sect = table_factors(3) + table_factors(18) + table_factors(31);
  C3H8_sect = table_factors(5) + table_factors(25) + table_factors(43);
  C2H2_sect = table_factors(36);
  ALK3_sect = table_factors(1) + table_factors(2) + table_factors(15) + table_factors(16);
  ALK4_sect = table_factors(9) + table_factors(10) + table_factors(12) + ...
              table_factors(20) + table_factors(24) + table_factors(26) + table_factors(38) + ...
              table_factors(40) + table_factors(49);
  ALK5_sect = table_factors(19) + table_factors(28) + table_factors(32) + ...
              table_factors(33) + table_factors(34) + table_factors(39) + table_factors(44) + ...
              table_factors(47) + table_factors(50);
  ETHENE_sect = table_factors(7);
  C3H6_sect = table_factors(21);
  % OLE1_sect = 0;
  OLE2_sect = table_factors(35) + table_factors(48);
  ARO1_sect = table_factors(6) + 0.295*table_factors(13) + table_factors(23) +...
              table_factors(45);
  ARO2_sect = table_factors(11) + table_factors(22) + table_factors(29) +...
              table_factors(30) + table_factors(42) + table_factors(46);
  HCHO_sect = table_factors(14);
  CCHO_sect = table_factors(37);
  % RCHO_sect = 0;
  ACET_sect = table_factors(8);
  MEK_sect = table_factors(17);
  % ISO_sect = 0;
  TERP_sect = table_factors(41);
  % SESQ_sect = 0;
  % PHEN_sect = 0;
  % CRES_sect = 0;
  MEOH_sect = table_factors(4);
  % GLY_sect = 0;
  % MGLY_sect = 0;
  % BACL_sect = 0;
  % ISOPROD_sect = 0;
  % METACHRO_sect = 0;
  % MVK_sect = 0;
  PROD_2_sect = table_factors(27);
  % BALD_sect = 0;
  % HCOOH_sect = 0;
  % CCO_OH_setc = 0;
  % RCO_OH_sect = 0;
  
  SPECIES_VOCS_SAPRC = {'C2H6', 'C3H8', 'C2H2', 'ALK3', 'ALK4',...
                        'ALK5', 'ETHENE', 'C3H6', 'OLE2', 'ARO1',...
                        'ARO2', 'HCHO', 'CCHO', 'ACET', 'MEK',...
                        'TERP', 'MEOH', 'PROD2'};
  
  table_cbmz_factors = [C2H6_sect, C3H8_sect, C2H2_sect, ALK3_sect, ALK4_sect,...
                        ALK5_sect, ETHENE_sect, C3H6_sect, OLE2_sect, ARO1_sect,...
                        ARO2_sect, HCHO_sect, CCHO_sect, ACET_sect, MEK_sect,...
                        TERP_sect, MEOH_sect, PROD_2_sect];...
  
  % inventory should be in g, if it is in kg the output will be in kmol
  if(any(strcmp(SPECIES_VOCS_SAPRC, SPECIES_NAME)))
      speciation_factor=table_cbmz_factors(strcmp(SPECIES_VOCS_SAPRC, SPECIES_NAME));
  else
      speciation_factor=0;
  end
  
  data_VOC_speciated=data_VOC*speciation_factor; %mol/cell or kmol/cell

end

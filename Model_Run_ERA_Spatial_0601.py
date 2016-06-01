# -*- coding: utf-8 -*-
"""

Kudryavtsev Model

     By Kang Wang 
     June 1, 2016

Modified from the original version on March 29, 2016  

Note: Performance improved by matrix calculation.
          
Input:

    (1) Location:
        input_lat: Latitude
        input_lon: Longitude
        
    (2) Climate : 
        Ta  : Mean annual air temperature (C)
        Aa  : Amplitude of air temperature (C)
        Hsn : Winter-Averaged Snow Depth (m)
        Rsn : Snow Density (kg/m3)
        vwc : Volumetric Water Content (m3 / m3)
        
    (3) Vegetation:
        Hvgf: Height of vegetation in frozen period (m)
        Hvgt: Height of vegetation in thawed period (m)
        Dvf : Thermal diffusivity of vegetation in frozen period (m2 s)
        Dvt : Thermal diffusivity of vegetation in thawed period (m2 s)

Output:
        1) Mean annual temperature on the top of permafrost (C)
        2) Active Layer Thickness (m)
    
References:
    
    Anisimov, O. A., Shiklomanov, N. I., & Nelson, F. E. (1997). 
        Global warming and active-layer thickness: results from transient general circulation models. 
        Global and Planetary Change, 15(3), 61-77.
    Romanovsky, V. E., & Osterkamp, T. E. (1997). 
        Thawing of the active layer on the coastal plain of the Alaskan Arctic. 
        Permafrost and Periglacial processes, 8(1), 1-22.
    Sazonova, T. S., & Romanovsky, V. E. (2003). 
        A model for regional‐scale estimation of temporal and spatial variability of active layer thickness and mean annual ground temperatures. 
        Permafrost and Periglacial Processes, 14(2), 125-139.
    Sturm, M., Holmgren, J., König, M., & Morris, K. (1997). 
        The thermal conductivity of seasonal snow. Journal of Glaciology, 43(143), 26-41.
    Ling, F., & Zhang, T. (2004). 
        A numerical model for surface energy balance and thermal regime of the active layer and permafrost containing unfrozen water. 
        Cold Regions Science and Technology, 38(1), 1-15.
    Wieder, W.R., J. Boehnert, G.B. Bonan, and M. Langseth. (2014). 
        Regridded Harmonized World Soil Database v1.2. Data set. 
        Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA.  http://dx.doi.org/10.3334/ORNLDAAC/1247  . 
        
"""

import function_20160601 as k_model
import numpy as np
#import matplotlib.pyplot as plt
from time import clock
import csv

print '#============================================'
print '#==== Prepare soil texture for each grid ===='
print '#============================================'

start1=clock()

# Import all soil texture database: ===

lonname    = 'lon';
latname    = 'lat';

input_file = 'Parameters/T_CLAY.nc4'
varname    = 'T_CLAY';
[lat_grid, lon_grid, Clay_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

input_file = 'Parameters/T_SAND.nc4'
varname    = 'T_SAND';
[lat_grid, lon_grid, Sand_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

input_file = 'Parameters/T_SILT.nc4'
varname    = 'T_SILT';
[lat_grid, lon_grid, Silt_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

input_file = 'Parameters/T_OC.nc4'
varname    = 'T_OC';
[lat_grid, lon_grid, Peat_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

input_file0 = 'Input/'+str(2015)+'.csv';

# Extract soil texture for each grid: ===

data00 = np.loadtxt(input_file0,skiprows=1,delimiter=',')

lat_list = data00[:,0]
lon_list = data00[:,1]

n_grid = len(lat_list);

p_clay_list = lat_list*0.;
p_sand_list = lat_list*0.;
p_silt_list = lat_list*0.;
p_peat_list = lat_list*0.;

for i in range(n_grid):
    
    input_lat   = lat_list[i]
    input_lon   = lon_list[i]
    
    [p_clay0, p_sand0, p_silt0, p_peat0] = k_model.Extract_Soil_Texture(input_lat, input_lon, 
                         lon_grid, lat_grid, 
                         Clay_percent, Sand_percent, Silt_percent, Peat_percent);

    p_clay_list[i] = p_clay0
    
    p_sand_list[i] = p_sand0

    p_silt_list[i] = p_silt0

    p_peat_list[i] = p_peat0                         

# Import typical thermal parameters of different soil types ===

input_file2 = 'Parameters/Typical_Thermal_Parameters.csv';

[Capacity_Clay, BDensity_Clay, K_t_Clay, K_f_Clay,
 Capacity_Peat, BDensity_Peat, K_t_Peat, K_f_Peat, 
 Capacity_Sand, BDensity_Sand, K_t_Sand, K_f_Sand, 
 Capacity_Silt, BDensity_Silt, K_t_Silt, K_f_Silt] = k_model.Read_typical_thermal_parameters(input_file2);

finish1=clock()
    
print'Prepare soil texture for each grid: {0:0.1f} s'.format((finish1-start1)) 

print '#==================================='
print '#===   START MAJOR CALCULATION   ==='
print '#==================================='
                         
for i in range(2014,2016):
    
    start1=clock()
    
    #============================
    ## === Inputs preparation ===
    #============================
    
    year_label = i;
    
    print'Runing for {0:d}'.format(year_label)
      
    input_file = 'Input/'+str(year_label)+'.csv';    
    
    data0 = np.loadtxt(input_file,skiprows=1,delimiter=',')
    
    lat_list = data0[:,0]
    lon_list = data0[:,1]
    Ta_list  = data0[:,2]
    Aa_list  = data0[:,3]
    Hsn_list = data0[:,4]
    Rsn_list = data0[:,5]
    vwc_list = data0[:,6]
    Hvgf_list = data0[:,7]
    Hvgt_list = data0[:,8]
    Dvf_list  = data0[:,9]
    Dvt_list  = data0[:,10]
    
    #================================
    ## === Parameters Calculation ===
    #================================
    
    # Estimate snow thermal parameters:
    
    [Csn, Ksn] = k_model.Estimate_Snow_Thermal_Parameters(Rsn_list);
    
    # Estimate Latent Heat:
    
    L = k_model.Estimate_Latent_Heat(vwc_list);

    # Estimate soil heat capacity:    
    
    [Ct, Cf] = k_model.Estimate_Soil_Heat_Capacity(p_clay_list, p_sand_list, p_silt_list, p_peat_list, vwc_list, 
                                Capacity_Clay, BDensity_Clay, 
                                Capacity_Peat, BDensity_Peat,
                                Capacity_Sand, BDensity_Sand,
                                Capacity_Silt, BDensity_Silt);
    
    # Estimate soil thermal conductivity
                            
    [Kt, Kf] = k_model.Estimate_Soil_Thermal_Conductivity(p_clay_list, p_sand_list, p_silt_list, p_peat_list, vwc_list,
                    K_t_Clay, K_f_Clay,  
                    K_t_Peat, K_f_Peat,                              
                    K_t_Sand, K_f_Sand,                                 
                    K_t_Silt, K_f_Silt)
     
    # Estimate length of cold and warm season
               
    [tao,tao1,tao2] = k_model.Estimate_Cold_Warm_Season_Length(Ta_list, Aa_list);
    
    #===========================
    ## === Start calculating ===
    #===========================

    [Tvg, Avg] = k_model.Estimate_Snow_Effect(Ta_list, Aa_list, Hsn_list, Csn, Ksn, Rsn_list, tao);
    
    [Tgs, Ags] = k_model.Estimate_Vegetation_Effect(Tvg, Avg, Hvgf_list, Hvgt_list, 
                                            Dvf_list, Dvt_list, tao1, tao2, tao);
                                            
    Tps_numerator = k_model.Calculate_Tps_numerator(Tgs, Ags, Kf, Kt);
            
    Tps = k_model.Calculate_Temperature_Top_PF(Tps_numerator, Kf, Kt, Cf, Ct);
    
    ALT_grid = k_model.Calculate_ALT(Tps_numerator, Tps, Ags, Kf, Kt, Cf, Ct, L, tao);
       
    finish1=clock()
    
    print'Calculating: {0:0.1f} s'.format((finish1-start1))
    
    #=====================
    ## === Writing out ===
    #=====================
    
    final = np.vstack((np.vstack((lat_list,lon_list)), ALT_grid));
    final = np.transpose(final);
 
    csvfile = file('Output/'+str(year_label)+'_Results_Parallel.csv', 'wb')
    writer = csv.writer(csvfile)
   
    writer.writerows(final)
   
    csvfile.close()
    
print '#=================='
print '#===   FINISH   ==='
print '#=================='
def Estimate_Snow_Thermal_Parameters(rho_sn):
    
    """ 
    The function is to estimate snow thermal conductivity and heat capacity ;
    INPUTs:
            rho_sn: snow density (kg m-3). 

    OUTPUTs:
            Csn: heat capacity of snow  (J kg-1 C-1)
            Ksn: thermal conductivity of snow  (W m-1 C-1)

    DEPENDENTs:
            None 
    """
    
    # Conductivity of snow: 
    #   eq-4, Sturm et al., 1997:    
    Ksn = (rho_sn/1000.)*(rho_sn/1000.)*3.233-1.01*(rho_sn/1000.)+0.138; # Unit: (W m-1 C-1)
    
    # Capacity of snow:
    #   eq-30, Ling et al., 2004; OR Table-1, Goodrich, 1982.
    Csn = 2090.0 + Ksn*0.;                                                # Unit: J m-3 C-1 
    
    return Csn, Ksn

def Read_typical_thermal_parameters(input_file):
    
    """ 
    The function is to import typical thermal parameters of different soil types ;
    
    INPUTs:
            input_file: file contains typical thermal parameters of clay, sand, silt and peat.

    OUTPUTs:
            Capacity_Clay: heat capacity of clay (J kg-1 C-1) 
            BDensity_Clay: Bulk density of clay (kg m-3)
            K_t_Clay     : thermal conductivity of thawed clay (W m-1 C-1)
            K_f_Clay     : thermal conductivity of frozen clay (W m-1 C-1)
            
            Capacity_Silt: heat capacity of silt (J kg-1 C-1) 
            BDensity_Silt: Bulk density of silt (kg m-3)
            K_t_Silt     : thermal conductivity of thawed silt (W m-1 C-1)
            K_f_Silt     : thermal conductivity of frozen silt (W m-1 C-1)
            
            Capacity_Sand: heat capacity of sand (J kg-1 C-1) 
            BDensity_Sand: Bulk density of sand (kg m-3)
            K_t_Sand     : thermal conductivity of thawed sand (W m-1 C-1)
            K_f_Sand     : thermal conductivity of frozen sand (W m-1 C-1)
            
            Capacity_Peat: heat capacity of peat (J kg-1 C-1) 
            BDensity_Peat: Bulk density of peat (kg m-3)
            K_t_Peat     : thermal conductivity of thawed peat (W m-1 C-1)
            K_f_Peat     : thermal conductivity of frozen peat (W m-1 C-1)

    DEPENDENTs:
            None 
    """ 

    import numpy as np    
    
    s_data = np.genfromtxt(input_file, names = True, 
                      delimiter=',', dtype=None)
                    
    Bulk_Density_Texture = s_data['Bulk_Density'];
    Heat_Capacity_Texture = s_data['Heat_Capacity'];
    Thermal_Conductivity_Thawed_Texture = s_data['Thermal_Conductivity_Thawed']
    Thermal_Conductivity_Frozen_Texture = s_data['Thermal_Conductivity_Frozen']
    
    # Clay:
    row = 2;
    Capacity_Clay = Heat_Capacity_Texture[row]
    BDensity_Clay  = Bulk_Density_Texture[row]
    K_t_Clay = Thermal_Conductivity_Thawed_Texture[row]
    K_f_Clay = Thermal_Conductivity_Frozen_Texture[row]
    
    # Sand:
    row = 1;
    Capacity_Sand = Heat_Capacity_Texture[row]
    BDensity_Sand  = Bulk_Density_Texture[row]
    K_t_Sand = Thermal_Conductivity_Thawed_Texture[row]
    K_f_Sand = Thermal_Conductivity_Frozen_Texture[row]    
    
    # Silt:
    row = 0;
    Capacity_Silt = Heat_Capacity_Texture[row]
    BDensity_Silt  = Bulk_Density_Texture[row]
    K_t_Silt = Thermal_Conductivity_Thawed_Texture[row]
    K_f_Silt = Thermal_Conductivity_Frozen_Texture[row]
    
    # Peat:
    row = 3;
    Capacity_Peat = Heat_Capacity_Texture[row]
    BDensity_Peat  = Bulk_Density_Texture[row]
    K_t_Peat = Thermal_Conductivity_Thawed_Texture[row]
    K_f_Peat = Thermal_Conductivity_Frozen_Texture[row]
    
    return Capacity_Clay, BDensity_Clay, K_t_Clay, K_f_Clay, Capacity_Peat, BDensity_Peat, K_t_Peat, K_f_Peat,Capacity_Sand, BDensity_Sand, K_t_Sand, K_f_Sand,Capacity_Silt, BDensity_Silt, K_t_Silt, K_f_Silt

def Estimate_Soil_Heat_Capacity(p_clay, p_sand, p_silt, p_peat, vwc, 
                                Capacity_Clay, BDensity_Clay, 
                                Capacity_Peat, BDensity_Peat,
                                Capacity_Sand, BDensity_Sand,
                                Capacity_Silt, BDensity_Silt):
    
    """ 
    The function is to estimate soil heat capacity ;
    INPUTs:
            input_file: file contains typical thermal parameters of clay, sand, silt and peat.
            p_clay: percent of clay (%)
            p_sand: percent of sand (%)
            p_silt: percent of silt (%)
            p_peat: percent of peat (%)
            vwc   : volumetric water content (m3 / m3)
            
            Capacity_Clay: heat capacity of clay (J kg-1 C-1) 
            BDensity_Clay: Bulk density of clay (kg m-3)
            
            Capacity_Silt: heat capacity of silt (J kg-1 C-1) 
            BDensity_Silt: Bulk density of silt (kg m-3)
            
            Capacity_Sand: heat capacity of sand (J kg-1 C-1) 
            BDensity_Sand: Bulk density of sand (kg m-3)
            
            Capacity_Peat: heat capacity of peat (J kg-1 C-1) 
            BDensity_Peat: Bulk density of peat (kg m-3)

    OUTPUTs:
            Ct: heat capacity of thawed soil (J m-3 C-1) 
            Cf: heat capacity of frozen soil (J m-3 C-1) 

    DEPENDENTs:
            None 
    """    
    
    # Adjust percent of sand, silt, clay and peat ==
    
    tot_percent = p_sand+p_clay+p_silt+p_peat;
    
    percent_sand = p_sand / tot_percent;
    percent_clay = p_clay / tot_percent;
    percent_silt = p_silt / tot_percent;
    percent_peat = p_peat / tot_percent;
    
    # Calculate heat capacity and bulk density of soil using exponential weighted.
                             
    Heat_Capacity =      Capacity_Clay*percent_clay + \
                         Capacity_Sand*percent_sand + \
                         Capacity_Silt*percent_silt + \
                         Capacity_Peat*percent_peat       # Unit: J kg-1 C-1 
                       
    Bulk_Density  =      BDensity_Clay*percent_clay + \
                         BDensity_Sand*percent_sand + \
                         BDensity_Silt*percent_silt + \
                         BDensity_Peat*percent_peat        # Unit: kg m-3
    
    # Estimate heat capacity for composed soil
    # based on the empirical approaches suggested by Anisimov et al. (1997)
        
    Ct = Heat_Capacity*Bulk_Density + 4190.*vwc; # eq-15, Anisimov et al. 1997; Unit: J m-3 C-1
    Cf = Heat_Capacity*Bulk_Density + 2025.*vwc; # eq-15, Anisimov et al. 1997; Unit: J m-3 C-1
    
    return Ct, Cf

def Estimate_Soil_Thermal_Conductivity(p_clay, p_sand, p_silt, p_peat, vwc,
                                       K_t_Clay, K_f_Clay, 
                                       K_t_Peat, K_f_Peat,
                                       K_t_Sand, K_f_Sand,
                                       K_t_Silt, K_f_Silt):

    """ 
    The function is to estimate soil thermal conductivity ;
    INPUTs:
            input_file: file contains typical thermal parameters of clay, sand, silt and peat.
            p_clay: percent of clay (%)
            p_sand: percent of sand (%)
            p_silt: percent of silt (%)
            p_peat: percent of peat (%)
            vwc   : volumetric water content (m3 / m3)
            
            K_t_Clay     : thermal conductivity of thawed clay (W m-1 C-1)
            K_f_Clay     : thermal conductivity of frozen clay (W m-1 C-1)
            
            K_t_Silt     : thermal conductivity of thawed silt (W m-1 C-1)
            K_f_Silt     : thermal conductivity of frozen silt (W m-1 C-1)
            
            K_t_Sand     : thermal conductivity of thawed sand (W m-1 C-1)
            K_f_Sand     : thermal conductivity of frozen sand (W m-1 C-1)
            
            K_t_Peat     : thermal conductivity of thawed peat (W m-1 C-1)
            K_f_Peat     : thermal conductivity of frozen peat (W m-1 C-1)

    OUTPUTs:
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 

    DEPENDENTs:
            None 
    """ 
                      
    # Adjust percent of sand, silt, clay and peat ==
        
    tot_percent = p_sand+p_clay+p_silt+p_peat;
    
    percent_sand = p_sand / tot_percent;
    percent_clay = p_clay / tot_percent;
    percent_silt = p_silt / tot_percent;
    percent_peat = p_peat / tot_percent;
                    
    
    # Estimate thermal conductivity for composed soil 
    
    Kt_Soil =  K_t_Silt**percent_silt * \
               K_t_Clay**percent_clay * \
               K_t_Sand**percent_sand * \
               K_t_Peat**percent_peat
              
    Kf_Soil =  K_f_Silt**percent_silt * \
               K_f_Clay**percent_clay * \
               K_f_Sand**percent_sand * \
               K_f_Peat**percent_peat
    
    # Consider the effect of water content on thermal conductivity
                  
    Kt = Kt_Soil**(1.0-vwc)*0.54**vwc #   Unit: (W m-1 C-1)
    Kf = Kf_Soil**(1.0-vwc)*2.35**vwc #   Unit: (W m-1 C-1)
    
    return Kt, Kf
    
    
def Estimate_Cold_Warm_Season_Length(Ta, Aa):
    
    """ 
    The function is to calculate length of cold season and warm season;
    
    INPUTs:
            Ta: mean annual air temperature (C);
            Aa: amplitude of annual air temperature (C);

    OUTPUTs:
            tao : period of a year (second)
            tao1: period of cold season; (second)
            tao2: period of warm season; (second)
            
    DEPENDENTs:
            None 
    """
    
    import numpy as np
    
    tao = Ta*0.+365.*24.*3600.; # seconds in a year
    
    tao1 = tao*(0.5 - 1./np.pi*np.arcsin(Ta/Aa)); # Cold Season , Page-129, Sazonova, 2003 
    tao2 = tao - tao1;                            # Warm Season , Page-129, Sazonova, 2003
    
    return tao, tao1, tao2

def Estimate_Latent_Heat(vwc):

    """ 
    The function is to calculate latent heat;
    
    INPUTs:
            vwc: volumetric water content (m3 / m3)

    OUTPUTs:
            L : volumetric latent heat of the water  (J/kg)

    DEPENDENTs:
            None 
    """
        
    L = 3.34E8*vwc;                  # Latent Heat, Unit: J kg-1 C-1;
                                     # eq-16, Anisimov et al. 1997
    
    return L
    
def Estimate_Snow_Effect(Ta, Aa, Hsn, Csn, Ksn, rho_sn, tao):

    """ 
    The function is to calculate snow effect;
    
    INPUTs:
            Ta: mean annual air temperature (C)
            Aa: amplitude of air temperature (C)
            Hsn: Snow depth (m)
            Csn: heat capacity of snow  (J m-3 C-1)
            Ksn: thermal conductivity of snow  (W m-1 C-1)
            rho_sn: snow density (kg m-3).             
            tao : period of a year (second)
            
    OUTPUTs:
            Tvg : Temperature on the top of vegetation (C)
            Avg : Amplitude of temperature on the top of vegetation (C)

    DEPENDENTs:
            None 
    """
    
    import numpy as np
    
    #   Estimating Snow Effects
    
    K_diffusivity = Ksn/(rho_sn*Csn)   
    
    deta_Tsn = Aa*(1.0 - np.exp(-1.0*Hsn*np.sqrt(np.pi/(tao*K_diffusivity)))); # eq-7, Anisimov et al. 1997
    deta_Asn = 2.0/np.pi*deta_Tsn; # eq-2, Sazonova et al., 2003
    
    # mean annual temperature and amplitude 
    # bellow snow OR top of vegetation
    
    Tvg = Ta + deta_Tsn;    # Page-129, Sazonova et al., 2003
    Avg = Aa - deta_Asn;    # Page-129, Sazonova et al., 2003
    
    return Tvg, Avg

def Estimate_Vegetation_Effect(Tvg, Avg, Hvgf, Hvgt, Dvf, Dvt, tao1, tao2, tao):
    
    """ 
    The function is to calculate vegetation effect;
    
    INPUTs:
            Tvg : Temperature on the top of vegetation (C)
            Avg : Amplitude of temperature on the top of vegetation (C)
            Hvgf: Height of vegetation in frozen period (m)
            Hvgt: Height of vegetation in thawed period (m)
            Dvf : Thermal diffusivity of vegetation in frozen period (m2 s)
            Dvt : Thermal diffusivity of vegetation in thawed period (m2 s)
            tao : period of a year (second)
            tao1: period of cold season; (second)
            tao2: period of warm season; (second)
            
    OUTPUTs:
            Tgs : Temperature on the ground surface (C)
            Ags : Amplitude of temperature on the ground surface (C)

    DEPENDENTs:
            None 
    """
    
    import numpy as np

    # Estimate vegetation effects
    
    # WINTER vegetation thermal effects:  
    #   eq-10, Anisimov et al. 1997
    deta_A1 = (Avg - Tvg) * \
              (1.-np.exp(-1.*Hvgf*np.sqrt(np.pi/(Dvf*2.*tao1))));
              
    # SUMMER vegetation thermal effects:  
    #   eq-11, Anisimov et al. 1997        
    deta_A2 = (Avg  + Tvg) * \
              (1.-np.exp(-1.*Hvgt*np.sqrt(np.pi/(Dvt*2.*tao2))));
              
    # Effects of vegetation on seasonal amplitude of temperature: 
    #   eq-8, Anisimov et al. 1997
    deta_Av = (deta_A1*tao1+deta_A2*tao2) / tao;
    
    # Effects of vegetation on annual mean temperature: 
    #   eq-9, Anisimov et al. 1997    
    deta_Tv = (deta_A1*tao1-deta_A2*tao2) / tao * (2. / np.pi);
    
    # mean annual temperature and amplitude 
    # on the ground surface
    Tgs = Tvg + deta_Tv;    # eq-13, Sazonova et al., 2003
    Ags = Avg - deta_Av;    # eq-14, Sazonova et al., 2003
    
    return Tgs, Ags
    
def Calculate_Tps_numerator(Tgs, Ags, Kf, Kt):
    
    """ 
    The function is to calculate Tps_Numerator;
    
    INPUTs:
            Tgs : Temperature on the ground surface (C)
            Ags : Amplitude of temperature on the ground surface (C)
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 
            
    OUTPUTs:
            Tps_numerator : Numerator in Equation of
                            Temperature on the top of permafrost

    DEPENDENTs:
            None 
    """    
    
    import numpy as np
        
    Tps_numerator = 0.5*Tgs*(Kf+Kt)\
                    +(Ags*(Kt-Kf)/np.pi\
                    *(Tgs/Ags*np.arcsin(Tgs/Ags)\
                    +np.sqrt(1.-(np.pi**2.0/Ags**2.0)))); # eq-14, Anisimov et al. 1997
                    
    return Tps_numerator

def Calculate_Temperature_Top_PF(Tps_numerator, Kf, Kt, Cf, Ct):
    
    """ 
    The function is to calculate temperature on the top of permafrost;
    
    INPUTs:
            Tps_numerator : Numerator in Equation of
                            Temperature on the top of permafrost
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 
            Ct: heat capacity of thawed soil (J m-3 C-1) 
            Cf: heat capacity of frozen soil (J m-3 C-1) 
            
    OUTPUTs:
            Tps : Temperature on the top of permafrost (C)

    DEPENDENTs:
            None 
    """ 
   
    import numpy as np
   
#   Calculating the temperature at the top of permafrost :
                        
    K_star = Kf;
    K_star[np.where(Tps_numerator>0.0)] = Kt;
               
    # Temperature at the top of permafrost    
       
    Tps = Tps_numerator/K_star;  # eq-14, Anisimov et al. 1997
    
    return Tps
    
def Calculate_ALT(Tps_numerator, Tps, Ags, Kf, Kt, Cf, Ct, L, tao):
    
    """ 
    The function is to calculate active layer thickness;
    
    INPUTs:
            Tps_numerator : Numerator in Equation of
                            Temperature on the top of permafrost
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 
            Ct: heat capacity of thawed soil (J m-3 C-1) 
            Cf: heat capacity of frozen soil (J m-3 C-1) 
            L : volumetric latent heat of the water  (J/kg)
            tao : period of a year (second)
            
    OUTPUTs:
            Zal : Active Layer thickness (m)

    DEPENDENTs:
            None 
    """ 
    
    import numpy as np
            
    C = Cf; C[np.where(Tps_numerator>0.0)] = Ct;
    
    K = Kf; K[np.where(Tps_numerator>0.0)] = Kt;
    
    Aps = (Ags - abs(Tps))/np.log((Ags+L/(2.*C)) / \
          (abs(Tps)+L/(2.*C))) - \
          L/(2.*C);                                                           # eq-4, Romanovsky et al. 1997
    
    Zc = (2.*(Ags - abs(Tps))*np.sqrt((K*tao*C)/np.pi)) / (2.*Aps*C + L);     # eq-5, Romanovsky et al. 1997
    
    Zal = (2.*(Ags - abs(Tps))*np.sqrt(K*tao*C/np.pi)\
        +(((2.*Aps*C*Zc+L*Zc)*L*np.sqrt(K*tao/(np.pi*C)))\
        /(2.*Ags*C*Zc + L*Zc +(2.*Aps*C+L)*np.sqrt(K*tao/(np.pi*C)))))\
        /(2.*Aps*C+L);                                                    # Active Layer Thickness, eq-3, Romanovsky et al. 1997
                
    Zal[np.where(Tps_numerator>0.0)] = -999.99
    Tps[np.where(Tps_numerator>0.0)] = -999.99
    
    return Zal
    
def Extract_Soil_Texture(input_lat, input_lon, 
                         lon_grid, lat_grid, 
                         Clay_percent, Sand_percent, Silt_percent, Peat_percent):
    
    """ 
    The function is to extract the soil texture according to input of latitude and longitude;
    INPUTs:
            input_lat: Latitude;
            input_lon: Longitude;
            
    OUTPUTs:
            p_clay: percent of clay (%)
            p_sand: percent of sand (%)
            p_silt: percent of silt (%)
            p_peat: percent of peat (%)  
            
    DEPENDENTs:
            function "Extract_Grid_Value"
            function "import_ncfile"
    """
    
    clay_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Clay_percent)
                         
    sand_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Sand_percent)
    
    silt_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Silt_percent)
                                   
    peat_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Peat_percent)
                         
    return clay_perc, sand_perc, silt_perc, peat_perc
    
def import_ncfile(input_file, lonname,  latname,  varname):
    
    """ 
    The function is to import a whole matrix from NetCDF file
    
    INPUTs:
            input_file: filename of netcdf file;
            lonname   : name of Longitude in netcdf file;
            latname   : name of Latitude in netcdf file;
            varname   : name of variable in netcdf file;
            
    OUTPUTs:
            lon_grid : Array of longitude
            lat_grid : Array of latitude
            p_data: whole value (Matrix)   
                    
    DEPENDENTs:
            None 
    """
    
    from netCDF4 import Dataset
        
    # Read the nc file 
        
    fh = Dataset(input_file, mode='r')
        
    # Get the lat and lon
        
    lon_grid = fh.variables[lonname][:]; 
    lat_grid = fh.variables[latname][:];
        
    p_data  = fh.variables[varname][:];
        
    return lat_grid,lon_grid,p_data
    
def Extract_Grid_Value(input_lat, input_lon, lon_grid, lat_grid, p_data): 

    """ 
    The function is to extract the grid value from matrix,
    according to input of latitude and longitude;
    
    INPUTs:
            input_lat: Latitude;
            input_lon: Longitude;
            lon_grid : Array of longitude
            lat_grid : Array of latitude
            p_data   : Matrix of data (from NetCDF file)
            
    OUTPUTs:
            q_data: grid value (SINGLE)   
                    
    DEPENDENTs:
            None 
    """
        
    import numpy as np
    
    lon_grid_scale = 0.05;
    lat_grid_scale = 0.05;
    
    lon_grid_top = lon_grid + lon_grid_scale / 2.0;
    lat_grid_top = lat_grid + lat_grid_scale / 2.0;
    
    lon_grid_bot = lon_grid - lon_grid_scale / 2.0;
    lat_grid_bot = lat_grid - lat_grid_scale / 2.0;
    
    # Get the index of input location acccording to lat and lon inputed
    
    idx_lon = np.where((input_lon <= lon_grid_top) & (input_lon >= lon_grid_bot))          
    idx_lat = np.where((input_lat <= lat_grid_top) & (input_lat >= lat_grid_bot))
    
    idx_lon = np.array(idx_lon)
    idx_lat = np.array(idx_lat)
    
    if np.size(idx_lon) >= 1 and np.size(idx_lat) >= 1:
        q_data  = p_data[idx_lat[0,0], idx_lon[0,0]]
    else:
        q_data  = np.nan;
    
    return q_data
    
    
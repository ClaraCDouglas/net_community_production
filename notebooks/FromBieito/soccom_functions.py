import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import netCDF4
import scipy.interpolate as intrp
import datetime
import gsw
import seawater as sw
import os
#from mpl_toolkits.basemap import Basemap
from cartopy import crs
import cartopy.feature as cfeat
import cmocean
#import pygamma
import copy
import glob
import xarray as xr
from holteandtalley import HolteAndTalley
import time

class grids_one_buoy():
    def __init__(self,filename,**kargs):
        if "den_ml_crit" in kargs:
            den_ml_crit = kargs["den_ml_crit"]
        else:
            den_ml_crit = 0.03
        if "DO_ml_crit" in kargs:
            DO_ml_crit = kargs["DO_ml_crit"]
        else:
            #DO_ml_crit = 1. #Kortzinger 2008 proportional to 0.03 kg/m3 if 0.125 kg/m-3 in kortzinger
            #DO_ml_crit = 5. #Kortzinger 2008
            DO_ml_crit = 2.5
        if "dz" in kargs:
            dz = kargs["dz"]
        else:
            dz = 5.
        if "dzLT" in kargs:
            dzLT = kargs["dzLT"]
        else:
            dzLT = 20.
        if "gridding" in kargs:
            gridding = kargs["gridding"]
        else:
            gridding = False
        if "display_info" in kargs:
            display_info = kargs["display_info"]
        else:
            display_info = False
        if "verbose" in kargs:
            verbose = kargs["verbose"]
        else:
            verbose = False
        if "clear_short" in kargs:
            #clears short cut propfiles at 950 m
            clear_short = kargs["clear_short"]
        else:
            clear_short = False

        nfc = netCDF4.Dataset(filename)
        metadata = nfc.__dict__["Comments"]
        if display_info:
            display(nfc)
        variables = list(nfc.variables.keys())
        #print(nfc)
        self.raw = dict()
        self.raw["depth"] = nfc["Depth"][:] #
        self.raw["Lat"] = nfc["Lat"][:]
        self.raw["Lon"] = nfc["Lon"][:]
        self.raw["Lon"][self.raw["Lon"]>180] = self.raw["Lon"][self.raw["Lon"]>180] - 360.
        
        #UOW CODE
        i0 = filename.rfind("/")+1
        i1 = filename.rfind("_")
        self.raw["code"]= filename[i0:i1]
        #WMO code
        WMO_str = "WMO ID:"
        i0 = metadata.find(WMO_str) + len(WMO_str) + 1
        i1 = metadata[i0:].find("\n") + i0
        self.raw["WMO_code"] = metadata[i0:i1]
        
        ref_date_str = nfc["REFERENCE_DATE_TIME"][:].tostring().decode("ascii")
        ref_date = datetime.datetime.strptime(ref_date_str,"%Y%m%d%H%M%S")
        self.raw["date"] = nfc["JULD"][:] + ref_date.toordinal() 
        self.raw["date_dt"] = convert_time_to_date(self.raw["date"])
        #reads the variables
        self.raw["depth"] = nfc["Depth"][:].T # .T is transpose function
        if np.ma.isMaskedArray(self.raw["depth"]):
            self.raw["depth"].mask = (self.raw["depth"].mask) | (nfc["Depth_QFA"][:].T == 8) | (self.raw["depth"]<0)
        else:
            self.raw["depth"] = np.ma.array(self.raw["depth"])
            self.raw["depth"].mask = (nfc["Depth_QFA"][:].T == 8)
            
        self.raw["Pressure"] = nfc["Pressure"][:].T
        if np.ma.isMaskedArray(self.raw["Pressure"]):
            self.raw["Pressure"].mask = (self.raw["Pressure"].mask) | (nfc["Pressure_QFA"][:].T == 8)
        else:
            self.raw["Pressure"] = np.ma.array(self.raw["Pressure"])
            self.raw["Pressure"].mask = (nfc["Pressure_QFA"][:].T == 8)
        
        self.raw["Temperature"] = nfc["Temperature"][:].T
        if np.ma.isMaskedArray(self.raw["Temperature"]):
            self.raw["Temperature"].mask = (self.raw["Temperature"].mask) | (nfc["Temperature_QFA"][:].T == 8)
        else:
            self.raw["Temperature"] = np.ma.array(self.raw["Temperature"])
            self.raw["Temperature"].mask =  (nfc["Temperature_QFA"][:].T == 8)

        self.raw["Salinity"] = nfc["Salinity"][:].T
        if np.ma.isMaskedArray(self.raw["Salinity"]):
            self.raw["Salinity"].mask = (self.raw["Salinity"].mask) | (nfc["Salinity_QFA"][:].T == 8)
        else:
            self.raw["Salinity"] = np.ma.array(self.raw["Salinity"])
            self.raw["Salinity"].mask =  (nfc["Salinity_QFA"][:].T == 8)


        #derived values
        self.raw["SA"] = gsw.SA_from_SP( self.raw["Salinity"], self.raw["Pressure"], self.raw["Lon"], self.raw["Lat"] ) #-10.1325 #absolute salinity from practical salinity
        self.raw["CT"] = gsw.CT_from_t(self.raw["SA"],self.raw["Temperature"],self.raw["Pressure"]) #-10.1325 #Conservative Temperature of seawater from in-situ temperature
        self.raw["Sigma_theta"]  = gsw.sigma0(self.raw["SA"],self.raw["CT"]) #Calculates potential density anomaly with reference pressure of 0 dbar,
#        self.raw["gamma_n"] = np.transpose(pygamma.gamma_n( self.raw["Salinity"].T, self.raw["Temperature"].T, self.raw["Pressure"].T, self.raw["Lon"], self.raw["Lat"]   )[0]) # neutral density calculations removed for the moment because I don't have fortran/coding skills to obtain np package
#        if not np.ma.isMaskedArray(self.raw["gamma_n"]):
#            self.raw["gamma_n"] = np.ma.array( self.raw["gamma_n"] )
        self.raw["gamma_n"] =self.raw["Sigma_theta"] # using potential density instead of neutral density
    
        #biogeochemical
        bg_vars = ["Oxygen","OxygenSat","Nitrate","DIC_LIAR","TALK_LIAR","pCO2_LIAR","Chla_corr","POC"]
        self.raw_bg = dict()
        if "Oxygen" in variables:
            self.raw_bg["Oxygen"] = nfc["Oxygen"][:].T
            if np.ma.isMaskedArray(self.raw_bg["Oxygen"]):
                self.raw_bg["Oxygen"].mask = (self.raw_bg["Oxygen"].mask) | (nfc["Oxygen_QFA"][:].T == 8)
            else:
                self.raw_bg["Oxygen"] = np.ma.array(self.raw_bg["Oxygen"])
                self.raw_bg["Oxygen"].mask =  (nfc["Oxygen_QFA"][:].T == 8)

        if "OxygenSat" in variables:
            self.raw_bg["OxygenSat"] = nfc["OxygenSat"][:].T
            if np.ma.isMaskedArray(self.raw_bg["OxygenSat"]):
                self.raw_bg["OxygenSat"].mask = (self.raw_bg["OxygenSat"].mask) | (nfc["OxygenSat_QFA"][:].T == 8)
            else:
                self.raw_bg["OxygenSat"] = np.ma.array(self.raw_bg["OxygenSat"])
                self.raw_bg["OxygenSat"].mask =  (nfc["OxygenSat_QFA"][:].T == 8)

        if "Nitrate" in variables:
            self.raw_bg["Nitrate"] = nfc["Nitrate"][:].T
            if np.ma.isMaskedArray(self.raw_bg["Nitrate"]):
                self.raw_bg["Nitrate"].mask = (self.raw_bg["Nitrate"].mask) | (nfc["Nitrate_QFA"][:].T == 8)
            else:
                self.raw_bg["Nitrate"] = np.ma.array(self.raw_bg["Nitrate"])
                self.raw_bg["Nitrate"].mask =  (nfc["Nitrate_QFA"][:].T == 8)

        if "DIC_LIAR" in variables:
            self.raw_bg["DIC_LIAR"] = nfc["DIC_LIAR"][:].T
            if np.ma.isMaskedArray(self.raw_bg["DIC_LIAR"]):
                self.raw_bg["DIC_LIAR"].mask = (self.raw_bg["DIC_LIAR"].mask)  | (nfc["DIC_LIAR_QFA"][:].T == 8)
            else:
                self.raw_bg["DIC_LIAR"] = np.ma.array(self.raw_bg["DIC_LIAR"])
                self.raw_bg["DIC_LIAR"].mask =  (nfc["DIC_LIAR_QFA"][:].T == 8)

        if "TALK_LIAR" in variables:
            self.raw_bg["TALK_LIAR"] = nfc["TALK_LIAR"][:].T
            if np.ma.isMaskedArray(self.raw_bg["TALK_LIAR"]):
                self.raw_bg["TALK_LIAR"].mask = (self.raw_bg["TALK_LIAR"].mask) | (nfc["TALK_LIAR_QFA"][:].T == 8)
            else:
                self.raw_bg["TALK_LIAR"] = np.ma.array(self.raw_bg["TALK_LIAR"])
                self.raw_bg["TALK_LIAR"].mask =  (nfc["TALK_LIAR_QFA"][:].T == 8)

        if "pCO2_LIAR" in variables:
            self.raw_bg["pCO2_LIAR"] = nfc["pCO2_LIAR"][:].T
            if np.ma.isMaskedArray(self.raw_bg["pCO2_LIAR"]):
                self.raw_bg["pCO2_LIAR"].mask = (self.raw_bg["pCO2_LIAR"].mask) | (nfc["pCO2_LIAR_QFA"][:].T == 8)
            else:
                self.raw_bg["pCO2_LIAR"] = np.ma.array(self.raw_bg["pCO2_LIAR"])
                self.raw_bg["pCO2_LIAR"].mask =  (nfc["pCO2_LIAR_QFA"][:].T == 8)

        if "Chl_a_corr" in variables:
            self.raw_bg["Chl_a"] = nfc["Chl_a_corr"][:].T
            if np.ma.isMaskedArray(self.raw_bg["Chl_a"]):
                self.raw_bg["Chl_a"].mask = (self.raw_bg["Chl_a"].mask) | (nfc["Chl_a_corr_QFA"][:].T == 8)
            else:
                self.raw_bg["Chl_a"] = np.ma.array(self.raw_bg["Chl_a"])
                self.raw_bg["Chl_a"].mask =  (nfc["Chl_a_corr_QFA"][:].T == 8)

        if "POC" in variables:
            self.raw_bg["POC"] = nfc["POC"][:].T
            if np.ma.isMaskedArray(self.raw_bg["POC"]):
                self.raw_bg["POC"].mask = (self.raw_bg["POC"].mask) | (nfc["POC_QFA"][:].T == 8)
            else:
                self.raw_bg["POC"] = np.ma.array(self.raw_bg["POC"])
                self.raw_bg["POC"].mask =  (nfc["POC_QFA"][:].T == 8)

        
    
        nt = self.raw["Temperature"].shape[1] # shape[1]=profiles

        #LT
        self.raw["LT_ov"] = np.full( self.raw["Temperature"].shape, np.nan )
        self.raw["size_ov"] = np.full( self.raw["Temperature"].shape, np.nan )
        
        #grids - making grids of depthxprofiles, filled with nans for the moment
        self.gr = dict()
        self.gr["depth"] = np.arange(0,2000+dz,dz) # default above is 5 m intervals
        nz = self.gr["depth"].size
        self.gr["date"] = np.copy(self.raw["date"])
        #self.gr["date_dt"] = convert_time_to_date(self.gr["date"])
        self.gr["Lon"] = np.copy(self.raw["Lon"])
        self.gr["Lat"] = np.copy(self.raw["Lat"])
        self.gr["code"] = copy.copy(self.raw["code"])
        self.gr["WMO_code"] = copy.copy(self.raw["WMO_code"])
        #gridded variables
        self.gr["Pressure"] = np.full((nz, nt), np.nan)
        self.gr["Temperature"] = np.full((nz, nt), np.nan)
        self.gr["Salinity"] = np.full((nz, nt), np.nan)
        self.gr["SA"] = np.full((nz, nt), np.nan)
        self.gr["CT"] = np.full((nz, nt), np.nan)
        self.gr["Sigma_theta"] = np.full((nz, nt), np.nan)
        self.gr["gamma_n"] = np.full((nz, nt), np.nan)
        self.gr["N2"] = np.full((nz, nt), np.nan)
        self.gr["PV"] = np.full((nz, nt), np.nan)

        #biogeochemical variables
        for var in bg_vars:
            self.gr[var] = np.full((nz, nt), np.nan)
        
        #mixing parameters
        self.gr["LT"] = np.full((nz, nt), np.nan)
        self.gr["mld"] = np.full(nt, np.nan) # one value per profile
        self.gr["mld_HT"] = np.full(nt, np.nan)
        #self.gr["gpa0"] = np.full(nt, np.nan) # not sure what this would be if it was uncommented?
        self.gr["mld_DO"] = np.full(nt, np.nan) 
        self.gr["LT_ml"] = np.full(nt, 0.)
        self.gr["LT_ov"] = np.full((nz,nt), 0.)
        self.gr["LT_largest_ov"] = np.full(nt, 0.)
        self.gr["size_largest_ov"] = np.full(nt, 0.)
        self.gr["h_largest_ov"] = np.full(nt, 0.)
        self.gr["h_no_ov"] = np.full(nt, 0.)
        for i in range(nt):
            if verbose:
                print("Float %s, profile: %d"%(self.raw["code"],i+1))
            #Interpolates temperature
            ii = np.argsort(self.raw["depth"][:,i]) # argsort creates an array the same size as input. Each cell in array is given the index of where to find the number in input array that should be located in that cell in order to have numerical order. This looks across profile (column) i, down all entries (rows) in column i
            z0 = self.raw["depth"][ii,i] # z0 is the depths sorted from shallowest to deepest, with blank values following
            #deletes profiles shorter than 950 m
            if clear_short and max(z0)<950:
                continue # if chosen to delete short profiles, then loop goes back up to start without loading data into grid, as below
            p0 = self.raw["Pressure"][ii,i] # p0 is the pressure values sorted according to the argsort above
            T0 = self.raw["Temperature"][ii,i] # same, T0 is the temp sorted according to the argsort above (so orders the data for the shallowest to deepest value, then puts all the nans at the end
            msk = ~((T0.mask) | (z0.mask))
            self.gr["Temperature"][:,i] = grids_interpolates(z0[msk], T0[msk], self.gr["depth"], dz, grid = gridding)

            #Pressure
            msk = ~((p0.mask) | (z0.mask))
            self.gr["Pressure"][:,i] = grids_interpolates(z0[msk], p0[msk], self.gr["depth"], dz, grid = gridding)

            #Interpolates potential temperature
            CT0 = self.raw["CT"][ii,i]
            msk = ~((CT0.mask) | (z0.mask))
            self.gr["CT"][:,i] = grids_interpolates(z0[msk], CT0[msk], self.gr["depth"], dz, grid = gridding)

            #Interpolates salinity
            S0 = self.raw["Salinity"][ii,i]
            msk = ~((S0.mask) | (z0.mask))
            self.gr["Salinity"][:,i] = grids_interpolates(z0[msk], S0[msk], self.gr["depth"], dz, grid = gridding)

            #Interpolates SA
            SA0 = self.raw["SA"][ii,i]
            msk = ~((SA0.mask) | (z0.mask))
            self.gr["SA"][:,i] = grids_interpolates(z0[msk], SA0[msk], self.gr["depth"], dz, grid = gridding)
            
            #Interpolates density
            Sigma_theta0 = self.raw["Sigma_theta"][ii,i]
            msk = ~((Sigma_theta0.mask) | (z0.mask))
            self.gr["Sigma_theta"][:,i] = grids_interpolates(z0[msk], Sigma_theta0[msk], self.gr["depth"], dz, grid = gridding)

            #Interpolates gamma_n
            gamma_n0 = self.raw["gamma_n"][ii,i]
            msk = ~((gamma_n0.mask) | (z0.mask))
            self.gr["gamma_n"][:,i] = grids_interpolates(z0[msk].T, gamma_n0[msk].T, self.gr["depth"], dz, grid = gridding)

            ##
            #interpolates the biogeochemical variables
            ##
            for var in bg_vars:
                if var in self.raw_bg.keys():
                    XX = self.raw_bg[var][ii,i]
                    msk = ~((XX.mask) | (z0.mask))
                    if np.nansum(msk)>10:
                        self.gr[var][:,i] = grids_interpolates(z0[msk], XX[msk],self.gr["depth"], dz, grid = gridding)

            #mixed layer depth from density
            msk = ~((Sigma_theta0.mask) | (z0.mask))
            self.gr["mld"][i] = mixed_layer_depth(z0[msk],np.sort(np.array([Sigma_theta0[msk]]).T), Dd = den_ml_crit)[0]
            
            #Mixed layer Holte and Talley
            Pgr = self.gr["Pressure"][:,i]
            CTgr = self.gr["CT"][:,i]
            SAgr = self.gr["SA"][:,i]
            STgr = self.gr["Sigma_theta"][:,i]
            msk = ~( np.isnan(Pgr+CTgr+SAgr+STgr))
            if np.sum(msk)>10:
                html = HolteAndTalley( Pgr[msk], CTgr[msk], SAgr[msk], STgr[msk] )
                self.gr["mld_HT"][i] = html.densityMLD
            
            #stratification
            #N2,pmid = gsw.Nsquared( self.gr["SA"][:,i], self.gr["CT"][:,i], self.gr["Pressure"][:,i]-10.1325  )
            ddendz = first_centered_differences( -self.gr["depth"], self.gr["Sigma_theta"][:,i] )
            self.gr["N2"][:,i] = -(1000+self.gr["Sigma_theta"][:,i])**-1*gsw.grav( self.gr["Lat"][i],self.gr["Pressure"][:,i] )*ddendz #-10.1325
            self.gr["PV"][:,i] = (1000+self.gr["Sigma_theta"][:,i])**-1*gsw.f( self.gr["Lat"][i] )*ddendz
            #self.gr["PV"][:,i] = sw.f( self.gr["Lat"][i] )*self.gr["N2"][:,i]
            """
            #geopotential anomaly
            msk = ~( (S0.mask) | (T0.mask) | (p0.mask) )
            if np.sum(msk)>10:
                self.gr["gpa0"][i] = geopotential_anomaly(CT0[msk],SA0[msk], p0[msk])
            """
            
            #calculates thorpe displacements and mean LT
            igood = np.where( ~((Sigma_theta0.mask) | (z0.mask) ))[0]
            if igood.size<10:
                continue
            Sigma_theta00 = Sigma_theta0[igood].data
            z00 = z0[igood].data
            isort = np.argsort( Sigma_theta00)
            disp = z00 - z00[isort]
            nz1000 = np.where( self.gr["depth"]<=1000 )[0][-1]
            for j in range(nz1000):
                if self.gr["depth"][j]>1000:
                    break
                jj = (z00>= self.gr["depth"][j]-dzLT) & (z00<= self.gr["depth"][j]+dzLT)
                self.gr["LT"][j,i] = np.nanmean(disp[jj]**2)**0.5
    
            #detection of Thorpe overturns
            ii1000 = (z00<=1000) & (np.isfinite(Sigma_theta00))
            zth,LT, ovsize, ovnum = calculates_thorpe_scale(z00[ii1000], Sigma_theta00[ii1000])
            self.raw["LT_ov"][:,i] = grids_interpolates(zth,LT,self.raw["depth"][:,i].data, dz, grid = gridding)
            self.raw["size_ov"][:,i] = grids_interpolates(zth,ovsize,self.raw["depth"][:,i].data,dz)
            self.gr["LT_ov"][:,i] = grids_interpolates(zth,LT,self.gr["depth"], dz, grid = gridding)

                        
            #mean thorpe displacement in the mixed layer
            jjmld = np.where(z00<=self.gr["mld"][i])[0]
            if jjmld.size>0:
                self.gr["LT_ml"][i] = np.nanmean( (disp[jjmld]-np.mean(disp[jjmld]))**2)**0.5
            else:
                self.gr["LT_ml"][i] = 0.
                
            #stores the size and LT of biggest overturn within the mixed layer
            jjml = np.where(zth<=self.gr["mld"][i])[0]
            if jjml.size:
                j_largest = jjml[ np.argmax(ovsize[jjml]) ]
                n_largest_ov = ovnum[ j_largest ]
                j_bot_largest = np.where(ovnum == n_largest_ov)[0][-1]
                if n_largest_ov>0:
                    self.gr["size_largest_ov"][i] = ovsize[0]
                    self.gr["LT_largest_ov"][i] = LT[0]
                    self.gr["h_largest_ov"][i] = zth[ j_bot_largest]

            #first depth with no overturn
            i_nov = np.where(ovsize==0.)[0]
            if i_nov.size>0:
                self.gr["h_no_ov"][i] = zth[ i_nov[0] ]
            else:
                self.gr["h_no_ov"][i] = zth[ -1 ]

            #mixed layer from oxygen
            if "Oxygen" in self.raw_bg.keys():
                XX = self.raw_bg["Oxygen"][ii,i]
                msk = ~XX.mask
                if np.nansum(msk)>5:
                    mld_DO_0 = mixed_layer_depth(z0[msk], -np.array([XX[msk]]).T, Dd = DO_ml_crit)[0]
                    mld_DO_1 = mixed_layer_depth(z0[msk], np.array([XX[msk]]).T, Dd = DO_ml_crit)[0]
                    self.gr["mld_DO"][i] = np.nanmin(np.array([mld_DO_0,mld_DO_1]))
                    #self.gr["mld_DO"][i] = mixed_layer_depth(z0[msk], -np.array([XX[msk]]).T, Dd = DO_ml_crit, crit = "DO")[0]

        self.gr["gpa"] = gsw.geo_strf_dyn_height(self.gr["SA"], self.gr["CT"], self.gr["Pressure"], interp_method = "linear", p_ref = 500.)
        self.gr["gpa_500_1500"] = np.full(nt, np.nan)
        for i in range(nt):
            try:
                j =  np.nanargmin(np.abs(self.gr["Pressure"][:,i]-1500. ))
            except:
                j = np.nan
            if np.isnan(j) or np.abs(self.gr["Pressure"][j,i]-1500)>100:
                continue
            self.gr["gpa_500_1500"][i] = -self.gr["gpa"][j,i]
        
        #other derived variables
        self.gr["AOU"] = 100*self.gr["Oxygen"]/self.gr["OxygenSat"]-self.gr["Oxygen"]
        ##calculates PT and SP 
        #self.gr["SP"] = gsw.SP_from_SA( self.gr["SA"], self.gr["Pressure"], self.gr["Lon"], self.gr["Lat"] )
        #self.gr["PT"] = gsw.pt_from_CT( self.gr["SA"], self.gr["CT"] )

    def calculates_carbon_framework(self,**kargs):
        #kargs: CO2file (file for xCO2 data), sp (surface pressure in Pa), timemet (meteo time for surface pressure)
        print("Carbon framework")
        if "CO2file" in kargs:
            CO2args = {"textfile": kargs["CO2file"]}
        else:
            CO2args = {}
        if "ML_zero" in kargs:
            ML_zero = kargs["ML_zero"]
        else:
            ML_zero = True
        intCO2 = reads_CO2_file_cape_grim(interpolation = "linear",plots = False, **CO2args)
        xCO2 = intCO2(self.gr["date"])
        if "sp" in kargs:
            if type(kargs["timemet"])==np.datetime64:
                kargs["timemet"] = convert_datetime64_to_time(kargs["timemet"])
            sp = np.full( self.gr["date"].size, np.nan )
            for i in range(self.gr["date"].size):
                if i == 0:
                    time0 = self.gr["date"][0]-5.
                    if self.gr["date"].size>1:
                        time1 = 0.5*(self.gr["date"][0]+self.gr["date"][1])
                    else:
                        time1 = self.gr["date"][0]+5.
                if i==self.gr["date"].size-1:
                    time0 = 0.5*(self.gr["date"][i-1]+self.gr["date"][i])
                    time1 = self.gr["date"][i]+5.
                else:
                    time0 = 0.5*(self.gr["date"][i-1]+self.gr["date"][i])
                    time1 = 0.5*(self.gr["date"][i]+self.gr["date"][i+1])
                ij = np.where( (kargs["timemet"]>=time0) & (kargs["timemet"]<=time1) )[0]
                if ij.size == 0:
                    continue
                sp[i] = np.nanmean(kargs["sp"]/101325.)
            nt = self.gr["date"].size
            nz = self.gr["depth"].size
            zM = np.tile(self.gr["depth"],(nt,1)).T
            mldM = np.tile(self.gr["mld"],(nz,1))
            ismld = zM<mldM
            Tml = np.copy(self.gr["CT"])
            Tml[~ismld] = np.nan
            Tml = np.nanmean(Tml, axis = 0)
            Sml = np.copy(self.gr["SA"])
            Sml[~ismld] = np.nan
            Sml = np.nanmean(Sml, axis = 0)
            pH2O = partial_pressure_water_vapour( Sml, Tml )
            pCO2atm = xCO2*(sp - pH2O)
        else:
            pCO2atm = np.copy(xCO2)
        self.gr["CF"] = carbon_framework(self.gr["DIC_LIAR"], self.gr["TALK_LIAR"], self.gr["SA"],\
                                         self.gr["CT"], self.gr["Pressure"], self.gr["Lon"], self.gr["Lat"], \
                                         self.gr["AOU"], pCO2atm,self.gr["depth"], mld = self.gr["mld"], ML_zero = ML_zero)
        self.gr["CF"]["pCO2atm"] = np.copy(pCO2atm)
    
    
    def calculates_CO2_O2_flux(self, met,**kargs):
        if type(met["time"][0]) == np.datetime64:
            met["time"] = convert_datetime64_to_time(met["time"])
        met["Wsp"],met["wind_dir"] = uv_to_wdir( met["u10"], met["v10"] )
        nt = self.gr["date"].size
        nz = self.gr["depth"].size
        zM = np.tile(self.gr["depth"],(nt,1)).T
        mldM = np.tile(self.gr["mld"],(nz,1))
        ismld = zM<mldM
        
        Tml = np.copy(self.gr["CT"])
        Tml[~ismld] = np.nan
        Tml = np.nanmean(Tml, axis = 0)
        iif = np.isfinite(Tml)
        if np.sum(iif)>2:
            intTml = intrp.interp1d( self.gr["date"][iif], Tml[iif], bounds_error = False )
            Tml_met = intTml( met["time"])
            iif = np.where(np.isfinite(Tml_met))[0]
            Tml_met[0:iif[0]] = Tml_met[iif[0]]
            Tml_met[iif[-1]+1:] = Tml_met[iif[-1]]
        else:
            Tml_met = np.nanmean(Tml[iif])*np.ones(met["time"].size)

        Sml = np.copy(self.gr["SA"])
        Sml[~ismld] = np.nan
        Sml = np.nanmean(Sml, axis = 0)
        iif = np.isfinite(Sml)
        if np.sum(iif)>2:
            intSml = intrp.interp1d( self.gr["date"][iif], Sml[iif], bounds_error = False )
            Sml_met = intSml( met["time"])
            iif = np.where(np.isfinite(Sml_met))[0]
            Sml_met[0:iif[0]] = Sml_met[iif[0]]
            Sml_met[iif[-1]+1:] = Sml_met[iif[-1]]
        else:
            Sml_met = np.nanmean(Sml[iif])*np.ones(met["time"].size)

        denml = np.copy(self.gr["Sigma_theta"])
        denml[~ismld] = np.nan
        denml = np.nanmean(denml, axis = 0)
        iif = np.isfinite(denml)
        if np.sum(iif)>2:
            intdenml = intrp.interp1d( self.gr["date"][iif], denml[iif], bounds_error = False )
            denml_met = intdenml( met["time"])
            iif = np.where(np.isfinite(denml_met))[0]
            denml_met[0:iif[0]] = denml_met[iif[0]]
            denml_met[iif[-1]+1:] = denml_met[iif[-1]]
        else:
            denml_met = np.nanmean(denml[iif])*np.ones(met["time"].size)
                

        AOUml = np.copy(self.gr["AOU"])
        AOUml[~ismld] = np.nan
        AOUml = np.nanmean(AOUml, axis = 0)
        iif = np.isfinite(AOUml)
        if np.sum(iif)>10:
            intAOUml = intrp.interp1d( self.gr["date"][iif], AOUml[iif], bounds_error = False )
            AOUml_met = intAOUml( met["time"])
            iif = np.where(np.isfinite(AOUml_met))[0]
            AOUml_met[0:iif[0]] = AOUml_met[iif[0]]
            if iif[-1]>= AOUml_met.size*3./4.:
                AOUml_met[iif[-1]+1:] = AOUml_met[iif[-1]]
        else:
            AOUml_met = np.full(met["time"].size, np.nan)
    
        pCO2ml = np.copy(self.gr["pCO2_LIAR"])
        pCO2ml[~ismld] = np.nan
        pCO2ml = np.nanmean(pCO2ml, axis = 0)
        iif = np.isfinite(pCO2ml)
        if np.sum(iif) > 10:
            intpCO2ml = intrp.interp1d( self.gr["date"][iif], pCO2ml[iif], bounds_error = False )
            pCO2ml_met = intpCO2ml( met["time"])
            iif = np.where(np.isfinite(pCO2ml_met))[0]
            pCO2ml_met[0:iif[0]] = pCO2ml_met[iif[0]]
            if iif[-1]>= pCO2ml_met.size*3./4.:
                pCO2ml_met[iif[-1]+1:] = pCO2ml_met[iif[-1]]
        else:
            pCO2ml_met = np.full(met["time"].size, np.nan)
    
        
        if "CO2file" in kargs:
            CO2args = {"textfile": kargs["CO2file"]}
        else:
            CO2args = {}
        intCO2 = reads_CO2_file_cape_grim(interpolation = "linear",plots = False, **CO2args)
        
        #interpolates CO2
        xCO2met = intCO2(met["time"])
        pH2Oatm = partial_pressure_water_vapour( Sml_met, Tml_met )
        pCO2atm = xCO2met*(met["sp"]/101325. - pH2Oatm)

        K0 = CO2_solubility(Sml_met, Tml_met)
        #gets the CO2 flux
        kwCO2 = kw_wanninkhof(met["Wsp"],Tml_met, gas = "CO2")/100*24. #m/d
        FCO2 = kwCO2*K0*(pCO2ml_met - pCO2atm )*365/1000.*(1000+denml_met)/1000 #umol/kg *m/d *365/1000 ~ mol m-2 y-1 
        #gets the oxygen flux
        kwO2 = kw_wanninkhof(met["Wsp"],Tml_met, gas = "O2")/100*24. #m/d
        FO2 = -kwO2*(AOUml_met)*365/1000.*(1000+denml_met)/1000 #umol/kg *m/d *365/1000~ mmol m-2 d-1  ~ mol m-2 y-1

        self.gr["FCO2"] = np.full(nt, np.nan)
        self.gr["FO2"] = np.full(nt, np.nan)
        for i in range(nt):
            ij = np.where( (np.abs( self.gr["date"][i] - met["time"] )<5.) )[0]
            if ij.size == 0:
                continue
            if np.isnan(pCO2ml[i]) or np.isnan(Tml[i]):
                continue
            #removes data with ice
            if Tml[i]<-1:
                if np.sum( np.isfinite(self.gr["CT"][0:2,i]) ) == 0:
                    continue
            
            self.gr["FCO2"][i] = np.nanmean(FCO2[ij])
            self.gr["FO2"][i] = np.nanmean(FO2[ij])
            
    
    def plots_all_mixing_profiles(self, save = True, show = False):
        nprf = self.raw["date"].size
        for i in range(nprf):
            print("Plot profile %d of %d"%(i+1, nprf))
            self.plots_mixing_layer_profile(i, save = save, show = show)
             
    def plots_mixing_layer_profile(self,pn, save = True, show = False):
        if save:
            if not os.path.exists('prof_ml'):
                os.makedirs('prof_ml')
        date0 = datetime.datetime.fromordinal(int(self.raw["date"][pn]))
        date_str = date0.strftime("%Y %b %d")
        if "Oxygen" in self.raw_bg.keys():
            nsbp = 4
        else:
            nsbp = 3
        xsize = int(np.round(nsbp*2.5))
        fig, ax = plt.subplots(1,nsbp, sharey = True, figsize = (xsize,4))
        ax[0].plot(self.gr["CT"][:,pn],self.gr["depth"],"k-", ms = 2)
        ax[0].plot(self.raw["CT"][:,pn],self.raw["depth"][:,pn],"ko", ms = 2, mfc = "w")
        ax[0].set_ylim(ax[0].get_ylim()[::-1])
        ax[0].set_xlabel("$\\Theta$ [$^{\\mathrm{o}}$C]")
        ax[0].set_ylabel("Depth [m]")
        ax0 = ax[0].twiny()
        ax0.plot(self.gr["SA"][:,pn],self.gr["depth"],"-", color = "gray")
        ax0.plot(self.raw["SA"][:,pn],self.raw["depth"][:,pn],"o", ms = 2, mfc = "w", mec = "gray")
        ax0.set_xlabel("$S_A$", color = "gray")
        
        ax[1].plot(self.gr["Sigma_theta"][:,pn],self.gr["depth"],"k-", ms = 2)
        ax[1].plot( self.raw["Sigma_theta"][:,pn], self.raw["depth"][:,pn],"ko", ms = 2, mfc = "w")
        ax[1].set_xlabel("$\\sigma_{\\theta}$ [kg m$^{-3}$]")

        ax[2].plot(self.raw["size_ov"][:,pn], self.raw["depth"][:,pn], color = "gray", lw = 1)
        ax[2].plot(self.raw["LT_ov"][:,pn], self.raw["depth"][:,pn], color = "k")
        
        ax[2].set_xlabel("$L_T$ (black), $l_{ov}$ (gray)")
        
        if "Oxygen" in self.raw_bg:
            ax[3].plot(self.gr["Oxygen"][:,pn],self.gr["depth"],"k-", ms = 2)
            ax[3].plot( self.raw_bg["Oxygen"][:,pn], self.raw["depth"][:,pn],"ko", ms = 2, mfc = "w")
            ax[3].set_xlabel("DO [$\\mu$mol kg$^{-1}$]")
            ax3 = ax[3].twiny()
            ax3.plot(self.gr["OxygenSat"][:,pn],self.gr["depth"],"-", ms = 2, color = "gray")
            ax3.plot( self.raw_bg["OxygenSat"][:,pn], self.raw["depth"][:,pn],"o", ms = 2, mfc = "w", mec = "gray")
            ax3.set_xlabel("% DO$_{sat}$", color = "gray")
        for ax0 in ax:
            l0 = ax0.axhline(self.gr["mld"][pn], color = cm.tab10(0))
            l1 = ax0.axhline(self.gr["mld_HT"][pn], color = cm.tab10(2))
            l2 = ax0.axhline(self.gr["mld_DO"][pn], color = cm.tab10(3))
            l3 = ax0.axhline(self.gr["h_no_ov"][pn], color = cm.tab10(4))
            l4 = ax0.axhline(self.gr["h_largest_ov"][pn], color = cm.tab10(5))
            l = (l0,l1,l2, l3,l4)
        ax[1].legend(l, ["mld$_{\\sigma_{\\theta}}$","mld$_{\\mathrm{HT}}$","mld$_{\\mathrm{DO}}$","$l_{ov}=0$ m","larg$^{\\mathrm{st}}$. eddy"] )
        fig.suptitle("Float %s, date %s\nLon: %1.2f Lat: %1.2f"%(self.raw["code"], date_str, self.raw["Lon"][pn], self.raw["Lat"][pn]))
        if save:
            date_str0 = date0.strftime("%Y%m%d")
            figname = "prof_ml/%s_%s.png"%(self.raw["code"],date_str0)
            fig.savefig(figname, dpi = 300, bbox_inches = "tight")
        if show:
            plt.show()
        else:
            plt.close(fig)
            
    def plots_map_main_variables(self, saves = True, shows = False,**kargs):
        if not os.path.exists('float_maps'):
            os.makedirs('float_maps')

        if self.raw["Temperature"].shape[1] == 1:
            print("Only one profile")
            return
        #if "weddel_file" in kargs

        
        fig = plt.figure(figsize = (14,8))
        proj = crs.LambertAzimuthalEqualArea(central_latitude=-90.0)
        ax0 = fig.add_axes([0.10,0.67,0.3,0.3], projection = proj)
        ax0.gridlines(draw_labels=False)
        ax0.set_extent([-180, 180, -90, -45], crs.PlateCarree()) # originally -25 as north extent, will shorten to -45
        ax0.stock_img()

        cc = ax0.scatter(self.raw["Lon"], self.raw["Lat"], 20, c = self.raw["date"],transform = crs.PlateCarree(),)#-self.raw["date"][0])
        loc = mdates.AutoDateLocator()
        fig.colorbar(cc, ticks=loc,
                 format=mdates.AutoDateFormatter(loc))
        
        

        """
        ax0 = fig.add_axes([0.10,0.67,0.3,0.3])
        width = 15e6; lon_0 = 0; lat_0 = -90
        m1 = Basemap(width=width,height=width,projection='aeqd',
                     lat_0=lat_0,lon_0=lon_0)
        m1.drawcoastlines()
        m1.fillcontinents()
        m1.drawmapboundary(fill_color='skyblue')
        m1.fillcontinents(color='#cc9966',lake_color='#99ffff')
        m1.drawparallels(np.arange(-80,-20,10),labels=[1,0,0,0])
        m1.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1])
        x,y = m1( self.raw["Lon"], self.raw["Lat"])
        #plt.scatter(x,y,10,T_gr[5,:])
        #plt.plot(x,y,color = "crimson")
        cc = plt.scatter(x,y,20, c = self.raw["date"])#-self.raw["date"][0])
        loc = mdates.AutoDateLocator()
        fig.colorbar(cc, ticks=loc,
                 format=mdates.AutoDateFormatter(loc))
        #cb = fig.colorbar(cc)
        #cb.set_label("Survey day")
        """

        ax1 = fig.add_axes([0.07,0.35,0.47,0.27])
        cfT=ax1.contourf(self.gr["date"], self.gr["depth"], self.gr["CT"],20, cmap = cmocean.cm.thermal)
        #ccT = ax1.contour(self.gr["date"], self.gr["depth"], self.gr["Temperature"],20, colors = "w", linewidths = 1)
        ax1.plot(self.gr["date"], self.gr["mld"], color = "w", lw = 1)
        ax1.plot(self.gr["date"], self.gr["mld_HT"], color = "w", lw = 1, ls = "dotted")
        ax1.plot(self.gr["date"], self.gr["mld_DO"], ls = "--", color = "w", lw = 1)
        ax1.plot(self.gr["date"],1990*np.ones(self.gr["date"].size),marker = "|", color = "k")
        cD = ax1.contour(self.gr["date"], self.gr["depth"], self.gr["gamma_n"],[26.80,27.23,27.50], colors = "skyblue", linewidths = 1)
        plt.clabel(cD, fmt = "%1.2f", fontsize = 6)
        cb = fig.colorbar(cfT)
        ax1.annotate("$\Theta$ [$^{\\mathrm{o}}$C]", xy = (0.02,0.05), xycoords = "axes fraction", fontweight = "bold", color = "k", bbox = dict(facecolor = "w", alpha =0.8))
        if "ylim" in kargs:
            yl = kargs["ylim"]
        else:
            yl = ax1.get_ylim()[::-1]
        ax1.set_ylim(yl)
        ax1.set_ylabel("Depth [m]")
        ax1.set_xticklabels([])
        #ax1.xaxis.set_major_formatter(mdates.AutoDateFormatter(loc))

        ax2 = fig.add_axes([0.07,0.05,0.47,0.27])
        cfT=ax2.contourf(self.gr["date"], self.gr["depth"], self.gr["SA"],20, cmap = cmocean.cm.haline)
        #ccT = ax2.contour(self.gr["date"], self.gr["depth"], self.gr["Salinity"],20, colors = "gray", linewidths = 1)
        ax2.plot(self.gr["date"], self.gr["mld"], color = "w", lw = 1)
        ax2.plot(self.gr["date"], self.gr["mld_DO"], ls = "--",color = "w", lw = 1)
        cb = fig.colorbar(cfT)
        ax2.annotate("$S_A$", xy = (0.02,0.05), xycoords = "axes fraction", fontweight = "bold", color = "k", bbox = dict(facecolor = "w", alpha =0.8) )
        ax2.set_ylim(yl)
        ax2.set_ylabel("Depth [m]")
        ax2.xaxis.set_major_formatter(mdates.AutoDateFormatter(loc))

        """
        ax3 = fig.add_axes([0.54,0.65,0.47,0.27])
        ccT = ax3.pcolor(self.gr["date"], self.gr["depth"], self.gr["LT"], cmap = cm.inferno)
        ax3.plot(self.gr["date"], self.gr["mld"], color = "w", lw = 1)
        ax3.plot(self.gr["date"], self.gr["mld_DO"], ls ="--",color = "w", lw = 1)
        plt.colorbar(ccT, ax = ax3)
        ax3.set_ylim(yl)
        ax3.set_ylabel("Depth [m]")
        ax3.set_xticklabels([])
        ax3.annotate("$L_T$ [m]", xy = (0.02,0.05), xycoords = "axes fraction", fontweight = "bold", color = "k", bbox = dict(facecolor = "w", alpha =0.8))
        ax3.set_title("Float: %s"%(self.raw["code"]))
        """
        
        if "Nitrate" in self.gr.keys():
            ax3 = fig.add_axes([0.54,0.65,0.47,0.27])
            ccT = ax3.contourf(self.gr["date"], self.gr["depth"], self.gr["Nitrate"], 20, cmap = cmocean.cm.matter)
            ax3.plot(self.gr["date"], self.gr["mld"], color = "w", lw = 1)
            ax3.plot(self.gr["date"], self.gr["mld_DO"], ls ="--",color = "w", lw = 1)
            plt.colorbar(ccT, ax = ax3)
            ax3.set_ylim(yl)
            ax3.set_ylabel("Depth [m]")
            ax3.set_xticklabels([])
            ax3.annotate("Nitrate [$\\mu$mol kg$^{-1}$]" , xy = (0.02,0.05), xycoords = "axes fraction", fontweight = "bold", color = "k", bbox = dict(facecolor = "w", alpha =0.8))
            ax3.set_title("Float: %s"%(self.raw["code"]))
            
        
        
        if "Oxygen" in self.gr.keys():
            ax4 = fig.add_axes([0.54,0.35,0.47,0.27])
            cfT=ax4.contourf(self.gr["date"], self.gr["depth"], self.gr["Oxygen"]-100*self.gr["Oxygen"]/self.gr["OxygenSat"],20, cmap = cmocean.cm.oxy)
            #ccT = ax2.contour(self.gr["date"], self.gr["depth"], self.gr["Salinity"],20, colors = "gray", linewidths = 1)
            ccT = ax4.contour(self.gr["date"], self.gr["depth"], self.gr["Oxygen"]-100*self.gr["Oxygen"]/self.gr["OxygenSat"],[0], colors = "blue", linewidths = 1)
            ax4.plot(self.gr["date"], self.gr["mld"], color = "k", lw = 1)
            ax4.plot(self.gr["date"], self.gr["mld_DO"], ls = "--", color = "k", lw = 1)
            cb = fig.colorbar(cfT)
            ax4.annotate("DO-DO$_{\\mathrm{sat}}$ [$\\mu$ mol kg$^{-1}$]", xy = (0.02,0.05), xycoords = "axes fraction", fontweight = "bold", color = "k", bbox = dict(facecolor = "w", alpha =0.8))
            ax4.set_ylim(yl)
            ax4.set_yticklabels([])
            ax4.set_xticklabels([])


        if "DIC_LIAR" in self.gr.keys():
            ax5 = fig.add_axes([0.54,0.05,0.47,0.27])
            cfT=ax5.contourf(self.gr["date"], self.gr["depth"], self.gr["DIC_LIAR"],20, cmap = cmocean.cm.ice_r)
            #ccT = ax2.contour(self.gr["date"], self.gr["depth"], self.gr["Salinity"],20, colors = "gray", linewidths = 1)
            #ccT = ax2.contour(self.gr["date"], self.gr["depth"], self.gr["DIC_LIAR"],[0], colors = "gray", linewidths = 1)
            ax5.plot(self.gr["date"], self.gr["mld"], color = "k", lw = 1)
            ax5.plot(self.gr["date"], self.gr["mld_DO"], ls = "--", color = "k", lw = 1)
            cb = fig.colorbar(cfT)
            ax5.annotate("DIC [$\\mu$ mol kg$^{-1}$]", xy = (0.02,0.05), xycoords = "axes fraction", fontweight = "bold", color = "k", bbox = dict(facecolor = "w", alpha =0.8))
            ax5.set_ylim(yl)
            ax5.set_yticklabels([])
            ax5.xaxis.set_major_formatter(mdates.AutoDateFormatter(loc))
        filename = "float_maps/%s_map.png"%(self.raw["code"])
        if saves:
            fig.savefig(filename)
            plt.close(fig)
        if shows:
            plt.show()

        
def grids_interpolates(x0,y0,x,dx, grid = False):
    y = np.full(x.size,np.nan)
    if grid:
        for i in range(x.size):
            jj = (x0>=x[i]-dx/2.) & (x0<=x[i]+dx/2.)
            if np.nansum(jj)>0:
                y[i] = np.mean(y0[jj])

        igood = np.isfinite(y)
        if np.sum(igood)>5:
            intt = intrp.interp1d( x[igood], y[igood], bounds_error = False)
            y[~igood] = intt(x[~igood])
    elif np.sum(np.isfinite(y0))>5:
        intt = intrp.interp1d( x0, y0, bounds_error = False)
        y = intt(x)
    return y


##############################
######### OTHER FUNCTIONS ####
##############################

def mixed_layer_depth(z0, den0, Dd = 0.03, crit = "diff", z_min = 30., intrp = True):
    #Mixed layer calculation
    if crit != "diff" and crit != "grad" and crit != "DO": # != is callled bang, and it is asking if something is not equal to something else
        crit = "diff"
        print("Incorrect criterion, set to diff")
    c,f = den0.shape
    MLD = np.full(f, np.nan)
    for i in range(f):
        if z0.ndim ==1:
            z = np.copy(z0)
        else:
            z = z0[:,i]
        #den = np.sort(den0[:,i])
        den = den0[:,i]
        iif = np.isfinite(den+z)
        if np.sum(iif)<=1:
            continue
        den = den[iif]
        z = z[iif]
        if np.min(z0)>z_min:
            continue

        if crit == "diff":
            sden = den[0]
            denp = den-sden
            imld = np.where( denp>=Dd )[0]
            if imld.size == 0:
                MLD[i] = np.max(z)
            elif imld[0]>0:
                imld = imld[0]
                z2 = z[imld]
                z1 = z[imld-1]
                denp2 = denp[imld]
                denp1 = denp[imld-1]
                if intrp:
                    MLD[i] = (z2-z1)/(denp2-denp1)*(Dd - denp1) + z1
                else:
                    MLD[i] = (z1+z2)*0.5
            else:
                MLD[i] = np.max(z)
                #MLD[i] = z0[0,i]
                
        elif crit == "grad":
            grden = np.abs(first_centered_differences(z, den))
            imld = np.where(grden>=Dd)[0]
            if imld.size == 0:
                MLD[i] = np.max(z)
            elif imld[0]>0:
                imld = imld[0]
                z2 = z[imld]
                z1 = z[imld-1]
                grd2 = grden[imld]
                grd1 = grden[imld-1]
                if intrp:
                    MLD[i] = (z2-z1)/(grd2-grd1)*(Dd - grd1) + z1
                else:
                    MLD[i] = 0.5*(z1+z2)
            else:
                MLD[i] = z[0]

        if crit == "DO":
            sden = den[0]
            denp = den-sden
            imld = np.where( np.abs(denp)>=Dd )[0]
            if imld.size == 0:
                MLD[i] = np.max(z)
            elif imld[0]>0:
                imld = imld[0]
                z2 = z[imld]
                z1 = z[imld-1]
                MLD[i] = z1
            else:
                MLD[i] = np.max(z)
                #MLD[i] = z0[0,i]
                
    return MLD


def calculates_thorpe_scale(z,dens,PLOT = False):
    #sorts for ascending depth
    ii = np.argsort(z)
    z = z[ii]
    dens = dens[ii]

    #sorts for ascending density
    jj = np.argsort(dens)
    disp = z - z[jj]

    nn = disp.size

    #Looks for individual overturns
    LT = np.zeros(nn)
    ov_size = np.zeros(nn)
    ov_num = np.zeros(nn)
    ovN0 = 1
    i = 0
    while True:
        #plt.plot(dens[i:]-dens[i])
        ii_lighter0 = np.where( (dens[i:]-dens[i])<=0 )[0]        
        if ii_lighter0.size>1:
            ii_lighter = np.arange(i,i+ii_lighter0[-1]+1)
            #print(ii_lighter0)
            dens_ov = dens[ii_lighter]
            z_ov = z[ii_lighter]
            jj = np.argsort(dens_ov)
            disp_ov = z_ov - z_ov[jj]
            #print(disp_ov)
            LT[ii_lighter] = np.nanmean(disp_ov**2)**0.5
            if LT[ii_lighter][0]>0:
                ov_size[ii_lighter] = np.max(z_ov)-np.min(z_ov)
                ov_num[ii_lighter] = ovN0
                ovN0+=1
            i = ii_lighter[-1]+1
        else:
            i+=1
        if i>=nn:
            break
    if PLOT == True:    
        fig, ax = plt.subplots(1,2, sharey = True)
        ax[0].plot(dens, z)
        ax[0].set_ylim(ax[0].get_ylim()[::-1])
        ax[0].set_xlabel("$\\sigma_{\\theta}$ [kg m$^{-3}$]")
        ax[0].set_ylabel("Depth [m]")
        ax[1].plot(np.abs(disp),z, lw = 1, color = "gray")
        ax[1].plot(LT,z, color = "k")
        #ax[1].plot(ov_size,z)
        ax[1].set_xlabel("$L_T$ [m]")
        plt.show()

    return z, LT, ov_size, ov_num


def geopotential_anomaly(CT,SA,p, pref = np.array([500.,1500.])):
    rho = gsw.rho(SA,CT,p)
    rho0 = gsw.rho(35.,0.,p)
    delta = rho**-1 - rho0**-1
    #delta = gsw.specvol_anom_standard(SA,CT,p+10)
    if np.max(p)<np.max(pref):
        return np.nan
    p_i = np.arange(pref[0], pref[1]+1.,1.)
    dp = 1.*1e4 #Pa
    intd = intrp.interp1d( p, delta, bounds_error = False )
    delta_i = intd( p_i )
    gpa = np.sum(dp*delta_i)
    return gpa



def FCD_2d(x, y, axis = 0):
    if x.ndim != 2 or y.ndim !=2:
        sys.exit("Invalid dimensions")
    if axis != 0 and axis != 1:
        sys.exit("Invalid axis")
    if axis == 1:
        x = x.T
        y = y.T
    dy = np.full(y.shape,np.nan)    
    for i in range(x.shape[1]):
        dy[:,i] = first_centered_differences(x[:,i], y[:,i])
        
    if axis == 1:
        dy = dy.T
    return dy

def first_centered_differences(x, y, fill = False):
    if x.size != y.size:
        print("first-centered differences: vectors do not have the same size")
    dy = np.full( x.size, np.nan )
    iif = np.where( (np.isfinite(x)) & (np.isfinite(y))) [0]
    if iif.size < 2:
        return dy
    x0 = x[iif]
    y0 = y[iif]
    dy0 = np.full( x0.size, np.nan )
    #calculates differences
    dy0[0] = (y0[1] - y0[0])/(x0[1]-x0[0])
    dy0[-1] = (y0[-1] - y0[-2])/(x0[-1]-x0[-2])
    dy0[1:-1] = (y0[2:] - y0[0:-2])/(x0[2:]- x0[0:-2])
    
    dy[iif] = dy0

    if fill:
        dy[0:iif[0]] = dy[iif[0]]
        dy[iif[-1]+1:] = dy[iif[-1]]
    return dy


def moving_average(x,n, window = "flat"):
    if n%2 == 0:
        n+=1
    N = x.size
    cx = np.full(x.size, np.nan)
    for i in range(N):
        ii = np.arange(i-n//2, i+n//2+1,1)
        if window == "flat":
            ww = np.ones(ii.size)
        elif window == "gauss":
            xx = ii - i
            
            ww = np.exp(- xx**2/(float(n)/4)**2 )
        elif window == "hanning":
            ww = np.hanning(ii.size)
        ww = ww[ (ii>=0) & (ii<N)]
        ii = ii[ (ii>=0) & (ii<N)]
        
        kk = np.isfinite(x[ii])
        if np.sum(kk)<0.25*ii.size:
            continue
        cx[i] = np.sum(x[ii[kk]]*ww[kk])/np.sum(ww[kk])
    return cx

#time conversion
def convert_time_to_date(time):
     date = [datetime.datetime.fromordinal(int(time0)) + datetime.timedelta(time0%1) for time0 in  time] 
     return date
def convert_date_to_time(date):
     N = len(date)
     time = np.full(N, np.nan)
     for i in range(N):
          time[i]=date[i].toordinal() + date[i].hour/24. + date[i].minute/24./60. + date[i].second/24./60./60. + date[i].microsecond/24./60./60./1e6
     return time
def convert_datetime64_to_date(date64):
    ts = (date64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    date = [datetime.datetime.utcfromtimestamp(ts0) for ts0 in ts]
    return date
def convert_datetime64_to_time(date64):
    ts = (date64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    date = [datetime.datetime.utcfromtimestamp(ts0) for ts0 in ts]
    time = convert_date_to_time(date)
    return time

####
### Meteo functions
###
#wind transformations
def wdir_to_uv(w,alpha):
    alpha = 270.-alpha
    alpha *=np.pi/180
    
    u = w*np.cos(alpha)
    v = w*np.sin(alpha)
    return u,v

def uv_to_wdir(u,v):
    w = (u**2+v**2)**0.5
    alpha = 180/np.pi*np.arctan2(v,u)
    alpha = 270.-alpha
    alpha[alpha>360]-=360
    #alpha[alpha>180] = 360 - alpha[alpha>180]
    return w, alpha

def cd_large_and_pond( U10 ):
    #drag coefficient from Large and Pond 1981
    CD = np.full(U10.size, np.nan)
    CD[U10<11.] = 1.2
    CD[U10>=11.] = 0.49 + 0.065*U10[U10>=11.]
    CD *=1e-3
    return CD

class ERAmeteo():
    def __init__(self, folder, **kargs):
        if "t_chunks" in kargs:
            t_chunks = kargs["t_chunks"]
        else:
            t_chunks = 24
        filelist = sorted(glob.glob(folder + "/*.nc"))
        self.DATASET = xr.open_mfdataset(filelist, parallel = True, chunks = {"time":t_chunks})#, "latitude": 28, "longitude": 144})#
        #display(self.DATASET)
        #self.time = self.DATASET.time.data
        
    
    def get_data(self, date_fl, lon_fl, lat_fl, VARS =  ['u10', 'v10', 't2m', 'mslhf', 'msnlwrf', 'msnswrf', 'msshf', 'sst', 'sp']):
        #transforms time coordinates
        dt_fl = mdates.num2date(date_fl)
        dt64_fl = np.array([np.datetime64(dt0) for dt0 in dt_fl])
        DATASET = self.DATASET.sel( time = slice(dt64_fl[0]-np.timedelta64(10,"D"), dt64_fl[-1]+np.timedelta64(1,"D")))
        DATASET = DATASET.sel(longitude = slice(np.nanmin(lon_fl)-0.5, np.nanmax(lon_fl)+0.5))
        DATASET = DATASET.sel(latitude = slice(np.nanmax(lat_fl)+0.5, np.nanmin(lat_fl)-0.5)) 
        display(DATASET)
        timeERA = DATASET.time.data
        ntE = timeERA.size
        ntF = date_fl.size

        self.ERAhr = dict()    
        for vv in VARS:
            self.ERAhr[vv] = np.full(ntE, np.nan)
        self.ERAhr["time"] = timeERA 
        
        self.ERAlr = dict()    
        for vv in VARS:
            self.ERAlr[vv] = np.full(ntF, np.nan)
        self.ERAlr["time"] = timeERA 
        #np.datetime64(datetime.utcnow()).astype(datetime)
        
        #interpolated coordinates
        for i in range(ntF):
            if i == 0:
                lon = lon_fl[i]
                lat = lat_fl[i]
                time1 = dt64_fl[i]
                time0 = dt64_fl[i] -np.timedelta64(10,"D")
            else:
                lon = 0.5*(lon_fl[i]+lon_fl[i-1])
                lat = 0.5*(lat_fl[i]+lat_fl[i-1])
                time1 = dt64_fl[i]
                time0 = dt64_fl[i-1]
                if time1-time0>np.timedelta64(15,"D"):
                    time0 = time1 - np.timedelta64(10,"D")
            time_count00 = time.time()
            print("\nREADING METEO FLOAT %d of %d (%1.2f %%)"%(i+1, ntF, (i)/float(ntF)*100))
            print("Float time: %s, Long: %1.2f, Lat: %1.2f"%( dt64_fl[i].astype(datetime.datetime).strftime("%Y/%m/%d %H:%M"), lon_fl[i], lat_fl[i] ))
            ii = np.where( (timeERA>time0) & (timeERA<=time1))[0]    
            DT = DATASET.sel(time = slice(time0,time1),expver = 1)
            print("Time for search: %s to %s"%( DT.time.data[0],DT.time.data[-1]))
            DT = DT.sel( longitude = lon, latitude = lat,  method = "nearest" )
            #DT = DT.compute()
            jj = np.where( (DT.time.data>time0) & (DT.time.data<=time1))[0]  
            for vv in VARS:
                #print(vv)
                #display(DT[vv])
                self.ERAhr[vv][ii] = DT[vv].compute().data[jj]
                self.ERAlr[vv][i] = np.nanmean(self.ERAhr[vv][ii])
                #print(self.ERAlr[vv][i])
            print("Elapsed time %1.1f s"%( time.time()-time_count00 ))
        print("READING ENDED")
        
##
## GAS FLUXES
##
def CO2_solubility(S,T):
    #CO2 [umol/kg/atm] solubility in seawater according to Weiss 1974
    #See McGillis and Wanninkhof (2006). Marine Chemistry 98:100-108
    Tk=T+273.15 # in Kelvin
    lnK0 = -60.2409+93.4517*(100/Tk)+23.3585*np.log(Tk/100)+ S*(0.023517-0.023656*(Tk/100)+0.0047036*(Tk/100)**2)
    K0 = np.exp(lnK0)
    return K0
 
def partial_pressure_water_vapour(S,T):
    #Partial pressure of water vapour [atm]
    #See McGillis and Wanninkhof (2006). Marine Chemistry 98:100-108
    #it is used to calculate pCO2 from dry air molecular fraction (X) as:
    # pCO2 = (P - pH2O)*X
    # see Woolf et al. (2016) J. Geophys. Res: Oceans, 121 (2) : 1229-1248
    Tk=T+273.15
    pH2O =  np.exp( 24.4543 - 67.4509*(100/Tk) - 4.8489*np.log(Tk/100) - 0.000544*S )
    return pH2O

def reads_CO2_file_cape_grim(textfile = "atm_CO2/CapeGrim_CO2.csv", interpolation = "linear", plots = False):
    ff = open(textfile)
    time = []
    XCO2 = []
    for line in ff.readlines():
        lineS = line.split(",")
        if "Y" not in lineS[0]:
            date0 = datetime.datetime(int(lineS[0]),int(lineS[1]),int(lineS[2]))
            time0 =date0.toordinal()
            XCO20 = float(lineS[4])
            time.append(time0)
            XCO2.append(XCO20)
    time = np.array(time)
    XCO2 = np.array(XCO2)
    if interpolation == "linear":
        intCO2 = intrp.interp1d( time, XCO2, bounds_error = False )
    elif interpolation == "spline":
        intCO2 = intrp.UnivariateSpline( time, XCO2 )
    if plots:    
        xtime = np.arange( np.nanmin(time) , np.nanmax(time) ,1. )
        fig, ax = plt.subplots()
        ax.plot(time, XCO2,".")
        ax.plot(xtime, intCO2(xtime.astype(float)))
        plt.show()
    return intCO2 
    
def kw_wanninkhof(U, T, gas = "CO2"):
    #wanninkhof 2014 piston velocity
    kw =  0.251*U**2
    Sc = Schmidt_number(T, gas)
    kw = kw * (Sc/660)**-0.5
    return kw
    
def Schmidt_number(T, gas = "CO2"):
    if gas == "CO2":
        Scp = np.poly1d( [0.0007555,-0.0923207, 4.7353, -136.25, 2116.8] )
    elif gas == "O2":
        Scp = np.poly1d( [0.00093777,-0.10939, 5.2122, -135.6, 1920.4] )
    Sc = Scp(T)
    return Sc
        
    
##
## Carbon framework
##
def carbon_framework(DIC, TALK, SA, CT, pres, lon, lat, AOU, pCO2, depth, **kargs):
    import PyCO2SYS as pyco2
    if "ML_zero" in kargs:
        ML_zero = kargs["ML_zero"]
    else:
        ML_zero = True
    if "prealk_eqn" in kargs:
        prealk_eqn = kargs["prealk_eqn"]
    else:
        prealk_eqn = "GLODAP"
    RNO = - 16./170
    RCO = - 106./170.

    CF = dict()
        
    #calculates PT and SP for pyco2
    SP = gsw.SP_from_SA( SA, pres, lon, lat )
    PT = gsw.pt_from_CT( SA, CT )
    
    #function for surface alkalinity
    if prealk_eqn == "GLODAP":
        ALKpre = 42.5036*SP +825.1583
    elif prealk_eqn == "SOCCOM":
        ALKpre = 2818.56 - 80.81*SA - 4.74 * CT  + 1.922 * SA**2 + 0.117 * CT**2
    ##"2818.56 - 80.81*SA - 4.74 * CT  + 1.922 * SA**2 + 0.117 * CT**2"
    #ALKpre = eval("%s"%(prealk_eqn))

    
    #preindustrial saturation
    results = pyco2.sys(ALKpre, 278., 1,4, salinity = SP, temperature = PT)
    CF["DICsat_prein"] = results["dic"]

    #results = pyco2.sys(ALKpre, pCO2, 1,4, salinity = SP, temperature = PT)
    #CF["DICsat_prealk"] = results["dic"]
    
    #present day saturation
    #results = pyco2.sys(ALKpre, pCO2, 1,4, salinity = SP, temperature = PT)
    results = pyco2.sys(ALKpre, pCO2, 1,4, salinity = SP, temperature = PT)
    CF["DICsat"] = results["dic"]
    
    #with local alkalinity
    results = pyco2.sys(TALK, pCO2, 1,4, salinity = SP, temperature = PT)
    CF["DICsat_talk"] = results["dic"]
    
    
    
    #soft tissue
    CF["DICsoft"] = - RCO*AOU
    
    #carbonate
    CF["DICcarb"] = 0.5*( TALK - ALKpre - RNO * AOU )

    if ML_zero and "mld" in kargs:
        mld = kargs["mld"]
        nt = mld.size
        nz = depth.size
        #gets indices for mixed layer
        zM = np.tile(depth,(nt,1)).T
        mldM = np.tile(mld,(nz,1))
        ismld = zM<mldM
        CF["DICsoft"][ismld] = 0.
        CF["DICcarb"][ismld] = 0.
    
    #DeltaC referenced to pre-industrial levels
    CF["DICdelta_prein"] = DIC - CF["DICsat_prein"] - CF["DICsoft"] - CF["DICcarb"]
    #DeltaC referenced to present day
    CF["DICdelta"] = DIC - CF["DICsat"] - CF["DICsoft"] - CF["DICcarb"]
    #Disequilibrium C preindustrial
    CF["DICdis_prein"] = DIC - CF["DICsat_prein"]
    #Disequilibrium C present day
    CF["DICdis"] = DIC - CF["DICsat"]
    #disequilibrium with local talk
    CF["DICdis_talk"] = DIC - CF["DICsat_talk"]
    
    CF["DIC"] = np.copy(DIC)
    CF["ALKpre"] = np.copy(ALKpre)
    
    return CF


###
### Net ecosystem production
def NEP_calculation(date, z, Lon, Lat, Nitrate, POC, SA, **kargs):
    ##FUNCTION TO CALCULATE NEP from nitrate depletion / POC accumulation
    if "PLOT" in kargs:
        PLOT = kargs["PLOT"]
    else:
        PLOT = False
        
    #first I convert the numerical date to a datetime format so I can get the month and year vectors
    RCN = 106/16. # Redfield ratio

    nt = date.size
    dateDT = convert_time_to_date( date )
    year = np.full( nt, np.nan )
    month = np.full(nt, np.nan)
    for i in range(nt):
        year[i] = dateDT[i].year
        month[i] = dateDT[i].month

    #integration depth
    if "mld" in kargs:
        H = min([np.nanmax(kargs["mld"]),500]) # calculates the maximum ML
        #print("Integration depth: %1.0f m"%(H))
    elif "H" in kargs:
        H = kargs["H"]
    else:
        H = 200.

    
    jh = np.where( z>= H)[0][0] # gets the depth index for the maxmum mixed layer

    #depth integrated nitrate
    dint_Nitrate = np.nanmean(Nitrate[:jh,:], axis = 0)*H*(1027/1e6)
    dint_POC = np.nanmean(POC[:jh,:], axis = 0)*H/1000.
    mSA = np.nanmean( SA[z>500,:], axis = 0 )
    #by multiplying by density ~1027 and dividing by 1e6 I get units mol m-2

    #for each year calculates the maximum and minimum
    Uyear = np.unique(year)
    nyr = Uyear.size
    date_nit_sum = np.full(nyr, np.nan)
    date_nit_win = np.full(nyr, np.nan)
    nit_win = np.full(nyr, np.nan)
    nit_sum = np.full(nyr, np.nan)
    nit_win_month_avg = np.full(nyr, np.nan)
    nit_sum_month_avg = np.full(nyr, np.nan)

    POC_win = np.full(nyr, np.nan)
    POC_sum = np.full(nyr, np.nan)
    POC_win_month_avg = np.full(nyr, np.nan)
    POC_sum_month_avg = np.full(nyr, np.nan)

    SA_win = np.full(nyr, np.nan)
    SA_sum = np.full(nyr, np.nan)
    Lat_win = np.full(nyr, np.nan)
    Lat_sum = np.full(nyr, np.nan)
    Lon_win = np.full(nyr, np.nan)
    Lon_sum = np.full(nyr, np.nan)
    flag_nit_NEP = np.full(nyr, False)
    for i, yr in enumerate(Uyear):
        #start_summer = datetime.datetime(int(yr),12,1,0,0).toordinal()
        #end_summer = datetime.datetime(int(yr)+1,4,1,0,0).toordinal()
        start_summer = datetime.datetime(int(yr)+1,1,1,0,0).toordinal()
        end_summer = datetime.datetime(int(yr)+1,4,1,0,0).toordinal()
        it_summer = np.where( (date>= start_summer) & (date<= end_summer) )[0]
        if it_summer.size > 0:
            if np.sum(np.isfinite(dint_Nitrate[it_summer]))>0:
                imin_nit = it_summer[ np.nanargmin( dint_Nitrate[it_summer] ) ]
                date_nit_sum[i] = date[imin_nit]  
                nit_sum[i] =np.nanmin( dint_Nitrate[it_summer])
                POC_sum[i] = dint_POC[imin_nit]
                #ii_sum_month = np.where( np.abs(date - date[imin_nit]  )<15 )[0]
                ii_sum_month = np.where( (month == month[imin_nit]) & (year == year[imin_nit])   )[0]
                nit_sum_month_avg[i] =np.nanmean( dint_Nitrate[ii_sum_month])
                POC_sum_month_avg[i] =np.nanmean( dint_POC[ii_sum_month])
                SA_sum[i] = mSA[imin_nit]
                Lat_sum[i] = Lat[imin_nit]
                Lon_sum[i] = Lon[imin_nit]

        #start_winter = datetime.datetime(int(yr),5,1,0,0).toordinal()
        #end_winter = datetime.datetime(int(yr),12,1,0,0).toordinal()
        start_winter = datetime.datetime(int(yr),8,1,0,0).toordinal()
        end_winter = datetime.datetime(int(yr),12,1,0,0).toordinal()
        it_winter = np.where( (date>= start_winter) & (date<= end_winter) )[0]
        if it_winter.size > 0:
            if np.sum(np.isfinite(dint_Nitrate[it_winter]))>0:
                imax_nit = it_winter[ np.nanargmax( dint_Nitrate[it_winter] ) ]
                date_nit_win[i] = date[imax_nit]  
                nit_win[i] = np.nanmax( dint_Nitrate[it_winter])
                POC_win[i] =  dint_POC[imax_nit]
                #ii_win_month = np.where( np.abs(date - date[imax_nit]  )<15 )[0]
                ii_win_month = np.where( (month == month[imax_nit]) &  (year == year[imax_nit])   )[0]
                nit_win_month_avg[i] =np.nanmean( dint_Nitrate[ii_win_month])
                POC_win_month_avg[i] =np.nanmean( dint_POC[ii_win_month])
                SA_win[i] = mSA[imax_nit]
                Lat_win[i] = Lat[imax_nit]
                Lon_win[i] = Lon[imax_nit]

    flag_NEP = (np.abs(date_nit_win-date_nit_sum)<8*30) & (np.abs(SA_win-SA_sum)<0.05) & (np.abs(Lon_win-Lon_sum)<8.) & (np.abs(Lat_win-Lat_sum)<5.)


    #calculates net ecosystem production (molC m-2 yr-1)
    NEP = (nit_win - nit_sum)*RCN
    #from the monthly means
    NEP_avg = (nit_win_month_avg - nit_sum_month_avg)*RCN

    NEP_POC = -(POC_win - POC_sum)
    NEP_POC_avg = -(POC_win_month_avg - POC_sum_month_avg)

    #gets the date around the depletion
    date_NEP = 0.5*(date_nit_sum +date_nit_win )
    Lon_NEP = 0.5*(Lon_win+Lon_sum)
    Lat_NEP = 0.5*(Lat_win+Lat_sum)
    

    if PLOT:        
        print( "\n-------------------------------------------------------------------------")
        print("YEAR\t    NEP Nit\t     <NEP Nit>\t      NEP POC\t      <NEP POC>" )
        print("\t\t\t\t [mol/m2/yr]")
        print( "-------------------------------------------------------------------------")
        for i in range(nyr):
            print("%d-%d\t %1.2f\t\t%1.2f\t\t%1.2f\t\t%1.2f"%(Uyear[i],Uyear[i]+1, NEP[i], NEP_avg[i], NEP_POC[i], NEP_POC_avg[i])  )
        print( "-------------------------------------------------------------------------")
        print("Mean     \t%1.2f\t\t%1.2f\t\t%1.2f\t\t%1.2f"%(np.nanmean(NEP), np.nanmean(NEP_avg),np.nanmean(NEP_POC), np.nanmean(NEP_POC_avg)))
        print( "-------------------------------------------------------------------------")


        #Plots the results 
        fig, ax = plt.subplots(3,1,figsize = (8,6), sharex = True)
        ax[0].plot( date, dint_Nitrate, "k" )
        l1,=ax[0].plot(date_nit_sum, nit_sum,"o", ms = 10, mec = "k", color = "goldenrod")
        l2,=ax[0].plot(date_nit_win, nit_win,"o", ms = 10, mec = "k", color = "green")
        for i in range(nyr):
            ax[0].plot([date_nit_sum[i]-15,date_nit_sum[i]+15], [nit_sum_month_avg[i],nit_sum_month_avg[i]], color = "k", zorder = -1)
            ax[0].plot([date_nit_win[i]-15,date_nit_win[i]+15], [nit_win_month_avg[i],nit_win_month_avg[i]], zorder = -1, color = "k")
        yl = ax[0].get_ylim()
        for i in range(nyr):
            ax[0].fill_between( [date_nit_sum[i]-15,date_nit_sum[i]+15], y1 = yl[0], y2 = yl[1], color = l1.get_color(), alpha = 0.3 )
            ax[0].fill_between( [date_nit_win[i]-15,date_nit_win[i]+15], y1 = yl[0], y2 = yl[1], color = l2.get_color(), alpha = 0.3 )

        ax[0].set_ylim(yl)
        ax[0].set_ylabel( "$\\int \\mathrm{Nitrate}\, \\rm d z$\n[mol m$^{-2}$]" )
        ax[0].grid(True)


        ax[1].plot( date, dint_POC, "k" )
        l1,=ax[1].plot(date_nit_sum, POC_sum,"o", ms = 10, mec = "k", color = "goldenrod")
        l2,=ax[1].plot(date_nit_win, POC_win,"o", ms = 10, mec = "k", color = "green")
        for i in range(nyr):
            ax[1].plot([date_nit_sum[i]-15,date_nit_sum[i]+15], [POC_sum_month_avg[i],POC_sum_month_avg[i]], color = "k", zorder = -1)
            ax[1].plot([date_nit_win[i]-15,date_nit_win[i]+15], [POC_win_month_avg[i],POC_win_month_avg[i]], zorder = -1, color = "k")
        yl = ax[1].get_ylim()
        for i in range(nyr):
            ax[1].fill_between( [date_nit_sum[i]-15,date_nit_sum[i]+15], y1 = yl[0], y2 = yl[1], color = l1.get_color(), alpha = 0.3 )
            ax[1].fill_between( [date_nit_win[i]-15,date_nit_win[i]+15], y1 = yl[0], y2 = yl[1], color = l2.get_color(), alpha = 0.3 )

        ax[1].set_ylim(yl)
        ax[1].set_ylabel( "$\\int \\mathrm{POC}\, \\rm d z$\n[mol m$^{-2}$]" )
        ax[1].grid(True)

        ax[2].bar( date_NEP[flag_NEP]-50, NEP[flag_NEP], width = 50, ec = "k", label = "Nit 1-prof"  )
        ax[2].bar( date_NEP[flag_NEP]-30, NEP_avg[flag_NEP], width = 50, ec = "k", label = "Nit month" )
        ax[2].bar( date_NEP[flag_NEP]+30, NEP_POC[flag_NEP], width = 50, ec = "k", label = "POC 1-prof"  )
        ax[2].bar( date_NEP[flag_NEP]+50, NEP_POC_avg[flag_NEP], width = 50, ec = "k", label = "POC month" )
        ax[2].set_ylabel("NEP\n[molC m$^{-2}$ y$^{-1}$]")
        ax[2].legend(loc = "center left", bbox_to_anchor = (1.01,0.5))
        formatter = mdates.DateFormatter("%Y") ### formatter of the date
        locator = mdates.YearLocator() ### where to put the labels
        ax[2].xaxis.set_major_locator(locator)
        ax[2].xaxis.set_major_formatter(formatter)
        ax[2].grid(True)

    return date_NEP, Lon_NEP, Lat_NEP, NEP_avg, NEP_POC_avg, flag_NEP


def oxygen_consumption_rate(date, z, Lon, Lat, Oxygen, SA, **kargs):
    if "PLOT" in kargs:
        PLOT = kargs["PLOT"]
    else:
        PLOT = False
    if "zmax" in kargs:
        zmax = kargs["zmax"]
    else:
        zmax = 500.
    if "zmin" in kargs:
        zmin = kargs["zmin"]
    else:
        zmin = 100.
    RCO = - 106./170.
    
    nt = date.size
    dateDT = convert_time_to_date( date )
    year = np.full( nt, np.nan )
    month = np.full(nt, np.nan)
    for i in range(nt):
        year[i] = dateDT[i].year
        month[i] = dateDT[i].month
        
    dz = z[1]-z[0]
    jh = np.where((z>=zmin) & (z<=zmax))[0]

    #depth integrated O2
    dint_O2 = moving_average( np.nanmean(Oxygen[jh,:], axis = 0),10)
    mSA = np.nanmean( SA[z>500,:], axis = 0 )

    O2 = np.copy(Oxygen)
    nz, nt = O2.shape
    for j in range(nz):
        O2[j,:] = moving_average(O2[j,:],10)

    if "mld" in kargs:
        zM = np.tile(z,(nt,1)).T
        mldM = np.tile(kargs["mld"],(nz,1))
        ismld = zM<mldM
        O2[ismld] = np.nan

    #for each year calculates the maximum and minimum
    Uyear = np.unique(year)
    nyr = Uyear.size
    date_O2_sum = np.full(nyr, np.nan)
    date_O2_win = np.full(nyr, np.nan)
    R_O2 = np.full(nyr, np.nan)

    SA_O2_win = np.full(nyr, np.nan)
    SA_O2_sum = np.full(nyr, np.nan)
    Lat_O2_win = np.full(nyr, np.nan)
    Lat_O2_sum = np.full(nyr, np.nan)
    Lon_O2_win = np.full(nyr, np.nan)
    Lon_O2_sum = np.full(nyr, np.nan)

    for i, yr in enumerate(Uyear):
        start_winter = datetime.datetime(int(yr),8,1,0,0).toordinal()
        end_winter = datetime.datetime(int(yr),12,1,0,0).toordinal()
        it_winter = np.where( (date>= start_winter) & (date<= end_winter) )[0]
        imax_O2 = np.nan
        imin_O2 = np.nan
        if it_winter.size > 0:
            if np.sum(np.isfinite(dint_O2[it_winter]))>0:
                imax_O2 = it_winter[ np.nanargmax( dint_O2[it_winter] ) ]

        start_summer = datetime.datetime(int(yr)+1,1,1,0,0).toordinal()
        end_summer = datetime.datetime(int(yr)+1,4,1,0,0).toordinal()
        it_summer = np.where( (date>= start_summer) & (date<= end_summer) )[0]
        if it_summer.size > 0:
            if np.sum(np.isfinite(dint_O2[it_summer]))>0:
                imin_O2 = it_summer[ np.nanargmin( dint_O2[it_summer] ) ]

        if np.isfinite(imin_O2) and np.isfinite(imax_O2) and (imin_O2>imax_O2):
            iiy = np.arange( imax_O2, imin_O2+1  )
            dO2dt = np.full(nz, 0.)
            for j in jh:
                ox = O2[j,iiy]
                time = date[iiy]
                iif = np.isfinite(ox)
                time = time[iif]
                ox = ox[iif]
                p = np.polyfit(time,ox,1)
                #print(p[0])
                dO2dt[j] = p[0]*(time[-1]-time[0])
            R_O2[i] = np.nansum( 0.5*(dO2dt[1:]+dO2dt[:-1])*1027*1e-6*RCO*(z[1:]-z[:-1]))
            date_O2_win[i] = date[imax_O2]  
            SA_O2_win[i] = mSA[imax_O2]
            Lat_O2_win[i] = Lat[imax_O2]
            Lon_O2_win[i] = Lon[imax_O2]
            date_O2_sum[i] = date[imin_O2]  
            SA_O2_sum[i] = mSA[imin_O2]
            Lat_O2_sum[i] = Lat[imin_O2]
            Lon_O2_sum[i] = Lon[imin_O2]

    flag_O2_R = (np.abs(date_O2_win-date_O2_sum)>3*30) & (np.abs(SA_O2_win-SA_O2_sum)<0.05) & (np.abs(Lon_O2_win-Lon_O2_sum)<8.) & (np.abs(Lat_O2_win-Lat_O2_sum)<5.)        
    date_R = 0.5*(date_O2_sum +date_O2_win )
    Lon_R = 0.5*(Lon_O2_win+Lon_O2_sum)
    Lat_R = 0.5*(Lat_O2_win+Lat_O2_sum)

    
    if PLOT:
        fig, ax = plt.subplots(2,1,figsize = (8,5), sharex = True)  
        ax[0].plot(date, dint_O2,"k")
        yl = ax[0].get_ylim()
        for i in range(date_O2_sum.size):
            ax[0].fill_between( [date_O2_win[i], date_O2_sum[i]], y1 = yl[0], y2 = yl[1],color = "gray" )
        ax[0].set_ylim(yl)
        ax[0].set_ylabel("$\\langle \\mathrm{Oxygen} \\rangle$ [$\\mu$mol kg $^{-1}$]")
        ax[1].bar( date_R[flag_O2_R], R_O2[flag_O2_R], width = 50, ec = "k" )
        ax[1].set_ylabel("R\n[molC m$^{-2}$ y$^{-1}$]")
        #ax[1].legend(loc = "center left", bbox_to_anchor = (1.01,0.5))
        formatter = mdates.DateFormatter("%Y") ### formatter of the date
        locator = mdates.YearLocator() ### where to put the labels
        ax[1].xaxis.set_major_locator(locator)
        ax[1].xaxis.set_major_formatter(formatter)
        ax[1].grid(True)
        
    return date_R, Lon_R, Lat_R, R_O2, flag_O2_R
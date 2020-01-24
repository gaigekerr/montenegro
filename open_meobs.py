#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Open and plot Montenegrin air quality observations.

Notes
-----
    Script opens daily observations of PM10, PM2.5, NO2, and O3 at observing 
    sites in Montenegro (openobs_singleyear) and concatenates into multi-year 
    files (openobs_multiyear). Timeseries of daily, weekly, or monthly averages 
    are plotted.

Author
------
Gaige Hunter Kerr (gaige.kerr@jhu.edu)

Revision History
----------------
    10012019 -- initial version created
    17012019 -- functions to plot daily, monthly time series of species added
    23012019 -- function plottimeseries added as a universal time series 
                plotter instead of having a separate one for daily, montly, or
                weekly timeseries 
"""
path_me = '/Users/ghkerr/montenegro/'
# Run first time 
import sys
if 'mpl' not in sys.modules:
    import matplotlib.font_manager
    # For unicode minus/negative sign implementation
    matplotlib.rcParams['axes.unicode_minus'] = False
    # Change width and thickness of ticks/spines
    matplotlib.rcParams['axes.linewidth'] = 1.0
    matplotlib.rcParams['xtick.major.width'] = 1.0
    matplotlib.rcParams['xtick.minor.width'] = 1.0
    matplotlib.rcParams['ytick.major.width'] = 1.0
    matplotlib.rcParams['ytick.minor.width'] = 1.0

def openobs_singleyear(year, path_me):
    """Read Montenegro MDA8 observations for a given year.
    
    This function is the main workhorse in this script and loads observations
    of NO2, PM10, PM2.5 and O3 for 7 in-situ monitors in Montenegro. No
    averaging is performed. 
    
    
    Parameters
    ----------
    year : int
        Year for which observations are desired
    path_me : str
        Path to directory which has a subdirectory "data" containing 
        observations
        
    Returns
    -------
    xl_ty : pandas.core.frame.DataFrame
        Trace gas/particulate observations; each column corresponds to a 
        different monitor location and species measured. 
    """
    import numpy as np
    import pandas as pd
    # path_me = '/Users/ghkerr/montenegro/'
    file = 'data/2015, 2016, 2017 PM, NO2 i O3.xlsx'    
    # Note that this file contains 'Osmoicasovni maksimumi(maximalne 
    # osmocasovne srednje dnevne vrijednosti).' This translates to 'Eight-hour
    # maximums (maximum eight-hour average daily values)'
    xl = pd.ExcelFile(path_me+file)
    # This is kludgey, and I blame Montenegro for the crazy data formats. 
    # Stations recording concentration are Bar, Niksic, Pljevlja, Nova Varos, 
    # Tivat, Gradina	, and Golubovci. Column names indicate both the station name 
    # AND the species (n.b., the units of ug/m3 in the data files are ommitted
    # from the column names and a dummy/blank column is include for each city) 
    dtype = {'Bar Datum': np.str, 
        'Bar NO2' : np.float64,
        'Bar PM10' : np.float64, 
        'Bar PM2.5' : np.float64,
        'Bar O3' : np.float64,
        'Bar Dummy' : np.float64, 
        'Niksic Datum' : np.str,
        'Niksic NO2' : np.float64,
        'Niksic PM10' : np.float64,
        'Niksic PM2.5' : np.float64,
        'Niksic O3' : np.float64,
        'Niksic Dummy' : np.float64,
        'Pljevlja Datum' : np.str,
        'Pljevlja NO2' : np.float64,
        'Pljevlja PM10' : np.float64, 
        'Pljevlja PM2.5' : np.float64, 
        'Pljevlja Dummy' : np.float64,
        'Nova Varos Datum' : np.str,  
        'Nova Varos NO2' : np.float64, 
        'Nova Varos PM10' : np.float64, 
        'Nova Varos Dummy' : np.float64, 
        'Tivat Datum' : np.str, 
        'Tivat PM2.5' : np.float64,
        'Tivat Dummy' : np.float64,
        'Gradina Datum' : np.str, 
        'Gradina NO2' : np.float64,
        'Gradina O3' : np.float64,
        'Gradina Dummy' : np.float64,
        'Golubovci Datum' : np.str,
        'Golubovci NO2' : np.float64,
        'Golubovci O3' : np.float64}
    # Read specific sheet (corresponding to year) to DataFrame; n.b. 
    # xl.sheet_names will list all available sheets
    xl_ty = xl.parse(str(year), dtype=dtype, index_col=None, header=None,
        sep='\s+', names=list(dtype.keys()), skiprows=4)
    # Delete unwanted columns 
    del xl_ty['Bar Dummy'], xl_ty['Niksic Dummy'], xl_ty['Gradina Dummy'], \
        xl_ty['Pljevlja Dummy'], xl_ty['Nova Varos Dummy'], \
        xl_ty['Tivat Dummy']
    return xl_ty

def openobs_multiyear(syear, eyear, path_me):
    """Open and concatenate yearly Montenegrin observations into a single, 
    multi-year DataFrame. 
    
    Parameters
    ----------
    syear : int
        Starting year of measuring period
    eyear : int
        Ending year of measuring period
    path_me : str
        Path to directory which has a subdirectory "data" containing 
        observations
        
    Returns
    -------
    me_all : pandas.core.frame.DataFrame
        Trace gas/particulate observations for multi-year measuring period; 
        each column corresponds to a different monitor location and species 
        measured.     
    """
    import numpy as np
    import pandas as pd
    years = np.arange(syear,eyear+1, 1)
    me_all = []
    # Loop through years in measuring period 
    for year in years: 
        me_ty = openobs_singleyear(year, path_me)
        me_all.append(me_ty)
    # See pd.concat documentation for more info
    me_all = pd.concat(me_all)
    # Every observing site has its own date column, which *should* all be
    # identical. Add a standardized column and delete individual columns.
    dr = pd.date_range(start='01-01-%s'%syear, end='12-31-%s'%eyear, freq='D')
    me_all['Date Local'] = dr
    me_all.index = me_all['Date Local']    
    del me_all['Bar Datum'], me_all['Niksic Datum'], \
        me_all['Pljevlja Datum'], me_all['Nova Varos Datum'], \
        me_all['Tivat Datum'], me_all['Gradina Datum'], \
        me_all['Golubovci Datum'], me_all['Date Local']
    return me_all 

def plottimeseries(me, species, site, freq, path_me, calcex=False):
    """Function plots time series of site-specific pollutant concentrations for
    the species of interest for the measuring period contained in the 
    Montenegrin observations. If desired, number of exceedances of WHO air 
    quality guidelines is printed to the terminal.    
    
    Parameters
    ----------
    me : pandas.core.frame.DataFrame
        Trace gas/particulate observations for multi-year measuring period; 
        each column corresponds to a different monitor location and species 
        measured.     
    species : str
        The species of interest (options are PM10, PM2.5, O3, or NO2)
    site : str
        The site of interest (options are Bar, Niksic, Pljevlja, Nova Varos, 
        Gradina, Tivat, or Golubovci, but note that observations of some 
        species do not exist for certain sites). 
    freq : str
        The temporal averaging that will be applied to the time series 
        (options are 1D (no averaging performed), 1W, or Y)
    path_me : str
        Path to directory which has a subdirectory "figs" containing 
        observations    
    calcex : bool (optional)
        If True, the number of daily or annual WHO air quality guideline 
        exceedences will be calculated, if applicable
        
    Returns
    -------
    None
    
    References
    ----------
    World Health Organization. Occupational and Environmental Health Team. 
    (‎2006)‎. WHO Air quality guidelines for particulate matter, ozone, nitrogen 
    dioxide and sulfur dioxide: global update 2005: summary of risk 
    assessment. World Health Organization
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    # Register the converters for matatplotlib plotting method (uncessary for 
    # certain versions of Pandas?)
    from pandas.plotting import register_matplotlib_converters
    register_matplotlib_converters()
    # Average to desired frequency (1D = daily, 1W = weekly, M = monthly)
    me_avg = me.resample(freq).mean()
    # Only continue if observations of species exist for site of interest
    if '%s %s'%(site, species) in me_avg.columns:
        me_avg = me_avg['%s %s'%(site, species)]
        # Initialize
        fig = plt.figure(figsize=(8,3))
        ax = plt.subplot2grid((1,1), (0,0), colspan=2)
        ax.plot(me_avg, lw=2, color='k')
        ax.text(0.03, 0.85, '%s'%site, ha='left', transform=ax.transAxes, 
            fontsize=14, zorder=20)
        # Modify date (x axis) ticks for Pandas DataFrame
        # Set xlimits
        ax.set_xlim(me_avg.index[0], me_avg.index[-1])
        # Set monthly locator 
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
        # Set formatter
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
        # Add ylabel
        if species == 'PM10':
            ax.set_ylabel('PM$_{\mathregular{10}}$ [$\mathregular{\mu}$g '+
            'm$^{\mathregular{-3}}$]', fontsize=14)
        elif species == 'PM2.5':
            ax.set_ylabel('PM$_{\mathregular{2.5}}$ [$\mathregular{\mu}$g '+
            'm$^{\mathregular{-3}}$]', fontsize=14)
        elif species == 'NO2':
            ax.set_ylabel('NO$_{\mathregular{2}}$ [$\mathregular{\mu}$g '+
                'm$^{\mathregular{-3}}$]', fontsize=14)    
        elif species == 'O3': 
            ax.set_ylabel('O$_{\mathregular{3}}$ [$\mathregular{\mu}$g '+
                'm$^{\mathregular{-3}}$]', fontsize=14)
        # Save figure with file naming convention 
        # *species*_*site*_*temporalavg*.png. This command assumes that there 
        # is a subdirectory "figs" within the work directory. Note that any 
        # white space or capitalization in the site name is removed along 
        # with periods in the species name (i.e., PM2.5)
        plt.savefig(path_me+'figs/'+'%s_%s_%s.png'
            %(species.replace('.',''), ''.join(site.split()).lower(), freq), 
            dpi=300)
        plt.show()
        # If optional argument calcex = True, the number of times WHO 
        # threshold limits for the species of interest was surpassed for daily 
        # data or annual averages. Note that the species for which we can 
        # calculate the number of exceedances is somewhat limited since we are 
        # working with MDA8 data (e.g., only O3 exceedances are 
        # determined basaed on 8 h, daily maximum values). An in-depth report
        # on WHO pollution targets is listed in the reference above, and a 
        # quick summary can be found at 
        # www.who.int/airpollution/publications/aqg2005/en/
        if calcex==True:
            if species=='O3':
                # 8 h daily maximum guideline level is 100 ug m-3
                exceed = np.where(me_avg > 100)[0].shape[0]
                print('WHO O3 MDA8 guideline of 100 ug m-3 exceeded %d '%exceed+\
                    'times from %s to %s at %s'%(me_avg.index[0], 
                    me_avg.index[-1], site))
    return     

me = openobs_multiyear(2015, 2017, path_me)
species_all = ['PM10', 'PM2.5', 'O3', 'NO2']
sites_all = ['Bar', 'Niksic', 'Pljevlja', 'Nova Varos', 'Gradina', 
              'Golubovci', 'Tivat']
#  Plot daily values 
for species in species_all:
    for site in sites_all:
        plottimeseries(me, species, site, '1D', path_me, calcex=True)
#  Plot weekly values 
for species in species_all:
    for site in sites_all:
        plottimeseries(me, species, site, '1W', path_me, calcex=True)
for species in species_all:
    for site in sites_all:
        plottimeseries(me, species, site, '1M', path_me, calcex=True)

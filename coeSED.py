# filttools.py

import os
import numpy as np

#%matplotlib inline
#%matplotlib notebook
import matplotlib.pyplot as plt
# https://matplotlib.org/tutorials/introductory/customizing.html
#plt.style.use('/Users/dcoe/p/matplotlibrc.txt')
plt.style.use('https://www.stsci.edu/~dcoe/matplotlibrc.txt')
#colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#import matplotlib.ticker as ticker

def getlamname(filt):
    return int(filt[1:4])

def getlam(filt):
    if filt in 'ch1 ch2 ch3 ch4'.split():
        lam = {'ch1':3600, 'ch2':4500, 'ch3':5800, 'ch4':8000}[filt]
    else:
        lam = getlamname(filt)
        if lam < 200:
            lam = 10 * lam
    return lam

def getinstr(filt):
    lam = getlam(filt)
    if lam < 400:
        instr = 'wfc3uvis'
    elif lam < 950:
        instr = 'acs'
    elif lam < 2000:
        instr = 'wfc3ir'
    else:
        instr = 'irac'
    return instr

fullinstrdict = {
    'acs': 'HST_ACS_WFC',
    'ir':  'HST_WFC3_IR',
    'wfc3ir':  'HST_WFC3_IR',
    'wfc3uvis':'HST_WFC3_UVIS',
    'uvis':'HST_WFC3_UVIS',
    'irac':'irac'
    }

def getfullinstr(filt):
    return fullinstrdict[getinstr(filt)]

def getfullfilt(filt):
    instr = getinstr(filt)
    fullinstr = fullinstrdict[instr]
    if filt[:2] != 'ch':
        filt = filt.upper()
    fullfilt = fullinstr + '_' + filt
    return fullfilt

def getfullfiltfile(filt):
    filter = getfullfilt(filt)
    return os.path.join('FILTER', filter+'.res')

def extract_id1(tbl, id, idlabel='id'):
    # Note input id must be same format as tbl[idlabel] (int or string)
    return tbl[tbl[idlabel] == int(id)][0]

def extract_id(cat, id, idlabel='id'):
    # Note input id must be same format as cat[idlabel] (int or string)
    # to create an array with mostly False entries and True for id
    # cat['id'] dtype='int64'
    # cat['id'] == 3: [False, False, True, False...]
    # cat['id'] == '3': False
    #
    # obj:          astropy.table.table.Table
    # obj[0]:       astropy.table.row.Row
    # obj['id']:    astropy.table.column.Column
    # obj[0]['id']: numpy.int64
    duck_duck_goose = cat[idlabel] == int(id)
    if len(duck_duck_goose):
        obj = cat[duck_duck_goose]
    return obj[0]

def extract_filters(cat, flux_suffix='_flux'):
    filters = []
    for label in cat.columns:
        if label.endswith(flux_suffix):
            filt = label[:-len(flux_suffix)]
            filters.append(filt)

    return(filters)

input_flux_units = 'uJy'

# output -> BAGPIPES
# can only have one input: id
# so either copy this into your notebook
# or pass catalog and filters into this module:
# coeSED.cat = cat
# coeSED.filters = filters
def load_phot(id, flux_suffix='_flux', mag_suffix='_mag'): # , flux_factors={}):
    #print('coeSED.load_phot')
    obj = extract_id(cat, id)
    fluxes = []
    fluxerrs = []
    for filt in filters:
        flux    = obj[filt+flux_suffix]
        fluxerr = obj[filt+flux_suffix+'err']
        mag     = obj[filt+mag_suffix]
        if mag < -90:
            flux = 0
            fluxerr = 1e30
        else:
            if input_flux_units == 'uJy':
                pass
            elif input_flux_units == 'nJy':
                flux    /= 1e3  # nJy -> uJy
                fluxerr /= 1e3  # nJy -> uJy

        fluxes.append(flux)
        fluxerr = np.min([fluxerr, 1e30])
        # otherwise Multinest outputs E+100 without the "E" (+100) and numpy crashes
        fluxerrs.append(fluxerr)

    #fluxes = np.array(fluxes)
    #print(fluxes, 'coeSED.load_phot fluxes')
    photometry = np.c_[fluxes, fluxerrs]
    return(photometry)

def load_phot_nJy(id, flux_suffix='_fluxnJy', mag_suffix='_mag'): # , flux_factors={}):
    obj = extract_id(cat, id)
    fluxes = []
    fluxerrs = []
    for filt in filters:
        flux    = obj[filt+flux_suffix]
        fluxerr = obj[filt+flux_suffix+'err']
        mag     = obj[filt+mag_suffix]
        if mag < -90:
            flux = 0
            fluxerr = 1e30
        else:
            flux    /= 1e3  # nJy -> uJy
            fluxerr /= 1e3  # nJy -> uJy

        fluxes.append(flux)
        fluxerr = np.min([fluxerr, 1e30])
        # otherwise Multinest outputs E+100 without the "E" (+100) and numpy crashes
        fluxerrs.append(fluxerr)

    #fluxes = np.array(fluxes)
    photometry = np.c_[fluxes, fluxerrs]
    return(photometry)
    
def hypotn(x):
    return np.sqrt(np.sum(x**2, axis=0))

def load_phot_sum(ids):
    all_fluxes = []
    all_fluxerrs = []
    for id in ids:
        fluxes, fluxerrs = load_phot(id).T
        all_fluxes.append(fluxes)
        all_fluxerrs.append(fluxerrs)

    all_fluxes = np.array(all_fluxes)
    all_fluxerrs = np.array(all_fluxerrs)
    
    flux_sum = np.sum(all_fluxes, axis=0)
    fluxerr_sum = hypotn(all_fluxerrs)
    photometry = np.c_[flux_sum, fluxerr_sum]
    return photometry

def plot_SED(ax, lams, fluxes, fluxerrs, mew=3, color='c', mec=None, ecolor=None):
    #ax.plot(lams, fluxes, 'o', ms=10, mew=3, mec='k', mfc='None', zorder=-10)    
    mec = mec or color
    ecolor = ecolor or 'k'
    for lam, flux, fluxerr in zip(lams, fluxes, fluxerrs):
        lam /= 1e3
        #fluxerr /= 1e4
        #print(lam, flux, fluxerr)
        if fluxerr < 1e30: # np.isfinite(fluxerr):
            ax.plot(lam, flux, 'o', ms=10, mew=mew, mec=mec, mfc='None', zorder=-10)
            if flux - fluxerr > 0:
                ax.errorbar(lam, flux, fluxerr, color=ecolor, mew=mew, capsize=2)
            else:
                yerr = np.array([flux, fluxerr])
                yerr.shape = 2, 1
                ax.errorbar([lam], [flux], yerr=yerr, color=ecolor, mew=mew, capsize=2)
        else:
            ax.plot(lam, 0, 'D', color='0.50', mfc='None')
            
    ax.axhline(0, c='0.50', lw=1, ls=':')
            
    ax.set_xlabel('Wavelength ($\mu$m)')
    ax.set_ylabel('Flux ($\mu$Jy)')
    #ax.set_ylabel('Flux (nJy)')

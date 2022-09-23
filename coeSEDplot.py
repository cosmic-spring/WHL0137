import numpy as np
import astropy.units as u
import math 

import os
from os.path import expanduser
home = expanduser("~")
import string

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
#%matplotlib notebook
# https://matplotlib.org/tutorials/introductory/customizing.html
try: plt.style.use('https://www.stsci.edu/~dcoe/matplotlibrc.txt') # 
except: plt.style.use(os.path.join(home, 'p', 'matplotlibrc.txt'))

mpl_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def extract_id(cat, id, idlabel='id'): # choose_object select_object
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

def roundint(x):
    return np.round(x).astype(int)
    
def between(lo, x, hi):
    return (x < lo) * (x < hi)
    
# Format exponential ticks to read as simple numbers unless too big / small
def fmtexp(x, pos):
    if between(1e-4, x, 1e4):
        s = '%g' % x
    else:
        p = int(np.log10(x))
        d = int(np.round(x / 10**p))
        #print d
        s = ''
        if d > 1:
            s += '%d$\\times$' % d
    
        s += '10${\\mathdefault{^{%d}}}$' % p
    return s

# Conversions for second axis
def AB2uJy(mAB):
    m = mAB * u.ABmag
    f = m.to(u.uJy)
    return f.value

def uJy2AB(F_uJy):
    f = F_uJy * u.uJy
    m = f.to(u.ABmag)
    return m.value

def AB2nJy(mAB):
    m = mAB * u.ABmag
    f = m.to(u.nJy)
    return f.value

def nJy2AB(F_nJy):
    f = F_nJy * u.nJy
    m = f.to(u.ABmag)
    return m.value

def slices_extent(x, y, dx, dy=0):
    dy = dy or dx
    xlo = roundint(x-dx)
    xhi = roundint(x+dx+1)
    ylo = roundint(y-dy)
    yhi = roundint(y+dy+1)
    xslice = slice(xlo, xhi)
    yslice = slice(ylo, yhi)
    slices = yslice, xslice
    extent = xlo, xhi, ylo, yhi
    return slices, extent

# Show color image and segment
def show_galaxy(id, figsize=(6,3)):
    global obj, xlo, xhi, ylo, yhi, cmap
    obj = extract_id(eazy_results, id)
    x = obj['x']
    y = obj['y']
    dx = dy = 100
    #xlo, xhi, ylo, yhi = x-dx, x+dx, y-dy, y+dy
    slices, extent = slices_extent(x, y, dx, dy)
    xlo, xhi, ylo, yhi = extent

    cmap = segm.make_cmap(seed=12345)
    cmap_colors = cmap.colors[:]
    cmap_colors[0]  = np.array([1, 1, 1])  # white background
    cmap_colors[id] = np.array([0, 0, 0])  # black selected object
    cmap = ListedColormap(cmap_colors)

    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(1,2,1)
    ax1.imshow(color_image[slices], extent=extent)
    plt.xlabel('')
    plt.ylabel('')
    ax1.axes.xaxis.set_visible(False)
    ax1.axes.yaxis.set_visible(False)

    ax2 = fig.add_subplot(1,2,2)
    plt.imshow(segm.data, cmap=cmap)
    plt.xlim(xlo,xhi)
    plt.ylim(ylo,yhi)
    plt.xlabel('')
    plt.ylabel('')
    ax2.axes.xaxis.set_visible(False)
    ax2.axes.yaxis.set_visible(False)
    fig.tight_layout()
    
    return fig

# coepipes.py
def plot_SED(id, fig, ax, plot_flux_units=u.uJy, plot_wavelength_units=u.um, rest_lam_interval=1000,
             plot_log_wavelength=True, plot_log_flux=False, lam_min=None, lam_max=None,
             obs_color='c', SED_color='brown', model_color='brown', model_lw=0.7, model_alpha=0.7):

    data = eazy_obj_data = eazy_run.show_fit(id, show_fnu=True, get_spec=True)
    z_phot = eazy_obj_data['z']  # photo-z
    #id = eazy_obj_data['id']
    #ix = eazy_obj_data['ix']
    #z = eazy_run.zgrid
    #pz = np.exp(eazy_run.lnp[ix]).flatten()

    #obj = extract_id(eazy_results, id)

    #eazy_flux_units = u.erg / u.s / u.cm**2 / u.AA

    data['flux_unit']
    data['wave_unit']
    
    wavelengths  = (data['pivot'] * data['wave_unit']).to(plot_wavelength_units)
    #fluxes       = (data['fobs']  * eazy_flux_units).to(plot_flux_units, u.spectral_density(wavelengths))
    fluxes       = (data['fobs']  * data['flux_unit']).to(plot_flux_units)
    flux_errs    = (data['efobs'] * data['flux_unit']).to(plot_flux_units)
    model_fluxes = (data['model'] * data['flux_unit']).to(plot_flux_units)
    #model_fluxes = data['model'] * u.uJy

    # Missing data
    missing  = (data['fobs']  < eazy_run.param['NOT_OBS_THRESHOLD']) 
    missing |= (data['efobs'] < 0)
    observed = ~missing

    snr_thresh = 1
    #snr_thresh = -1e30  # NEVER!
    SNR = data['fobs'] / data['efobs']
    #detected     = (~missing) & (SNR >  snr_thresh)
    #not_detected = (~missing) & (SNR <= snr_thresh)
    snr0         = (~missing) & (SNR < 1)
    snr1         = (~missing) & (SNR > 1) & (SNR < 2)
    snr2         = (~missing) & (SNR > 2) & (SNR < 3)
    snr3         = (~missing) & (SNR > 3)
    detected = snr1 + snr2 + snr3
    not_detected = snr0
    
    #print('snr0', wavelengths[snr0])
    #print('snr1', wavelengths[snr1])
    #print('snr2', wavelengths[snr2])
    #print('snr3', wavelengths[snr3])

    # Plot y scale
    goodflux = fluxes > flux_errs

    flux_min = np.min(fluxes[detected].value)
    flux_max = np.max(fluxes[goodflux].value)

    model_flux_max = np.max(model_fluxes.value)
    model_flux_min = np.min(model_fluxes.value)

    if np.any(not_detected):
        upper_limits = flux_errs[not_detected].value
        upper_limit_min = np.min(upper_limits)
        flux_min = np.min([flux_min, upper_limit_min])

    flux_max = np.max([flux_max, model_flux_max])
    #if plot_full_spec:
    #    flux_max = np.max([flux_max, spec_max])

    if plot_log_flux:
        yhi = flux_max * 1
        ylo = flux_min / 2
    else: # plot linear
        flux_min = np.min([0, flux_min])
        flux_margin = 0.1 * (flux_max - flux_min)
        yhi = flux_max + flux_margin * 2
        ylo = -0.1 * flux_max
        #ylo = flux_min - flux_margin
        #print("flux_min", flux_min, ylo)

    #fig, ax = plt.subplots(figsize=(8, 5))

    # Observed flux error bars
    plt.errorbar(wavelengths[observed], fluxes[observed], flux_errs[observed], ms=8, marker='o', mfc=obs_color, c='k', lw=3, alpha=0.5, ls='none')
    #plt.errorbar(wavelengths, fluxes, flux_errs, ms=8, marker='o', mfc=obs_color, c='k', lw=3, alpha=0.5, ls='none')

    #plt.errorbar(wavelengths[snr3], fluxes[snr3], flux_errs[snr3], ms=8, marker='o', mfc=obs_color, c='k', lw=3, alpha=0.5, ls='none')
    #plt.errorbar(wavelengths[snr2], fluxes[snr2], flux_errs[snr2], ms=8, marker='o', mfc=obs_color, c='k', lw=3, alpha=0.4, ls='none')
    #plt.errorbar(wavelengths[snr1], fluxes[snr1], flux_errs[snr1], ms=8, marker='o', mfc=obs_color, c='k', lw=3, alpha=0.3, ls='none')
    #plt.errorbar(wavelengths[snr0], fluxes[snr0], flux_errs[snr0], ms=8, marker='o', mfc=obs_color, c='k', lw=3, alpha=0.2, ls='none')

    # Observed fluxes
    #plt.plot(wavelengths[detected], fluxes[detected],       'o', color=obs_color,   ms=8,  
    #         label='Input fluxes', zorder=10)#, scaley=False)

    # Model fluxes
    plt.plot(wavelengths[observed], model_fluxes[observed], 's', color=model_color, ms=10, mfc='None', label='Model fluxes', zorder=10)#, scaley=False)
    #plt.plot(wavelengths, model_fluxes, 's', color=model_color, ms=10, mfc='None', label='Model fluxes', zorder=10)#, scaley=False)
    #plt.plot(wavelengths[detected], model_fluxes[detected], 's', color=model_color, ms=10, mfc='None', 
    #         label='Model fluxes', zorder=10)#, scaley=False)

    #plt.semilogx()
    if plot_log_flux:
        plt.semilogy()
    else:
        plt.axhline(0, c='0.50', lw=1, ls=':')
        #ax.set_ylim(ylo, yhi)

    if plot_log_wavelength:
        ax.semilogx()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter("%g"))
                
    #plt.margins(1)
    
    #plt.xlim(1,5)
        
    ax.autoscale(False)
    
    # Non detections
    #plt.plot(wavelengths[not_detected], flux_errs[not_detected],
    #         ms=8, marker='v', mfc=obs_color, c='k', lw=3, alpha=0.5, ls='none')

    # SED
    SED_wavelength = (data['templz'] * u.AA).to(u.um)
    SED_flux = data['templf'] * u.uJy
    plt.plot(SED_wavelength, SED_flux, '-', color=model_color, lw=model_lw, alpha=model_alpha, zorder=-10)

    #imin = SED_wavelength.value.searchsorted(xmin)
    #imax = SED_wavelength.value.searchsorted(xmax)
    #plt.plot(SED_wavelength[imin:imax], SED_flux[imin:imax], '-', color=model_color, lw=model_lw, alpha=model_alpha,
    #plt.plot(SED_wavelength, SED_flux, '-', scalex=False, scaley=False)

    xmin, xmax = plt.xlim()
    xmin = lam_min or xmin
    xmax = lam_max or xmax
    plt.xlim(xmin, xmax)
    plt.ylim(ylo, yhi)

    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('Flux ($\mu$Jy)')

    #ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmtexp))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%g"))

    if plot_log_flux:
        secax = ax.secondary_yaxis('right', functions=(uJy2AB, AB2uJy))
        #secax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
        secax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        secax.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%g"))
        secax.set_ylabel('Magnitude (AB)')
    else:
        secax = add_magnitude_axis(ax, plot_flux_units)

    topax = add_rest_wavelength_axis(ax, z_phot, rest_lam_interval=rest_lam_interval)
    
    
def add_rest_wavelength_axis(ax, z, rest_units=1000*u.AA, rest_lam_interval=1000, obs_lam_interval=1):
    factor = u.um.to(rest_units)
    def bottom2top(x): return x / (1+z) * factor
    def top2bottom(x): return x * (1+z) / factor
    ax.xaxis.set_major_locator(ticker.MultipleLocator(obs_lam_interval))  # Original bottom axis (microns)
    topax = ax.secondary_xaxis('top', functions=(bottom2top, top2bottom))
    tick_locator = (rest_lam_interval*u.AA / rest_units).value
    topax.xaxis.set_major_locator(ticker.MultipleLocator(tick_locator))  # New top axis
    xticks = topax.get_xticks()
    xmin1, xmax1 = ax.get_xlim()
    xmax2 = bottom2top(xmax1)
    if xmax2 > 25:
        tick_locator *= (xmax2 // 25 + 1)
        topax.xaxis.set_major_locator(ticker.MultipleLocator(tick_locator))  # New top axis
    #rest_label = 'Rest Wavelength ($1000 \AA$)'
    rest_label = 'Rest Wavelength ($%d \\rm{\AA}$)' % rest_units.value
    if rest_units.value == 1:
        rest_label = 'Rest Wavelength ($\\rm{\AA}$)' % rest_units.value
    topax.set_xlabel(rest_label, fontsize=14)
    topax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    #topax.xaxis.set_minor_formatter(ticker.FormatStrFormatter("%g"))
    return topax
    
# Add AB magnitudes as secondary x-axis at right
# https://matplotlib.org/gallery/subplots_axes_and_figures/secondary_axis.html#sphx-glr-gallery-subplots-axes-and-figures-secondary-axis-py

# sedplotnorm.py
def add_magnitude_axis(ax, flux_units=u.nJy, plothuge=True):
    ylo, yhi = plt.ylim() * flux_units
    maghi = yhi.to(u.ABmag).value
    ytx1 = np.ceil(maghi * 10) / 10.  # 24.101 -> 24.2
    ytx2 = np.ceil(maghi)  # 24.1 -> 25

    fpart = ytx1 - int(ytx1)  # 0.2
    if np.isclose(fpart, 0) or np.isclose(fpart, 0.9):
        ytx1 = []
    elif np.isclose(fpart, 0.1) or np.isclose(fpart, 0.2):
        ytx1 = np.array([ytx1, ytx2-0.7, ytx2-0.5, ytx2-0.3])  # 24.1, 24.3, 24.5, 24.7
    elif np.isclose(fpart, 0.3) or np.isclose(fpart, 0.4):
        ytx1 = np.array([ytx1, ytx2-0.5, ytx2-0.3])  # 24.3, 24.5, 24.7
    elif np.isclose(fpart, 0.5):
        ytx1 = np.array([ytx1, ytx2-0.3])  # 24.5, 24.7
    elif np.isclose(fpart, 0.6):
        ytx1 = np.array([ytx1, ytx2-0.2])  # 24.6, 24.8

    if isinstance(ytx1, float):
        ytx1 = np.array([ytx1])

    if plothuge:
        ytx3 = ytx2 + np.array([0, 0.2, 0.5, 1, 2])
    else:
        ytx3 = ytx2 + np.array([0, 0.2, 0.5, 1, 1.5, 2, 3])

    ytx2 = np.array([ytx2])
    ytx = np.concatenate([ytx1, ytx3])
    yts = ['%g' % mag for mag in ytx]

    ytx = (ytx * u.ABmag).to(flux_units).value
    
    ax2 = ax.twinx()
    ax.yaxis.set_label_position('left')
    ax2.yaxis.set_label_position('right')
    ax.yaxis.tick_left()
    ax2.yaxis.tick_right()

    #plt.yticks(ytx, yts)
    ax2.set_yticks(ytx)
    ax2.set_yticklabels(yts)
    ax2.set_ylabel('Magnitude (AB)')
    ax2.set_ylim(ylo.value, yhi.value)
    ax2.set_zorder(-100)  # interactive cursor will output left axis ax
    return ax2

# Plot log(1+z)
def z2log(z):
    return np.log10(1+z)

def log2z(log1z):
    return 10 ** log1z - 1
    
#from memory_profiler import profile

#@profile
def show_galaxy_properties(id, hist_color='C0', figsize=(12, 6), save=False, alt_id=None, dpi=100,
    lam_min=None, lam_max=None):
    data = eazy_obj_data = eazy_run.show_fit(id, show_fnu=True,get_spec=True)
    #id = eazy_obj_data['id']
    ix = eazy_obj_data['ix']
    z = eazy_run.zgrid
    pz = np.exp(eazy_run.lnp[ix]).flatten()

    obj = extract_id(eazy_results, id)

    # Load data
    x = obj['x']
    y = obj['y']
    dx = dy = 100
    #xlo, xhi, ylo, yhi = x-dx, x+dx, y-dy, y+dy
    slices, extent = slices_extent(x, y, dx, dy)
    xlo, xhi, ylo, yhi = extent
    cmap = segm.make_cmap(seed=12345)
    cmap_colors = cmap.colors[:]
    cmap_colors[0]  = np.array([1, 1, 1])  # white background
    cmap_colors[id] = np.array([0, 0, 0])  # black selected object
    cmap = ListedColormap(cmap_colors)

    # Create figure and gridspec
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 6, wspace=0)  # nrows, ncolumns

    # Add text to figure with ID and redshift (RA, Dec)
    coord_str = 'id #%d' % id
    if alt_id: coord_str += ' (%s)' % alt_id
    coord_str += '\n'
    coord_str += 'z = %.2f' %  obj['z_phot']
    coord_str += ' [%.2f $-$ %.2f]' %  (obj['z025'], obj['z975'])  #  (95%% CL)

    fig.text(0.01, 0.99, coord_str, va='top')
    
    
    # Add text to figure with coordinates (RA, Dec)
    space = '$~~~$'
    if obj['dec'] < 0: dec_prefix = ' $-$'
    else: dec_prefix = space

    coord_str = 'RA   =' + space
    coord_str += '%10.7f' %  obj['ra']
    coord_str += ' = ' + space
    coord_str += RAdeg2hms(obj['ra'],  ':', precision=2)
    coord_str += '\n'
    coord_str += 'Dec ='
    if np.abs(obj['dec']) < 10:
        coord_str += '  ' # space
    coord_str += dec_prefix + '%.7f' %  np.abs(obj['dec'])
    coord_str += ' =' 
    if np.abs(obj['dec']) < 10:
        coord_str += '' # space
    coord_str += dec_prefix + deg2hms(  np.abs(obj['dec']), ':', precision=1)

    fig.text(0.18, 0.99, coord_str, va='top')
    
    # Color image
    #ax1 = fig.add_subplot(3,6,1)  # nrows, ncolumns, index
    ax1 = fig.add_subplot(gs[0,0])
    #ax1.imshow(color_image)
    ax1.imshow(color_image[slices], extent=extent)
    #plt.xlim(xlo, xhi)
    #plt.ylim(ylo, yhi)
    ax1.axes.xaxis.set_visible(False)
    ax1.axes.yaxis.set_visible(False)

    # Segmentation map
    #ax2 = fig.add_subplot(3,6,2)  # nrows, ncolumns, index
    ax2 = fig.add_subplot(gs[0,1])
    plt.imshow(segm.data, cmap=cmap, interpolation='none')
    plt.xlim(xlo, xhi)
    plt.ylim(ylo, yhi)
    ax2.axes.xaxis.set_visible(False)
    ax2.axes.yaxis.set_visible(False)

    # Redshift
    gs3 = fig.add_gridspec(100, 3, wspace=0.2)  # nrows, ncolumns
    ax3 = fig.add_subplot(gs3[45:50,0])
    #plt.hist(fit.posterior.samples['redshift'], ec=colors[0], lw=3)
    #plt.plot(z, pz, color=colors[0], lw=3)
    #plt.fill_between(z, pz, pz*0, color=hist_color)
    plt.fill_between(z2log(z), pz, pz*0, color=hist_color)

    ztx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20]
    zmax = np.max(z)
    #ztx_minor = 0.2, 0.4, 0.6, 0.8, 9, 11, 12, 13
    ztx_minor = np.arange(0.2,1,0.2)
    ztx_minor = set(ztx_minor) | set(np.arange(zmax)).difference(ztx)
    ztx_minor = list(ztx_minor)
    plt.xlim(0, z2log(zmax))
    #ax3.xaxis.tick_bottom()
    plt.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=False, labelleft=False)
    sec_xax = ax3.secondary_xaxis('top', functions=(log2z, z2log))
    sec_xax.set_xticks(ztx)
    sec_xax.set_xticks(ztx_minor, minor=True)
    sec_xax.set_xlabel('Redshift', fontsize=16)        


    # specific SFR (sSFR)
    ax3b = fig.add_subplot(gs3[70:75,0])
    ssfr_percentiles = extract_id(eazy_results, id)['ssfr_p']
    log_ssfr_Gyr_p = np.log10(ssfr_percentiles) + 9
    #plt.hist(log_sSFR, ec=colors[0], lw=3)
    #for log_ssfr_Gyr in log_ssfr_Gyr_p:
    #    plt.axvline(log_ssfr_Gyr, color=colors[0], lw=1)
    #print(log_ssfr_Gyr_p[0], log_ssfr_Gyr_p[-1])
    if np.isnan(log_ssfr_Gyr_p[-1]):
        log_ssfr_Gyr_p[0] = log_ssfr_Gyr_p[-1] = -4
    log_ssfr_min, log_ssfr_max = -4, 1  # Fit in plot range
    log_ssfr_Gyr_p = np.clip(log_ssfr_Gyr_p, log_ssfr_min+0.03, log_ssfr_max-0.03)
    #print(log_ssfr_Gyr_p)
    plt.fill_betweenx([0,1], log_ssfr_Gyr_p[0], log_ssfr_Gyr_p[-1], color=hist_color, linewidth=3)
    plt.xlim(log_ssfr_min, log_ssfr_max)
    plt.ylim(0,1)
    ax3b.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax3b.xaxis.tick_top()
    ax3b.axes.yaxis.set_visible(False)
    #ax3b.set_title('log sSFR ( / Myr)', fontsize=16)
    ax3b.set_title('log sSFR (Gyr$^{-1}$)', fontsize=16)

    # Dust
    ax3c = fig.add_subplot(gs3[95:100, 0])
    dust_percentiles = extract_id(eazy_results, id)['Av_p']
    max_dust = 5
    dust_percentiles = np.clip(dust_percentiles, 0, max_dust)
    plt.fill_betweenx([0,1], dust_percentiles[0], dust_percentiles[-1], color=hist_color, linewidth=3)
    #for dust in dust_percentiles:
    #    plt.axvline(dust, color='c', lw=1)
    #plt.xlim(dust["Av"][0], dust["Av"][1])
    plt.xlim(0, max_dust)
    plt.ylim(0, 1)
    ax3c.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax3c.xaxis.tick_top()
    ax3c.axes.yaxis.set_visible(False)
    ax3c.set_title('dust $A_V$ (mag)', fontsize=16)

    # SED
    gs4 = fig.add_gridspec(1, 3)  # nrows, ncolumns
    ax4 = fig.add_subplot(gs4[0,1:])
    #plt1, fig1, ax1, secax = 
    plot_SED(id, fig, ax4, lam_min=lam_min, lam_max=lam_max)

    if save:
        fig.savefig('plots/eazy_%d.png' % id, dpi=dpi)
        plt.close(fig)
    
    return fig
    


# previously h2hms
def deg2hms(x, format=(), precision=3, lead0=True, RA=False):
    """
    CONVERTS decimal degrees/hours to degrees/hours : minutes : seconds
    deg2hms(13.52340987)
    deg2hms(13.52340987, ':')
    deg2hms(13.52340987, 'hms')
    deg2hms(13.52340987, 'dms')
    deg2hms(13.52340987, 'dms', 1)
    """
    if RA:
        x = x / 15  # 360 degrees -> 24 hours
    f, i = math.modf(x)
    i = int(i)
    m = 60 * f
    s, m = math.modf(m)
    m = int(m)
    s = 60 * s
    if m:
        s = abs(s)
    if i:
        m = abs(m)
    if type(format) == str:
        if precision == None:
            s = '%f' % s
        else:
            if lead0:
                nd = precision + 3
                fmt = '%%0%d.%df' % (nd, precision)
            else:
                fmt = '%%.%df' % precision
            s = fmt % s
        if s.startswith('60'):  # rounded up
            s = '0'
            m = m + 1
        if lead0:
            m = '%02d' % m
        else:
            m = '%d' % m
        if m == '60':  # rounded up
            m = '0'
            i += 1
        if lead0:
            i = '%02d' % i
        else:
            i = '%d' % i
        ims = (i,m,s)
        if len(format) == 1:
            out = i + format + m + format + s
        elif len(format) == 3:
            out = i+format[0] + m+format[1] + s+format[2]
    else:
        out = (i, m, s)
    return out

h2hms = deg2hms

def RAdeg2hms(x, format=(), precision=3, lead0=True):
    return deg2hms(x, format, precision, lead0, RA=True)


    # Redshift
    gs3 = fig.add_gridspec(100, 3, wspace=0.2)  # nrows, ncolumns
    ax3 = fig.add_subplot(gs3[45:50,0])
    #plt.hist(fit.posterior.samples['redshift'], ec=colors[0], lw=3)
    #plt.plot(z, pz, color=colors[0], lw=3)
    plt.fill_between(z, pz, pz*0, color=hist_color)
    ax3.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax3.xaxis.tick_top()
    ax3.axes.yaxis.set_visible(False)
    ax3.set_title('Redshift', fontsize=16)

    # specific SFR (sSFR)
    ax3b = fig.add_subplot(gs3[70:75,0])
    ssfr_percentiles = extract_id(eazy_results, id)['ssfr_p']
    log_ssfr_Gyr_p = np.log10(ssfr_percentiles) + 9
    #plt.hist(log_sSFR, ec=colors[0], lw=3)
    #for log_ssfr_Gyr in log_ssfr_Gyr_p:
    #    plt.axvline(log_ssfr_Gyr, color=colors[0], lw=1)
    #print(log_ssfr_Gyr_p[0], log_ssfr_Gyr_p[-1])
    if np.isnan(log_ssfr_Gyr_p[-1]):
        log_ssfr_Gyr_p[0] = log_ssfr_Gyr_p[-1] = -4
    plt.fill_betweenx([0,1], log_ssfr_Gyr_p[0], log_ssfr_Gyr_p[-1], color=hist_color)
    plt.xlim(-4, +1)
    plt.ylim(0,1)
    ax3b.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax3b.xaxis.tick_top()
    ax3b.axes.yaxis.set_visible(False)
    #ax3b.set_title('log sSFR ( / Myr)', fontsize=16)
    ax3b.set_title('log sSFR (Gyr$^{-1}$)', fontsize=16)

    # Dust
    ax3c = fig.add_subplot(gs3[95:100, 0])
    dust_percentiles = extract_id(eazy_results, id)['Av_p']
    plt.fill_betweenx([0,1], dust_percentiles[0], dust_percentiles[-1], color=hist_color)
    #for dust in dust_percentiles:
    #    plt.axvline(dust, color='c', lw=1)
    #plt.xlim(dust["Av"][0], dust["Av"][1])
    plt.xlim(0, 3)
    plt.ylim(0, 1)
    ax3c.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax3c.xaxis.tick_top()
    ax3c.axes.yaxis.set_visible(False)
    ax3c.set_title('dust $A_V$ (mag)', fontsize=16)

    # SED
    gs4 = fig.add_gridspec(1, 3)  # nrows, ncolumns
    ax4 = fig.add_subplot(gs4[0,1:])
    #plt1, fig1, ax1, secax = 
    plot_SED(id, fig, ax4)

    if save:
        fig.savefig('EAZY_%d.png' % id)
    
    return fig
    
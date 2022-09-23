# run_id = "tan-efolds_bpass_dust3-eta2_logU41_Zall"

from copy import deepcopy

def set_fit_instructions(run_id):
    dust = {}                         # Dust component
    
    if 'SMC' in run_id:
        dust["type"] = "Salim"         # Define the shape of the attenuation curve    
        dust["delta"] = -0.45          # SMC
        dust["B"] = 0                  # SMC
    else:
        dust["type"] = "Calzetti"         # Define the shape of the attenuation curve
    
    if 'dust5' in run_id:
        dust["Av"] = 0, 5                # magnitudes
    elif 'dust3' in run_id:
        dust["Av"] = 0, 3                # magnitudes
    else:
        dust["Av"] = 0, 1                # magnitudes

    if '1dust' in run_id:
        dust["eta"] = 1                  # Extra dust for young stars: multiplies Av
    elif 'eta2' in run_id:
        dust["eta"] = 2                  # Extra dust for young stars: multiplies Av
    elif 'eta1-3' in run_id:
        dust["eta"] = 1, 3               # Extra dust for young stars: multiplies Av
    else:
        dust["eta"] = 1, 10              # Extra dust for young stars: multiplies Av

    nebular = {}                      # Nebular emission component
    #nebular["logU"] = -3.             # log_10(ionization parameter)
    #if 'logU-1' in run_id:  # grid only runs from -4 to -2 (config.py)
    #    nebular["logU"] = (-4., -1.)             # log_10(ionization parameter)
    if 'logU41' in run_id:
        nebular["logU"] = (-4., -1.)             # log_10(ionization parameter)
    else:
        nebular["logU"] = (-4., -2.)             # log_10(ionization parameter)

    fit_instructions = {}                   # The model components dictionary
    if 'zall' in run_id:
        fit_instructions["redshift"] = 0,12      # Observed redshift  
    else:
        fit_instructions["redshift"] = 4,11      # Observed redshift  

    sfh = {}
    if 'Z02' in run_id:
        sfh['metallicity'] = 0.2  # Z / Zsun
    else:
        sfh['metallicity_prior'] = 'log_10'
        if 'Zall' in run_id:
            sfh['metallicity'] = 0.005, 5  # Z / Zsun; allowed range of mdoels
        else:
            sfh['metallicity'] = 0.005, 2.5
        
    sfh['massformed'] = 6, 12
    
    if 'delayed' in run_id:
        sfh['age'] = 0.001, 2         # Gyr: Time since star formation began
        sfh['tau'] = -1, -0.001        # Gyr: Timescale of change (negative rising; positive falling)
        fit_instructions["delayed"] = sfh
    elif 'const-burst' in run_id:
        sfh_const = deepcopy(sfh)
        sfh_const['age_min'] = 0.01, 2  # 10 Myr - 2 Gyr
        sfh_const['age_max'] = 0.01, 2  # 10 Myr - 2 Gyr
        #const['age_max'] = 0, 'age_of_universe'
        fit_instructions["constant"] = sfh_const

        sfh_burst = deepcopy(sfh)
        sfh_burst['age'] = 0.01, 0.1  # 10 - 100 Myr
        fit_instructions["burst"] = sfh_burst
    elif 'const' in run_id:
        sfh['age_min'] = 0  # doesn't stop
        #sfh['age_max'] = 0, 'age_of_universe'
        sfh['age_max'] = 0, 2
        fit_instructions["constant"] = sfh
    elif 'declining' in run_id:
        sfh['age'] = 0.001, 2         # Gyr: Time since star formation began
        sfh['tau'] = 0, 1       # Gyr: Timescale of change (negative rising; positive falling)
        fit_instructions["delayed"] = sfh
    elif 'efolds' in run_id:
        sfh['age'] = 0.001, 2         # Gyr: Time since star formation began
        #sfh['tau'] = -1, -0.001        # Gyr: Timescale of change (negative rising; positive falling)
        #sfh['efolds'] = -500, 500
        sfh['tanefolds'] = -1, 1
        fit_instructions["exponential"] = sfh

    fit_instructions["dust"] = dust
    fit_instructions["t_bc"] = 0.01         # Lifetime of birth clouds (Gyr)
    fit_instructions["nebular"] = nebular
    
    return fit_instructions


import astropy
from astropy.io import ascii
import astropy.units as u
import os
import numpy as np

import matplotlib.pyplot as plt
# https://matplotlib.org/tutorials/introductory/customizing.html
#plt.style.use('/Users/dcoe/p/matplotlibrc.txt')
plt.style.use('https://www.stsci.edu/~dcoe/matplotlibrc.txt')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
import matplotlib.ticker as ticker

import bagpipes as pipes


filt_root_dir = '/Users/dcoe/bpz2/'

def get_full_filter_name(filt):
    filt_subdir = 'FILTER'
    filt_dir = os.path.join(filt_root_dir, filt_subdir)
    filt = filt.upper()
    
    filt_names = []
    filt_names.append('JWST_NIRCAM_%s.res')
    filt_names.append('HST_WFC3_IR_%s.res')
    filt_names.append('HST_ACS_WFC_%s.res')
    filt_names.append('HST_WFC3_UVIS_%s.res')
    
    for filt_name in filt_names:
        filt_file_basename = filt_name % filt
        filt_file = os.path.join(filt_dir, filt_file_basename)
        if os.path.exists(filt_file):
            return os.path.join(filt_subdir, filt_file_basename)

    return None

#from coeETC import get_wavelength_color
#import coeETC


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
        
# coeSED
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
    return topax

# coeETC
def get_wavelength_color(lam, lam_min=0.3, lam_max=5, cmap='nipy_spectral'):
    cm = plt.get_cmap(cmap)
    cNorm  = colors.Normalize(vmin=lam_min, vmax=lam_max)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    color = scalarMap.to_rgba(lam)
    return color

# Plot SEDs:
#  "best" single best
#  "age range" 1-sigma range of mass-weighted ages

#def plot_SED(fit, cat, id, choice='best', age_label='mass_weighted_age', flux_units=u.nJy, lam_max=5, plot_full_spec=False, save_plot='', plot_title=None, overwrite=False):

def plot_SED(fit, id, magnification=1, choice='best', age_label='mass_weighted_age', flux_units=u.nJy, fmin=None, fmax=None, ferrmax=1e29, lam_min=0, lam_max=5, lam_log=False, plot_full_spec=False, fontsize=14, save_plot='', plot_title=None, plot_params=True, close_plot=False, overwrite=False, show_plot=True, create_fig=True, input_figax=[],
obs_color='c', obs_colors=[], obs_ecolor='c', SED_color='brown', model_color='brown', text_color='brown', model_lw=0.7,
more_text=True, less_text=False, no_text=False, show_full_SED_flux_limits=False):
# obs_color=None, obs_colors=[], obs_ecolor='k', SED_color='navajowhite', model_color='darkorange', text_color='chocolate', model_lw=0.5):
    #print('plot_title = ', plot_title)
    if save_plot:
        if save_plot == True:
            # fit.fname 'pipes/posterior/delayed-tau_zphot/469_'
            plot_file = fit.fname.replace('posterior/', 'plots/')
            plot_file = os.path.dirname(plot_file)
            plot_file = os.path.join(plot_file, 'bagpipes_sed_%d.png' % id)
            #plot_file = os.path.join(plot_file, 'SED_%d.png' % id)
        else:
            plot_file = save_plot
        if os.path.exists(plot_file) and not show_plot:
            if not overwrite:
                print(plot_file, 'EXISTS')
                return None, None, None, None
    
    age_label_form, sfh_function = extract_age_sfh(fit)
    n = len(fit.posterior.samples[age_label])
    #if create_fig:
    #    fig, ax1 = plt.subplots(1, figsize=(8,6), dpi=100)
    #else:
    #    fig = plt.gcf()
    #    ax1 = plt.gca()
    if len(input_figax):
        fig, ax1 = input_figax
    else:
        fig, ax1 = plt.subplots(1, figsize=(8,6), dpi=100)

    if choice == 'best':
        imodel = np.argmin(fit.posterior.samples['chisq_phot'])
        imodels = [imodel]
        SED_colors = [SED_color]
        SED_labels = [None]
        obs_color = obs_color or 'c'
    elif choice == 'age range':
        ages = fit.posterior.samples[age_label]
        imodelsort = np.argsort(ages)
        percentiles = np.array([16, 50, 84])
        isorts = percentiles/100. * len(ages)
        isorts = isorts.astype(int)
        imodels = imodelsort[isorts]

        #SED_colors = ['c', 'navajowhite', 'r']
        #SED_colors = ['c', 'y', 'r']
        SED_colors = ['b', (0.6,0.6,0), 'r']
        #obs_color = 'lime'
        #obs_color = 'g'
        #obs_color = 0,0.7,0
        #obs_color = 'gray'
        #obs_color = 'None'
        obs_color = 'k'

        SED_ages = fit.posterior.samples[age_label][imodels]
        SED_labels = []
        SED_labels.append('%3d Myr (1-sigma younger)'   % (1000*fit.posterior.samples[age_label][imodels[0]]))
        SED_labels.append('%3d Myr (median age mass-weighted)' % (1000*fit.posterior.samples[age_label][imodels[1]]))
        SED_labels.append('%3d Myr (1-sigma older)'     % (1000*fit.posterior.samples[age_label][imodels[2]]))
    elif choice == 'age young old':
        ages = fit.posterior.samples[age_label]
        imodelsort = np.argsort(ages)
        percentiles = np.array([16, 84])
        isorts = percentiles/100. * len(ages)
        isorts = isorts.astype(int)
        imodels = imodelsort[isorts]

        SED_colors = ['c', 'r']
        #SED_colors = ['c', 'y', 'r']
        obs_color = 'y'

        SED_ages = fit.posterior.samples[age_label][imodels]
        SED_labels = []
        SED_labels.append('%3d Myr (1-sigma young mass-weighted age)'   % (1000*fit.posterior.samples[age_label][imodels[0]]))
        SED_labels.append('%3d Myr (1-sigma old mass-weighted age)'     % (1000*fit.posterior.samples[age_label][imodels[2]]))
    elif choice == 'best ages':
        ages = fit.posterior.samples[age_label]
        median_age = np.median(ages)
        iall = np.arange(n)
        iyoung = [i for i in iall if ages[i] <  median_age]
        iold   = [i for i in iall if ages[i] >= median_age]
        ibest_young = iyoung[np.argmin(fit.posterior.samples['chisq_phot'][iyoung])]
        ibest_old   = iold[np.argmin(fit.posterior.samples['chisq_phot'][iold])]
        imodels = ibest_young, ibest_old
        print(imodels)
        
        SED_colors = ['c', 'r']
        obs_color = 'y'
        
        #SED_ages = fit.posterior.samples[age_label][imodels]
        SED_labels = []
        SED_labels.append('%3d Myr (best young mass-weighted)'   % (1000*fit.posterior.samples[age_label][imodels[0]]))
        SED_labels.append('%3d Myr (best old mass-weighted)'     % (1000*fit.posterior.samples[age_label][imodels[1]]))
    else:
        imodel = choice
        choice = 'one'
        imodels = [imodel]
        SED_colors = [SED_color]
        SED_labels = [None]
        obs_color = obs_color or 'c'

    #print(obs_color, 'obs_color')

    fit_keys = fit.posterior.samples.keys()
    fixed_redshift = 'redshift' not in fit_keys
    if fixed_redshift:
        redshift = fit.fit_instructions["redshift"]
    
    model_spectra = fit.posterior.samples['spectrum_full']
    spec_max = 0
    for iSED, imodel in enumerate(imodels):
        if not fixed_redshift:
            redshift = fit.posterior.samples['redshift'][imodel]

        model_lams = fit.posterior.model_galaxy.wavelengths * (1.+redshift)
        model_lams = (model_lams * u.AA).to(u.micron)
        ilam_max = model_lams.value.searchsorted(lam_max)

        color = SED_colors[iSED]
        spec_flam = model_spectra[imodel]
        #spec_flam /= 1e3; print("Dividing flux by 1000")
        #spec_flam *= 1e9; print("Multiplying flux by 1000")
        # If you alter spec_flam, it alters fit in place
        spec_fnu = (spec_flam * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(model_lams))
        #spec_fnu /= 1e3; print("Dividing flux by 1000")
        #plt.plot(model_lams[:ilam_max], spec_fnu[:ilam_max], lw=0.5, color=color, alpha=0.7, label=SED_labels[iSED])
        ax1.plot(model_lams[:ilam_max], spec_fnu[:ilam_max], lw=model_lw, color=color, alpha=0.7)
        #print(np.max(spec_flam[:ilam_max]), 'MAX FLUX spec_flam')
        #print(np.max(spec_fnu[:ilam_max]), 'MAX FLUX spec_fnu')
        if choice not in ('best', 'one'):
            ax1.plot(model_lams[:ilam_max].value*1e30, spec_fnu[:ilam_max].value*1e30, lw=1, color=color, alpha=1, label=SED_labels[iSED])
        if between(0, plot_full_spec, 1):
            spec_max1 = np.percentile(spec_fnu[:ilam_max].value, 100*plot_full_spec)
        else:
            spec_max1 = np.max(spec_fnu[:ilam_max].value)
        spec_max = np.max([spec_max, spec_max1])

    lam = fit.galaxy.filter_set.eff_wavs
    lam = (lam * u.AA).to(u.micron)

    modelflams = fit.posterior.samples['photometry']
    #modelflams /= 1e3; print("Dividing flux by 1000")
    #modelflams *= 1e9; print("Multiplying flux by 1000")
    modelfnus = (modelflams * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(lam))
    #modelfnus /= 1e3; print("Dividing flux by 1000")

    # Observed photometry (excluding unobserved filters)
    obslam  = fit.galaxy.filter_set.eff_wavs

    lam, flam, flamerr = fit.galaxy.photometry.T
    #flam    = np.array(flam) / 1e3; print("Dividing flux by 1000")
    #flamerr = np.array(flamerr) / 1e3; print("Dividing flux by 1000")
    lam    = (lam * u.AA).to(u.micron)
    fnu    = (flam    * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(lam))
    fnuerr = (flamerr * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(lam))
    #print(fnuerr)
    observed = fnuerr.value < ferrmax
    #print(ferrmax, np.sum(observed), len(observed))
    #print(flamerr)
    lam, flam, flamerr = fit.galaxy.photometry[observed].T
    #flam    = np.array(flam) / 1e3; print("Dividing flux by 1000")
    #flamerr = np.array(flamerr) / 1e3; print("Dividing flux by 1000")
    modelflams = modelflams[:,observed]
    modelfnus  = modelfnus[:,observed]
    lam    = (lam * u.AA).to(u.micron)
    fnu    = (flam    * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(lam))
    #fnu    = (0.001 * flam    * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(lam))
    #print("Dividing flux by 1000")
    fnuerr = (flamerr * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(lam))
    #fnuerr = (0.001 * flamerr * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(lam))
    #print("Dividing flux by 1000")
    #fnu    = np.array(fnu)    / 1e3; print("Dividing flux by 1000")
    #fnuerr = np.array(fnuerr) / 1e3; print("Dividing flux by 1000")
    
    #ax1.errorbar(lam.value, fnu.value, yerr=fnuerr.value, fmt='.', marker='None', color='k', zorder=10, capsize=6, mew=2.5)
    ax1.errorbar(lam.value, fnu.value, yerr=fnuerr.value, fmt='.', color='k', zorder=10, capsize=6, mew=2.5)
    if len(obs_colors):
        clam_min, clam_max = obs_colors
        for i in range(len(lam)):
            obs_color = coeETC.get_wavelength_color(lam[i], lam_min=clam_min, lam_max=clam_max)
            ax1.plot(lam[i], fnu[i], 'o', color=obs_color, markersize=8, mec=obs_ecolor)
    else:
        ax1.plot(lam, fnu, 'o', color=obs_color, markersize=8, mec=obs_ecolor)
        
    #plt.plot(lam, fnu, 'o', color='k', markersize=8, alpha=1)
    #plt.plot(lam, fnu, 'o', color='None', markersize=8, mec='k', mew=2)

    fnu_min = flux_min = np.min(fnu.value)
    #flux_max = np.max(fnu.value)
    goodflux = fnu > fnuerr
    fnu_max = flux_max = np.max(fnu[goodflux].value)

    # Model fluxes (only for *observed* photometry)
    for iSED, imodel in enumerate(imodels):
        if choice == 'best':
            color = model_color
        else:
            color = SED_colors[iSED]

        #xlam = np.array(lam) * 1
        #xlam = xlam + 0.01 * (iSED - 0.5)
        plt.plot(lam, modelfnus[imodel], 's', mec=color, mfc='None', markersize=10, mew=2, alpha=0.8)

    #fnu_min = np.min(fnu.value)
    #fnu_max = np.max(fnu.value)
    
    model_flux_max = np.max(modelfnus.value.flat)
    #model_flux_max = np.percentile(modelfnus.value.flat, 95)
        
    flux_max = np.max([fnu_max, model_flux_max])
    if plot_full_spec:
        flux_max = np.max([flux_max, spec_max])
    
    flux_min = np.min([fnu_min, 0])

    flux_margin = 0.1 * (flux_max - flux_min)
    yhi = flux_max + flux_margin * 2
    ylo = flux_min - flux_margin
    
    #if choice == 'best':
    if choice in ('best', 'one'):
        if plot_params:
            if show_full_SED_flux_limits:
                ylo, yhi = ax1.set_ylim()
            #xtext = 0.2
            xtext = 0.2 * ((lam_max - lam_min) / 5.) + lam_min
            if lam_log:
                xtext = lam_min * 1.1
            ytext = yhi - 0.05 * (yhi - ylo)
            ytext = yhi - 0.02 * (yhi - ylo)

            formation_age = fit.posterior.samples[age_label][imodel] * u.Gyr  # time ago
            age_of_universe = pipes.utils.cosmo.age(redshift).to(u.Gyr)
            if formation_age > age_of_universe:
                formation_age = age_of_universe
            formation_time = age_of_universe - formation_age
            print('z=', redshift, '; age of the universe =', age_of_universe)

            if not no_text:
                text = ''
                #if not plot_title:
                #print('plot_title = ', plot_title)
                if plot_title == False:
                    text += 'id = %d\n' % id
                text += 'z = %.2f\n' % redshift

                more_text = more_text and not less_text
            
                if more_text:
                    #age_form = fit.posterior.samples[age_label][imodel] * u.Gyr
                    #if age_form > age_of_universe:
                    #    age_form = age_of_universe
                    text += 'time$_{form}$ = %d Myr\n' % (formation_time.to(u.Myr).value)
                    text += 'age$_{form}$ = %d Myr\n' % (formation_age.to(u.Myr).value)
                if less_text:
                    text += 'age = %d Myr\n' % (fit.posterior.samples[age_label][imodel] * 1000)
                else:
                    text += 'age$_{mass}$ = %d Myr\n' % (fit.posterior.samples['mass_weighted_age'][imodel] * 1000)

                #if 'mu' in cat.columns:
                #    obj = cat[cat['id'] == id]
                #    magnification = obj['mu'][0]
                if magnification != 1:
                    if magnification >= 10:
                        text += 'magnification $\mu$ = %d\n' % magnification
                    else:
                        text += 'magnification $\mu$ = %.1f\n' % magnification
                    text += 'mass($\mu$) = %.1e $M_\odot$\n' % (10**fit.posterior.samples['stellar_mass'][imodel] / magnification)
                    #text += 'SFR($\mu$)  = %d $M_\odot$/yr\n' % (fit.posterior.samples['sfr'][imodel] / magnification)
                    if not less_text:
                        text += 'SFR($\mu$)  = %.1f $M_\odot$/yr\n' % (fit.posterior.samples['sfr'][imodel] / magnification)
                else:
                    #magnification = 1        
                    text += 'mass = %.1e $M_\odot$\n' % (10**fit.posterior.samples['stellar_mass'][imodel])
                    if not less_text:
                        text += 'SFR  = %d $M_\odot$/yr\n' % (fit.posterior.samples['sfr'][imodel])
        
                if not less_text:
                    #text += 'sSFR = %d / Gyr\n' % (10 ** (fit.posterior.samples['ssfr'][imodel] + 9))
                    text += 'sSFR = %d / Myr\n' % (10 ** (fit.posterior.samples['ssfr'][imodel] + 12))
                if more_text:
                    if sfh_function:
                        if sfh_function+':metallicity' in fit_keys:
                            text += 'metallicity = %.2f $Z_\odot$\n' % fit.posterior.samples[sfh_function+':metallicity'][imodel]
                if 'dust:Av' in fit_keys:
                    text += 'dust $A_V$ = %.2f mag\n' % fit.posterior.samples['dust:Av'][imodel]
                if 'dust:eta' in fit_keys:
                    text += 'birth dust: %.1f$\\times$\n' % fit.posterior.samples['dust:eta'][imodel]
                else:
                    if 'eta' in fit.fit_instructions['dust'].keys():
                        eta = fit.fit_instructions['dust']['eta']
                        text += 'birth dust: %g$\\times$\n' % eta
                if 'nebular:logU' in fit_keys:
                    text += 'log(U) = %.2f\n' % fit.posterior.samples['nebular:logU'][imodel]
                if 'EW_OIII+Hbeta' in fit_keys:
                    text += 'EW [OIII]+H$\\beta$ = %d$\AA$\n' % fit.posterior.samples['EW_OIII+Hbeta'][imodel]
                text += '$\chi^2$ = %.2f\n' % fit.posterior.samples['chisq_phot'][imodel]

                plt.text(xtext, ytext, text, fontsize=fontsize, color=text_color, va='top')
    else:
        plt.legend(loc=2)  # 8
    
    plt.xlabel('Wavelength ($\mu$m)')

    if magnification != 1:
        ylabel = 'Observed Flux'
    else:
        ylabel = 'Flux'

    if flux_units == u.nJy:
        ylabel += ' (nJy)'
    elif flux_units == u.uJy:
        ylabel += ' ($\mu$Jy)'

    ax1.set_ylabel(ylabel)

    ax1.set_xlim(lam_min, lam_max)
    if lam_log:
        ax1.semilogx()
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        ax1.xaxis.set_minor_formatter(ticker.FormatStrFormatter("%g"))
        
    if show_full_SED_flux_limits:
        ylo, yhi = plt.ylim()
    if fmin != None:
        ylo = fmin

    if fmax != None:
        yhi = fmax
        
    plt.ylim(ylo, yhi)
    #plt.ylim(-3000, 6000)

    ax1.axhline(0, c='0.50', lw=1, ls=':')

    #ax1 = plt.gca()
    secax = add_magnitude_axis(ax1, flux_units)
    #secax = plt.gca().secondary_yaxis('right', functions=(nJy2AB, AB2nJy))
    #secax.set_ylabel('magnitude (AB)')
    
    topax = add_rest_wavelength_axis(ax1, redshift)

    if plot_title != False:
        if plot_title == None:
            #plot_title = '#' + id
            plot_title = str(id)
            if choice == 'best':
                plot_title += ' best SED fit'
            elif choice == 'age range':
                plot_title += ' SED fits'
                #plot_title += ' SED fits for age range'

        ax1.set_title(plot_title)
    
    #plt.legend(loc=2)  # 8
    
    if show_plot:
        fig.show()
    
    if save_plot:
        print('SAVING', plot_file)
        fig.savefig(plot_file, dpi=100)

    if close_plot:
        fig.close()
    else:
        return plt, fig, ax1, secax
        #return plt

params_file = 'BAGPIPES parameters.csv'
params_file = os.path.join(pipes.utils.install_dir, params_file)
params_labels = ascii.read(params_file)
params_labels['title label'] = astropy.table.column.MaskedColumn(params_labels['title label'], dtype='S256')
#params_labels['title label'].dtype = 'S20'
for ikey, key in enumerate(params_labels['key']):
    if 0: # ('age' in key.split('_')): # or ('sfr' in key.split('_')) or ('ssfr' in key.split('_')):
        params_labels['plot label'][ikey]  = 'log ' + params_labels['plot label'][ikey]
        params_labels['title label'][ikey] = 'log ' + params_labels['title label'][ikey]
    
catalog_keys_delayed = '''
redshift
magnification
stellar_mass
sfr
ssfr
sfh_function:age
formation_time
sfh_function:tau
mass_weighted_age
tform
sfh_function:metallicity
dust:Av
dust:eta
nebular:logU
chisq_phot
'''.split('\n')[1:-1]
#tquench


plot_ranges_dict = {}
plot_ranges_dict['stellar_mass'] = 6, 11
plot_ranges_dict['sfr'] = -4, 4
plot_ranges_dict['ssfr'] = 0, 2.5
plot_ranges_dict['dust:Av'] = 0, 2
plot_ranges_dict['redshift'] = 4, 11
plot_ranges_dict['tanefolds'] = -1, 1
plot_ranges_dict['formation_time'] = 0, 1000
plot_ranges_dict['age'] = 0, 2000
plot_ranges_dict['tau'] = 0, 2000
plot_ranges_dict['metallicity'] = -2.3, 0.7  # log metallicity Z / Zsun (0.005 - 5)
plot_ranges_dict['logU'] = -4, -1  # log(U) ionization parameter
plot_ranges_dict['EW_OIII+Hbeta'] = 0, 4  # log EW
#plot_ranges_dict['mass_weighted_age'] = 0, 3.3  # log age Myr
plot_ranges_dict['mass_weighted_age'] = 0, 1000

# corner plots:
#nsigma = np.arange(0.5, 2.1, 0.5)  # 0.5, 1.0, 1.5, 2.0-sigma
nsigma = np.arange(1,3.1)  # 1, 2, 3-sigma
levels = 1.0 - np.exp(-0.5 * nsigma ** 2)

def set_model_type(run_id):
    if 'bpass' in run_id:
        pipes.config.set_model_type('bpass')
    else:
        pipes.config.set_model_type('bc03_miles')

    #print(pipes.config.model_type)

def extract_age_sfh(fit):
    age_label = None
    sfh_function = None
    fit_keys = fit.posterior.samples.keys()
    for key in fit_keys:
        if ':age' in key:
            age_label = key  # delayed:age, constant:age_max
            sfh_function = key.split(':')[0]
    return age_label, sfh_function

def equivalent_width_total(model_lams, spec_flam, nlines=1):
    nspec = len(model_lams)
    EWs = []
    for ispec in range(1,nspec-1):
        Flam_line = spec_flam[ispec]
        flams_cont = spec_flam[ispec-1], spec_flam[ispec+1]
        flams_cont[1] / flams_cont[0]  # check they're similar
        Flam_cont = (flams_cont[1] + flams_cont[0]) / 2.
        dlam = (model_lams[ispec+1] - model_lams[ispec-1]) / 2.
        dlam = dlam.to(u.AA)
        Flam_flux = (Flam_line - Flam_cont) * dlam
        equivalent_width = EW = Flam_flux / Flam_cont
        #EW = EW / (1+z)  # lam already rest wavelength
        EWs.append(EW.value)

    EWs = np.array(EWs)
    
    EW_total = np.sum(np.sort(EWs)[-nlines:])
    return EW_total
    
#def save_results(run_id, id, cat, calc_EWs=[], sort='chisq_phot', overwrite=False):
def save_results(fit, run_id, id, magnification=1, calc_EWs=[], overwrite=False):
    #results_directory = os.path.join('results', run_id)
    results_directory = os.path.join('pipes/cats', run_id)
    os.makedirs(results_directory, exist_ok=True)
    data_file = os.path.join(results_directory, str(id)+'.cat')
    if os.path.exists(data_file):
        print(data_file, 'EXISTS')
        if not overwrite:
            return

    #fit = load_results(run_id, id, cat)  # fit or load
    fit_keys = fit.posterior.samples.keys()
    nsamples = len(fit.posterior.samples[list(fit_keys)[0]])
            
    fixed_redshift = 'redshift' not in fit_keys
    if fixed_redshift:
        redshift = fit.fit_instructions["redshift"]
        redshift = [redshift] * nsamples
        redshift = np.array(redshift)
    else:
        redshift = fit.posterior.samples['redshift']
            
    age_label, sfh_function = extract_age_sfh(fit)
    #print(age_label, sfh_function, 'asfh')
    if age_label:
        formation_age = fit.posterior.samples[age_label] * u.Gyr  # time ago
        #age_of_universe = pipes.utils.cosmo.age(fit.posterior.samples['redshift']).to(u.Gyr)
        age_of_universe = pipes.utils.cosmo.age(redshift).to(u.Gyr)
        formation_time = age_of_universe - formation_age
        #log_age_yr = np.log10(formation_age.to(u.yr).value)
        
    if not sfh_function:
        sfh_function = ''

    save_data = []
    catalog_labels = []
    formats = {}

    catalog_keys = catalog_keys_delayed + []
    if sfh_function == 'constant':
        i = catalog_keys.index('sfh_function:age')
        catalog_keys[i] = 'sfh_function:age_max'
        catalog_keys.remove('sfh_function:tau')
        catalog_keys.remove('dust:eta')

    #if 'mu' in cat.columns:
    #    obj = cat[cat['id'] == id]
    #    magnification = obj['mu'][0]
    #else:
    #    magnification = 1
    #    catalog_keys.remove('magnification')
    if magnification == 1:
        catalog_keys.remove('magnification')

    if calc_EWs:
        nmodels = len(fit.posterior.samples['chisq_phot'])
        model_spectra = fit.posterior.samples['spectrum_full']
        for spec_line in calc_EWs:
            print("Calculating Equivalent Widths for %s..." % spec_line)
            if spec_line == 'OIII+Hbeta':
                lam_min = 0.48 * u.um
                lam_max = 0.51 * u.um
                nlines = 3
            
            model_lams = fit.posterior.model_galaxy.wavelengths #* (1.+redshift)
            model_lams = model_lams * u.AA

            ilam_min = model_lams.searchsorted(lam_min)
            ilam_max = model_lams.searchsorted(lam_max)
            model_lams = model_lams[ilam_min:ilam_max]            

            EWs = []
            for imodel in range(nmodels):
                spec_flam = model_spectra[imodel] 
                spec_flam = spec_flam * u.erg / u.s / u.cm**2 / u.AA
                spec_flam = spec_flam[ilam_min:ilam_max]
                EW_total = equivalent_width_total(model_lams, spec_flam, nlines=nlines)
                EWs.append(EW_total)

            EWs = np.array(EWs)
            fit_key = catalog_key = 'EW_'+spec_line
            fit.posterior.samples[catalog_key] = EWs            
            catalog_keys.append(catalog_key)
            print('Done: max EW = %d A' % np.max(EWs))            
            
    for key in catalog_keys:
        if not sfh_function:
            if key == 'formation_time':
                continue
                
        #print(key, 'key')
        #key = key.replace('sfh_function:', '')  # delayed / constant
        key = key.replace('sfh_function', sfh_function)  # delayed / constant
        
        if key not in fit_keys:
            if ':tau' in key:
                key = key.replace(':tau', ':efolds')
                if key not in fit_keys:
                    key = key.replace(':efolds', ':tanefolds')
                    if key not in fit_keys:
                        print ('ERROR:', sfh_function, 'has no time scale parameter')

        if key not in fit_keys:
            if key == 'dust:eta':
                print(key, 'not in fit keys; skipping...')
                continue
                
        key_data = None
        print(key)
        if key in fit_keys:
            key_data = fit.posterior.samples[key] + 0
        elif key == 'magnification':
            key_data = [magnification] * nsamples
        elif key == 'formation_time':
            key_data = formation_time + 0 #.to(u.Myr).value + 0
        elif key == 'redshift':
            key_data = redshift + 0
        else:
            print ('ERROR:', key, 'has no data')
            continue
            #break

        key = key.replace(sfh_function+':', '')  # delayed / constant
        key = key.replace('nebular:', '')

        if key == 'ssfr':  # log /Gyr
            key_data += 9
        elif key == 'sfr':
            key_data /= magnification
            key_data = np.log10(key_data)
        elif key in 'stellar_mass massformed'.split():
            key_data -= np.log10(magnification)
        elif 0: # key in 'mass_weighted_age tform tquench'.split():
            key_data = key_data * u.Gyr
            key_data = key_data.to(u.Myr).value + 0

        #print(len(key_data), key_data)
        save_data.append(key_data)
        
        catalog_label = key
        for ikey, param_key in enumerate(params_labels['key']):
            catalog_key = key.replace('sfh_function', '')
            if param_key == catalog_key:
                catalog_label = params_labels['catalog label'][ikey]
                catalog_units = params_labels['catalog units'][ikey]
                if catalog_units:
                    catalog_label += '_(%s)' % catalog_units

        catalog_labels.append(catalog_label)
        formats[catalog_label] = '% 10.4f'

        #print(key, key_data[:5])

    save_data = np.array(save_data).T
    
    #print('data', save_data)
    
    # sort chisq
    
    #print(save_data.shape)
    #print(len(catalog_labels))
    #print(catalog_labels)    
    
    print('SAVING', data_file)
    ascii.write(save_data, data_file, names=catalog_labels, format='fixed_width', delimiter='')
    #ascii.write(save_data, data_file, names=catalog_labels, format='fixed_width', delimiter='', formats=formats)
    
    #return fit
    
def corner_plot_samples(run_id, cat_id, plot_keys):
    #results_directory = os.path.join('results', run_id)
    results_directory = os.path.join('pipes/cats', run_id)
    data_file = os.path.join(results_directory, str(cat_id)+'.cat')  # .dat
    cat = ascii.read(data_file)
    
    plot_labels = []
    title_labels = []
    catalog_labels = []
    for plot_key in plot_keys:
        plot_label = plot_key
        title_label = plot_key
        catalog_label = plot_key  # default
        for ikey, key in enumerate(params_labels['key']):
            if key == plot_key:
                plot_label = params_labels['plot label'][ikey]
                title_label = params_labels['title label'][ikey]
                catalog_label = params_labels['catalog label'][ikey]
                catalog_units = params_labels['catalog units'][ikey]
                if catalog_units:
                    catalog_label += '_(%s)' % catalog_units
            
        plot_labels.append(plot_label)
        title_labels.append(title_label)
        catalog_labels.append(catalog_label)
            
    #print(catalog_labels)
    #print(plot_labels)
    plot_samples = [cat[key] for key in catalog_labels]
    plot_samples = np.array(plot_samples).T
    for ikey, key in enumerate(plot_keys):
        #print(ikey, key)
        #if ('age' in key):  # or ('tau' in key) or ('form' in key):
        if ('age' in key) or ('tau' in key) or ('form' in key):
            plot_samples[:,ikey] *= 1000  # Gyr -> Myr
        if ('EW' in key) or ('metal' in key): # ('age' in key)
            plot_samples[:,ikey] = np.log10(plot_samples[:,ikey])
        
        if key == 'tanefolds':
            plot_samples[:,ikey] = -plot_samples[:,ikey]

    return plot_samples, plot_labels, title_labels, catalog_labels

import corner

#def corner_plot(run_id, id, cat, plot_keys, plot_ranges=None, color='darkorange', fig=None, run_title=None, save_plot=None, close_plot=False, overwrite=False):
def corner_plot(run_id, cat_id, plot_keys, plot_ranges=None, color='darkorange', fig=None, run_title=None, plot_title=True, save_plot=None, close_plot=False, overwrite=False):
    if save_plot:
        if os.path.exists(save_plot):
            if not overwrite:
                print(save_plot, 'EXISTS')
                return

    plot_samples, plot_labels, title_labels, catalog_labels = corner_plot_samples(run_id, cat_id, plot_keys)
            
    fig = corner.corner(plot_samples, labels=plot_labels, quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_kwargs={"fontsize": 12}, title_labels=title_labels,
                    fig=fig, # overplot on previous figure
                    #title_fmt='.2f',
                    smooth=1., smooth1d=1., bins=25, 
                    range=plot_ranges,
                    max_n_ticks=6,
                    #color='darkorange',
                    color=color,
                    #hist2d_kwargs={"levels":levels},
                    levels=levels)
    
    if plot_title != False:
        if plot_title == True:
            plot_title = cat_id
        fig.suptitle(plot_title, fontsize=30)

    if run_title:
        fig.text(0.95, 0.88, run_title, color='k', fontsize=25, transform=fig.transFigure, ha='right')

    if save_plot:
        print('SAVING', save_plot)
        plt.savefig(save_plot)

    if close_plot:
        plt.close()
    else:
        return fig

#plot_keys = 'redshift tanefolds mass_weighted_age'.split()
#corner_fig = corner_plot()

# Plot SFH input vs. model
# Time before observed

def fmtexp(x, pos):
    if 1e-4 < x < 1e4:
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

def between(lo, x, hi):
    return (lo < x) & (x < hi)
    
sfr_min_max_plotted = 1e30

def plot_sfh(fit, id, magnification=1, age_label='mass_weighted_age', sfh_choices='age_1-sigma', sfr_min=None, sfr_max=None, plot_tuniv=False, plot_tlog=True, t_max=13800, plot_title=None, save_plot=None, show_plot=True, close_plot=False, overwrite=False):
    global sfr_min_max_plotted, sfr_max_max_plotted

    if save_plot:
        if save_plot == True:
            # fit.fname 'pipes/posterior/delayed-tau_zphot/469_'
            plot_file = fit.fname.replace('posterior/', 'plots/')
            plot_file = os.path.dirname(plot_file)
            plot_file = os.path.join(plot_file, 'bagpipes_sfh')
            if plot_tuniv:
                plot_file += '-tuniv'
            else:
                plot_file += '-tago'
            plot_file += '_%d.png' % id
            #print('plot_file', plot_file)
        else:
            plot_file = save_plot
            
        if os.path.exists(plot_file) and not show_plot:
            if not overwrite:
                print(plot_file, 'EXISTS')
                return None

    sfh_samples = fit.posterior.samples["sfh"] / magnification

    #age_of_universe = fit.posterior.sfh.age_of_universe
    #now = age_of_universe / 1e6
    fit_keys = fit.posterior.samples.keys()
    fixed_redshift = 'redshift' not in fit_keys
    nsamples = len(fit.posterior.samples[list(fit_keys)[0]])
    if fixed_redshift:
        redshifts = fit.fit_instructions["redshift"]
        redshifts = [redshifts] * nsamples
        redshifts = np.array(redshifts)
    else:
        redshifts = fit.posterior.samples['redshift']
        
    zmin = np.min(redshifts)
    max_age_of_universe = pipes.utils.cosmo.age(zmin)
    max_age_of_universe = max_age_of_universe.to(u.Myr).value

    sfrmin, sfrlo, sfrmed, sfrhi, sfrmax = np.percentile(sfh_samples, (0, 16, 50, 84, 100), axis=0)
    time_ago = fit.posterior.sfh.ages  # time ago
    time_ago = time_ago / 1e6  # yr -> Myr

    #fig = plt.figure(figsize=(9, 4), dpi=100)
    fig = plt.figure(figsize=(8, 4), dpi=100)
    #fig = plt.figure(figsize=(8, 6), dpi=100)
    #ax = plt.subplot()

    if not plot_tuniv:
        plt.fill_between(time_ago, sfrmin, sfrmax, alpha=0.5, color='0.80')
    
    #plt.fill_between(t, sfrlo, sfrhi, alpha=1, color='navajowhite')
    #plt.plot(t, model.sfh.sfh, 'b', label='input')
    #plt.plot(t, sfrmed, 'orange', label='median')

    #sfrmax_plotted = 0
    sfr_min_max_plotted = 1e30
    #print(sfr_min_max_plotted, 'sfr_min_max_plotted0')
    
    def plot_sfh_sample(i, color='b', lw=1, label=None, zorder=1, plot_age=None):
        global sfr_min_max_plotted, sfr_max_max_plotted
        #print(sfr_min_max_plotted, 'sfr_min_max_plotted11')
        #print('sfh_sample', i, color)
        sfr = sfh_samples[i]
        sfr_max = np.max(sfr)
        #print(sfr_max, 'sfr_max1')

        if plot_tuniv:
            z = redshifts[i]
            age_of_universe = pipes.utils.cosmo.age(z).to(u.Myr).value
            t_univ = age_of_universe - time_ago
            t_plot = t_univ + 0
            plt.fill_between(t_plot, 0, sfr, color=color, alpha=0.15, lw=lw)
            plt.plot((t_plot[0], t_plot[0]), (0, sfr[0]), color=color, lw=lw, zorder=zorder) # vertical line at end (right hand side)
        else:
            t_plot = time_ago + 0
        
        plt.plot(t_plot, sfr, color=color, lw=lw, label=label, zorder=zorder)
        t_range = np.compress(sfr, t_plot)
        #t_lo = t_range[0]
        #t_hi = t_range[-1]
        #t_med = (t_lo + t_hi) / 2.
        #print(t_hi, t_lo, t_med, age_of_universe - t_med)
        #print(t_med)
        #print(age_of_universe - t_med)
        #print(age_of_universe - t_hi)
        if plot_age:
            j = np.searchsorted(time_ago, plot_age)
            plt.plot(t_plot[j], sfr[j], 'o', color=color)
        
        #print(sfr_min_max_plotted, 'sfr_min_max_plotted1')
        try:
            sfr_max_max_plotted = np.nanmax([sfr_max_max_plotted, sfr_max])
            sfr_min_max_plotted = np.nanmin([sfr_min_max_plotted, sfr_max])
        except:
            #print('EXCEPT')
            sfr_max_max_plotted = sfr_max + 0       
            sfr_min_max_plotted = sfr_max + 0
        #print(sfr_min_max_plotted, 'sfr_min_max_plotted')
        #print(sfr_max_max_plotted, 'sfr_max_max_plotted')
        
    ages = fit.posterior.samples[age_label]
    imodelsort = np.argsort(ages)
    percentiles = np.array([16, 50, 84])
    isorts = percentiles/100. * len(ages)
    isorts = isorts.astype(int)
    imodels = imodelsort[isorts]
    
    # Plot all SFH samples
    if plot_tuniv:
        for imodel in imodelsort:
            plot_sfh_sample(imodel, '0.90', zorder=-100)

    if sfh_choices == 'age_1-sigma':
        SED_colors = ['c', 'navajowhite', 'r']
        SED_colors = ['c', 'y', 'r']
        SED_colors = ['b', (0.6,0.6,0), 'r']
        #SED_colors = ['b', 'gold', 'r']
        SED_ages = fit.posterior.samples[age_label][imodels]
        SED_labels = []
        SED_labels.append('%3d Myr (1-sigma younger)'          % (1000*SED_ages[0]))
        SED_labels.append('%3d Myr (median age mass-weighted)' % (1000*SED_ages[1]))
        SED_labels.append('%3d Myr (1-sigma older)'            % (1000*SED_ages[2]))
        for i in range(3):
            #print()
            #print(sfr_min_max_plotted, 'sfr_min_max_plotted00')

            #print(i, imodels[i], SED_colors[i], SED_labels[i])
            plot_sfh_sample(imodels[i], color=SED_colors[i], label=SED_labels[i], plot_age=1000*SED_ages[i])

    elif sfh_choices == 'best':
        SED_color = 'm'  # 'gold'  # (0.6,0.6,0)
        imodel = np.argmin(fit.posterior.samples['chisq_phot'])
        SED_age = fit.posterior.samples[age_label][imodel]
        SED_label = '%3d Myr (best)'          % (1000*SED_age)
        plot_sfh_sample(imodel, color=SED_color, label=SED_label, plot_age=1000*SED_age)
        print('max', fit.posterior.samples['constant:age_max'][imodel])
        print('min', fit.posterior.samples['constant:age_min'][imodel])
        
    elif sfh_choices == 'full_range':
        plot_sfh_sample(np.argmin(fit.posterior.samples['ssfr']), 'g', label='highest sSFR')
        plot_sfh_sample(np.argmin(sfh_samples[:,0]), 'darkred', label='lowest current SFR')
        plot_sfh_sample(np.argmax(sfh_samples[:,0]), 'c', label='highest current SFR')

        ntime_samples = sfh_samples.shape[1]
        #print(ntime_samples)
        plot_sfh_sample(np.argmax(sfh_samples[:,ntime_samples//2]), 'm', label='highest SFR midway')
        plot_sfh_sample(np.argmax(sfh_samples[:,int(ntime_samples*0.6)]), 'k', label='highest early SFR')

        plot_sfh_sample(np.argmin(fit.posterior.samples[age_label]), 'b', label='youngest')
        plot_sfh_sample(np.argmax(fit.posterior.samples[age_label]), 'r', label='oldest')

        n = len(fit.posterior.samples[age_label])
        i = fit.posterior.samples[age_label].argsort()[n//2]
        plot_sfh_sample(i, 'orange', label='median age', lw=2)

        i = np.argmin(fit.posterior.samples['chisq_phot'])
        plot_sfh_sample(i, 'lime', label='best fit', lw=3)

    #plt.xlim(0, now)
    #if plot_tuniv:
    #    plt.xlim(0, max_age_of_universe)
    #else:
    #    plt.xlim(max_age_of_universe, 0)
    #plt.ylim(1e-3, 100)

    sfrmax_all = np.nanmax(sfh_samples.flat)
    #print(sfr_max_max_plotted, 'sfr_max_max')
    #print(1.2*sfrmax_all, '1.2*sfrmax_all')
    yhi = np.nanmax([1.2*sfrmax_all, 10 * sfr_max_max_plotted])
    ylo = sfr_min_max_plotted / 10.
    #print(sfr_min_max_plotted, 'sfr_min_max')
    #print(ylo, yhi)
    
    if sfr_min:
        ylo = sfr_min
    if sfr_max:
        yhi = sfr_max

    #print(ylo, yhi)

    plt.ylim(ylo, yhi)
    
    #plt.ylim(0.01, yhi)
    #plt.ylim(1e-3*yhi, 10*yhi)
    #plt.ylim(1e-2, 300)
    #plt.xlim(120, 1)
    if plot_tuniv:
        plt.xlim(0, t_max)
        plt.semilogy()
        plt.xlabel('Age of the Universe (Myr)')
    else:
        if plot_tlog:
            plt.xlim(t_max, 1)  # 2000
            plt.loglog()
        else:
            plt.xlim(t_max, 0)  # 2000
            plt.semilogy()
        plt.xlabel('Time before observed (Myr)')
        
    plt.ylabel('SFR ($M_\odot$ / yr)')
    #plt.ylim(1e-8, 300)
    #ylo, yhi = plt.ylim()
    
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(fmtexp))

    plt.legend(loc=2)
    if plot_title == None:
        #plot_title = '#' + id
        plot_title = str(id)
        plot_title += ' star formation history'
        
    plt.title(plot_title)

    if magnification != 1:
        if magnification >= 10:
            magnif_str = '%d' % magnification
        else:
            magnif_str = '%.1f' % magnification
        magnification_text = 'corrected for %s$\\times$ lensing magnification' % magnif_str
        fig.text(0.935, 0.835, magnification_text, color='k', fontsize=9, transform=fig.transFigure, ha='right')

    if show_plot:
        fig.show()

    #print('save_plot', save_plot, plot_file)
    if save_plot:
        if overwrite or not os.path.exists(plot_file):
            print('SAVING', plot_file)
            plt.savefig(plot_file)

    if close_plot:
        plt.close()
    else:
        return plt
        
def get_model_fluxes(model, magnification=1):    
    speclam, specflam = model.spectrum.T
    speclam = (speclam * u.AA).to(u.micron)
    specfnu = (specflam * u.erg / u.s / u.cm**2 / u.AA).to(u.nJy, u.spectral_density(speclam))
    specfnu = specfnu * magnification

    photlam  = model.filter_set.eff_wavs
    photflam = model.photometry
    photlam = (photlam * u.AA).to(u.micron)
    photfnu = (photflam * u.erg / u.s / u.cm**2 / u.AA).to(u.nJy, u.spectral_density(photlam))    
    photfnu = photfnu * magnification

    return speclam, specfnu, photlam, photfnu


def get_best_fit_model_fluxes(fit, flux_units=u.nJy, wave_units=u.micron):
    imodel = np.argmin(fit.posterior.samples['chisq_phot'])  # simply the best

    redshift = fit.posterior.samples['redshift'][imodel]    
    spec_lam = fit.posterior.model_galaxy.wavelengths * (1.+redshift)
    spec_lam = (spec_lam * u.AA).to(wave_units)

    spec_flam = fit.posterior.samples['spectrum_full'][imodel]
    spec_fnu = (spec_flam * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(spec_lam))

    phot_lam = fit.galaxy.filter_set.eff_wavs
    phot_lam = (phot_lam * u.AA).to(wave_units)

    phot_flam = fit.posterior.samples['photometry']
    phot_fnu = (phot_flam * u.erg / u.s / u.cm**2 / u.AA).to(flux_units, u.spectral_density(phot_lam))
    
    return spec_lam, spec_fnu, phot_lam, phot_fnu, redshift

def extract_best_model(fit):
    model_components = deepcopy(fit.fit_instructions)
    imodel = np.argmin(fit.posterior.samples['chisq_phot'])
    chisq = fit.posterior.samples['chisq_phot'][imodel]
    redshift = fit.posterior.samples['redshift'][imodel]

    input_keys = fit.fit_instructions.keys()
    output_keys = fit.posterior.samples.keys()
    for key1 in input_keys:
        value1 = fit.fit_instructions[key1]
        if type(value1) == dict:
            for key2 in value1.keys():
                value2 = value1[key2]
                key = key1+':'+key2
                if key in output_keys:
                    model_components[key1][key2] = fit.posterior.samples[key][imodel] + 0
        else:
            if key1 in output_keys:
                model_components[key1] = fit.posterior.samples[key1][imodel] + 0

    return model_components, chisq
#!/usr/bin/env python
# coding: utf-8

# To generate plots for all the galaxies,
# don't run the .ipynb notebook (it will get bogged down).a
# Instead, un this .py script from the terminal

# # Exploring EAZY results
# 
# * Load and print results
# * Galaxy card compiling color image and plots

# Earendel + Sunrise Arc lensed by WHL0137:  
# https://cosmic-spring.github.io/earendel.html  
# 
# Interactive color image + catalog explorer:  
# https://cosmic-spring-jwst.herokuapp.com  
# 
# WHL0137 images and galaxy catalogs processed by grizli and EAZY:  
# https://s3.amazonaws.com/grizli-v2/JwstMosaics/v4/index.html
# 
# Color images:  
# https://www.easyzoom.com/imageaccess/126a44c2acae47df81bee5c85c98d4c7  
# https://stsci.box.com/s/cq992a7cfd13wlsmu784zkioho1pga65  
# 
# EAZY:  
# https://eazy-py.readthedocs.io

# 

# ## Imports

# In[1]:


#pip install -U memory_profiler


# In[2]:


from memory_profiler import profile


# In[3]:


import eazy
import eazy.hdf5


# In[4]:


from astropy.table import Table
import astropy.units as u
import astropy

# for labeling color image
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.wcs as wcs


# In[5]:


# to show segmentation map
import photutils


# In[6]:


import numpy as np
from glob import glob
import string
from importlib import reload

import os
from os.path import expanduser
home = expanduser("~")


# In[7]:


# color images
import PIL
from PIL import Image, ImageDraw, ImageFont
PIL.Image.MAX_IMAGE_PIXELS = 933120000  # allow to load large images avoiding DecompressionBombError


# In[8]:


# plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
#%matplotlib inline   # non-interactive (easier for notebook scrolling)
#get_ipython().run_line_magic('matplotlib', 'notebook')
#plt.style.use(os.path.join(home, 'p', 'matplotlibrc.txt')) # https://matplotlib.org/tutorials/introductory/customizing.html
plt.style.use('https://www.stsci.edu/~dcoe/matplotlibrc.txt') # https://matplotlib.org/tutorials/introductory/customizing.html
mpl_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter, FuncFormatter, MultipleLocator


# ## Helper functions

# In[9]:


import coeSEDplot


# In[10]:


def find_object_in_catalog(RA, Dec):
    object_coordinates = SkyCoord(ra=RA*u.deg, dec=Dec*u.deg)
    
    # First look in segmentation map:
    x, y = image_wcs.world_to_pixel(object_coordinates)
    x = roundint(x)
    y = roundint(y)
    id = segm.data[y,x]
    
    # If nothing there, check for nearby object
    if not id:
        catalog_coordinates = SkyCoord(ra=eazy_results['ra'], dec=eazy_results['dec'])  # *u.deg
        idx, d2d, d3d = object_coordinates.match_to_catalog_sky(catalog_coordinates)
        id = eazy_results['id'][idx]
    
    return id


# In[11]:


def roundint(x):
    return np.round(x).astype(int)
    
def extract_id(cat, id, idlabel='id'): # choose_object select_object
    duck_duck_goose = cat[idlabel] == int(id)
    if len(duck_duck_goose):
        obj = cat[duck_duck_goose]
    return obj[0]


# # START HERE

# # Load EAZY results

# In[12]:


inroot = 'sunrise-grizli-v4.0-fix'
field = 'whl0137'


# In[13]:


eazy_file = inroot + '.eazypy.zout.fits'
eazy_results = Table.read(eazy_file)
eazy_results[:2]


# In[14]:


eazy_phot = Table.read(inroot + '_phot_apcorr.fits')
eazy_phot[:2]


# In[15]:


#eazy_data = fits.open(inroot + '.eazypy.data.fits')
#eazy_data.info()


# In the .h5 file, EAZY saved local paths to templates and FILTER.RES.latest.   
# Link to them on your machine so we can find them now! For example:
# 
# ln /Users/dcoe/miniconda3/envs/erophot/lib/python3.10/site-packages/eazy/data/templates   
# ln /Users/dcoe/miniconda3/envs/erophot/lib/python3.10/site-packages/eazy/data/filters/FILTER.RES.latest
# 
# New Carnall templates from Github sfhz:  
# https://github.com/gbrammer/eazy-photoz/tree/master/templates/sfhz  
# cd /Users/dcoe/miniconda3/envs/erophot/lib/python3.10/site-packages/eazy/data/templates  
# ln sfhz xspline_templates

# In[16]:


# Load detailed results from h5 file
h5file  = inroot + '.eazypy.h5'
eazy_run = eazy.hdf5.initialize_from_hdf5(h5file=h5file)


# In[17]:


eazy_run.cat[:2]


# # Combined output file: photometry + photo-z's

# In[18]:


outfile = field + '_phot-eazy.ecsv'
if os.path.exists(outfile):
    print('Loading', outfile)
    catalog = astropy.io.ascii.read(outfile)


# In[19]:


catalog[:2]


# # Load color image and segmentation map

# In[20]:


# Color image
field = 'whl0137'
color_image_file = '../color/%s_v4_bright.png' % field
im = Image.open(color_image_file)
color_image = np.asarray(im)
color_image = color_image[::-1]  # flip top-bottom
color_image_file


# In[21]:


color_image.shape


# In[22]:


# Segmentation map
#segm_file = os.path.join('../phot', field+'_total_detections_segm.fits.gz')
#segm_file = os.path.join('../catalogs', 'sunrise-grizli-v2-ir_seg.fits')
#segm_file = inroot + '-ir_seg.fits.gz'
segm_file = 'sunrise-grizli-v4.0-ir_20mas_seg.fits'
segm_data = fits.open(segm_file)[0].data
segm = photutils.segmentation.SegmentationImage(segm_data)


# In[23]:


segm_data.shape


# In[24]:


#image_files = glob('../images/*.fits')
#infile = image_files[0]
hdu = fits.open(segm_file)
idata = 0
image_wcs = wcs.WCS(hdu[idata].header, hdu)
#image_wcs


# In[25]:


# Calculate x,y -- to be used for image stamps
#catalog_coordinates = SkyCoord(ra=eazy_results['ra']*u.deg, dec=eazy_results['dec']*u.deg)
catalog_coordinates = SkyCoord(ra=eazy_results['ra'], dec=eazy_results['dec'])  # *u.deg
eazy_results['x'], eazy_results['y'] = image_wcs.world_to_pixel(catalog_coordinates)
eazy_run.cat['x'], eazy_run.cat['y'] = image_wcs.world_to_pixel(catalog_coordinates)


# In[26]:


eazy_results[:2]


# In[27]:


coeSEDplot.eazy_results = eazy_results
coeSEDplot.eazy_run = eazy_run
coeSEDplot.segm = segm
coeSEDplot.color_image = color_image


# # Look at a galaxy

# Say you're browsing the interactive image + catalog https://cosmic-spring-jwst.herokuapp.com  
# and find an interesting galaxy you want to look at. It will give you the ID number that you can use here.  
# Alternatively, you might have the RA, Dec coordinates, which you can use further below.

# # Show galaxy

# In[28]:


#get_ipython().run_line_magic('matplotlib', 'inline')


# In[29]:


reload(coeSEDplot)


# In[30]:


id = 8
obj = extract_id(eazy_results, id)
# z_phot -1 is junk
obj


# In[31]:


obj['ssfr_p']


# In[32]:


reload(coeSEDplot)
#fig = coeSEDplot.show_galaxy_properties(id, figsize=(12, 6), lam_min=0.4, lam_max=5, save=True); 
#, alt_id='1.31');


# In[33]:


#id = 180
#fig = coeSEDplot.show_galaxy_properties(id, figsize=(12, 6), lam_min=0.4, lam_max=5, save=True);
#plt.close(fig)


# # Plot and save all SED fits

# In[ ]:


for id in catalog['id']:
    #plt.close()
    #plt.close('all')
    #if id in [261,1363]:
    #    continue
    try:
        obj = extract_id(eazy_results, id)
        if obj['z_phot'] < 0:
            print('Object #%d seems like junk (z_phot < 0)' % id)
        else:
            outfile = 'plots/eazy_%d.png' % id
            if os.path.exists(outfile):
                print(outfile, 'EXISTS')
            else:
                print('Object #%d plotting now...' % id)
                fig = coeSEDplot.show_galaxy_properties(id, figsize=(12, 6), lam_min=0.4, lam_max=5, save=True);
                plt.close(fig)
                #, alt_id='1.31');
    except:
        print('FAILED')
        pass

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import subprocess
import numpy as np
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import pickle as pkl
from astropy.io import fits
import astropy.units as u

data_dir = '/home/stella/1st-year-project/data'
default_config = '/home/stella/1st-year-project/sex/config.sex'
dict_directory = '/home/stella/1st-year-project/dict'
dualsex_dir = '/home/stella/1st-year-project/dualsex_NIRCam'

#dictionaries of filter properties
zp_dict = {'f277w':26.47802928633102-0.061928578858356786, 'f356w': 26.47802928633102+0.08875112502714799 ,'f444w':26.47802928633102+0.12791695886634002, 'f770w':25.265227862373404, 'f1000w':25.265227862373404, 'f1500w':25.265227862373404, 'f1800w':25.265227862373404}
pixelsize_dict = {'f444w':0.062912, 'f770w':0.110915, 'f1000w':0.110915, 'f1500w':0.110915, 'f1800w':0.110915}
FWHM_dict = {'f444w':0.14, 'f770w':0.25, 'f1000w':0.32,'f1500w':0.48, 'f1800w':0.58}

#detection limits
detectMult = 1
MINAREA_dict = {'f444w':4, 'f770w':4, 'f1000w':6,'f1500w':15, 'f1800w':15}
THRESH_dict = {'f444w':1*detectMult, 'f770w':1*detectMult, 'f1000w':0.67*detectMult,'f1500w':0.26*detectMult, 'f1800w':0.26*detectMult}

##########
# CONFIG #
##########
filter_list = ['f444w','f356w','f277w']
detection_filter_id = 0
measurement_filter_ids = [1,2]

#############
# ADDRESSES #
#############
filter_dir = {}
for observation in glob.glob('%s/*' %data_dir):
    if 'niriss' not in observation:
        filter_dir[observation.split('_')[-1].split('clear-')[-1]] = observation

sci_array = []
wht_array = []
for filter_name in filter_list:
    contents = glob.glob('%s/*' %filter_dir[filter_name])
    for content in contents:
        print(content)
        # if 'dual' in content:
        if 'Acorrected' in content:
            sci_array.append(content)
            wht_array.append('%s/wht.fits' %filter_dir[filter_name])

###########
# FRAMING #
###########
filter_id = 0
fits_file = fits.open(sci_array[filter_id])
shape_y, shape_x = np.shape(fits_file[0].data)
print('old shape of ', filter_list[filter_id], shape_y, shape_x)
fits_file[0].data = fits_file[0].data[1:-1,1:-1]
shape_y, shape_x = np.shape(fits_file[0].data)
print('new shape of ', filter_list[filter_id], shape_y, shape_x)
fits_file.writeto('%s/nircam_sci.fits' %filter_dir[filter_list[filter_id]], overwrite=True)

wht_file = fits.open(wht_array[filter_id])
wht_file[0].data = wht_file[0].data[1:-1,1:-1]
wht_file.writeto('%s/nircam_wht.fits' %filter_dir[filter_list[filter_id]], overwrite=True)

#
filter_id = 1
fits_file = fits.open(sci_array[filter_id])
shape_y, shape_x = np.shape(fits_file[0].data)
print('old shape of ', filter_list[filter_id], shape_y, shape_x)
fits_file[0].data = fits_file[0].data[0:-2,2:-1]
shape_y, shape_x = np.shape(fits_file[0].data)
print('new shape of ', filter_list[filter_id], shape_y, shape_x)
fits_file.writeto('%s/nircam_sci.fits' %filter_dir[filter_list[filter_id]], overwrite=True)

wht_file = fits.open(wht_array[filter_id])
wht_file[0].data = wht_file[0].data[0:-2,2:-1]
wht_file.writeto('%s/nircam_wht.fits' %filter_dir[filter_list[filter_id]], overwrite=True)

#
filter_id = 2
fits_file = fits.open(sci_array[filter_id])
shape_y, shape_x = np.shape(fits_file[0].data)
print('old shape of ', filter_list[filter_id], shape_y, shape_x)
fits_file[0].data = fits_file[0].data[1:-2,1:-2]
shape_y, shape_x = np.shape(fits_file[0].data)
print('new shape of ', filter_list[filter_id], shape_y, shape_x)
fits_file.writeto('%s/nircam_sci.fits' %filter_dir[filter_list[filter_id]], overwrite=True)

wht_file = fits.open(wht_array[filter_id])
wht_file[0].data = wht_file[0].data[1:-2,1:-2]
wht_file.writeto('%s/nircam_wht.fits' %filter_dir[filter_list[filter_id]], overwrite=True)

#############
# ADDRESSES #
#############
filter_dir = {}
for observation in glob.glob('%s/*' %data_dir):
    if 'niriss' not in observation:
        filter_dir[observation.split('_')[-1].split('clear-')[-1]] = observation

sci_array = []
wht_array = []
for filter_name in filter_list:
    contents = glob.glob('%s/*' %filter_dir[filter_name])
    for content in contents:
        # if 'dual' in content:
        if 'nircam_sci' in content:
            sci_array.append(content)
            wht_array.append('%s/nircam_wht.fits' %filter_dir[filter_name])

#############
# DETECTION #
#############
#this is for the detection filter only
f = open(default_config,'r')
lines = f.readlines()
f.close()

lines[6]  = 'CATALOG_NAME     %s/dual_%s.cat           # name of the output catalog\n' %(dualsex_dir, filter_list[detection_filter_id])
lines[14] = 'DETECT_MINAREA   %s                       # min. # of pixels above threshold\n' %MINAREA_dict[filter_list[detection_filter_id]]
lines[15] = 'DETECT_THRESH    %s                       # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n' %THRESH_dict[filter_list[detection_filter_id]]
lines[16] = 'ANALYSIS_THRESH  %s                       # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n' %THRESH_dict[filter_list[detection_filter_id]]
lines[33] = 'WEIGHT_IMAGE     %s/wht.fits              #<detection>,<measurement> weight-image(s)\n' %filter_dir[filter_list[detection_filter_id]]
lines[34] = 'WEIGHT_TYPE      MAP_WEIGHT               # weighting scheme (for single image, or\n'
lines[49] = 'PHOT_APERTURES   %s                       # MAG_APER aperture diameter(s) in pixels\n' %(1.0/pixelsize_dict[filter_list[detection_filter_id]])
lines[59] = 'MAG_ZEROPOINT    %s                       # magnitude zero-point\n' %zp_dict[filter_list[detection_filter_id]]
lines[63] = 'PIXEL_SCALE      %s                       # size of pixel in arcsec (0=use FITS WCS info)\n' %pixelsize_dict[filter_list[detection_filter_id]]
lines[67] = 'SEEING_FWHM      %s                       # stellar FWHM in arcsec\n' %FWHM_dict[filter_list[detection_filter_id]]
lines[84] = 'CHECKIMAGE_NAME  %s/dual_%s.fits          # Filename for the check-image\n' %(dualsex_dir, filter_list[detection_filter_id])

sex_address = '%s/dual_%s.sex' %(dualsex_dir, filter_list[detection_filter_id])
f = open(sex_address,'w')
f.writelines(lines)
f.close()

sex_command = 'sex %s -c %s -WEIGHT_IMAGE %s' %(sci_array[detection_filter_id], sex_address, wht_array[detection_filter_id])
p = subprocess.Popen(sex_command, stdout=subprocess.PIPE, shell=True)
print(p.communicate())

##############
# MEASURMENT #
##############
for measurement_filter_id in measurement_filter_ids:
    f = open(default_config,'r')
    lines = f.readlines()
    f.close()

    lines[6]  = 'CATALOG_NAME     %s/dual_%s.cat           # name of the output catalog\n' %(dualsex_dir, filter_list[measurement_filter_id])
    lines[14] = 'DETECT_MINAREA   %s                       # min. # of pixels above threshold\n' %MINAREA_dict[filter_list[detection_filter_id]]
    lines[15] = 'DETECT_THRESH    %s                       # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n' %THRESH_dict[filter_list[detection_filter_id]]
    lines[16] = 'ANALYSIS_THRESH  %s                       # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n' %THRESH_dict[filter_list[detection_filter_id]]
    lines[33] = 'WEIGHT_IMAGE     %s/wht.fits, %s/wht.fits #<detection>,<measurement> weight-image(s)\n' %(filter_dir[filter_list[detection_filter_id]], filter_dir[filter_list[measurement_filter_id]])
    lines[34] = 'WEIGHT_TYPE      MAP_WEIGHT               # weighting scheme (for single image, or\n'
    lines[49] = 'PHOT_APERTURES   %s                       # MAG_APER aperture diameter(s) in pixels\n' %(1.0/pixelsize_dict[filter_list[detection_filter_id]])
    lines[59] = 'MAG_ZEROPOINT    %s                       # magnitude zero-point\n' %(zp_dict[filter_list[measurement_filter_id]])
    lines[63] = 'PIXEL_SCALE      %s                       # size of pixel in arcsec (0=use FITS WCS info)\n' %pixelsize_dict[filter_list[detection_filter_id]]
    lines[67] = 'SEEING_FWHM      %s                       # stellar FWHM in arcsec\n' %FWHM_dict[filter_list[detection_filter_id]]
    lines[84] = 'CHECKIMAGE_NAME  %s/dual_%s.fits          # Filename for the check-image\n' %(dualsex_dir, filter_list[measurement_filter_id])

    sex_address = '%s/dual_%s.sex' %(dualsex_dir, filter_list[measurement_filter_id])
    f = open(sex_address,'w')
    f.writelines(lines)
    f.close()

    sex_command = 'sex %s,%s -c %s -WEIGHT_IMAGE %s,%s' %(sci_array[detection_filter_id], sci_array[measurement_filter_id], sex_address, wht_array[detection_filter_id], wht_array[measurement_filter_id])
    p = subprocess.Popen(sex_command, stdout=subprocess.PIPE, shell=True)
    print(p.communicate())

##############
# CATALOGING #
##############
dual_catalogue = {}
for filter_name in filter_list:
    dual_catalogue[filter_name] = {}

    cat_address = '%s/dual_%s.cat' %(dualsex_dir, filter_name)
    f = open(cat_address,'r')
    lines = f.readlines()
    f.close()

    header_list = []
    count = 0
    for line in lines[:40]:
        if line.split()[0] == '#':
            header_list.append(line.split()[2])
            dual_catalogue[filter_name][line.split()[2]] = []
            first_cat_line = count + 1
        count = count + 1

    for line in lines[first_cat_line:]:
        values = line.split()
        for i in range(len(header_list)):
            if i == 0:
                dual_catalogue[filter_name][header_list[i]].append(int(values[i]))
            else:
                dual_catalogue[filter_name][header_list[i]].append(values[i])

with open('%s/dual_catalogue_NIRCAMONLY.pickle' %dict_directory, 'wb') as handle:
    pkl.dump(dual_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

############
# MATCHING #
############
objects_array = dual_catalogue['f444w']['NUMBER']
matched_catalogue = {}
for object in objects_array:

    matched_catalogue[object] = {}
    i = object - 1

    for filter_name in filter_list:
        matched_catalogue[object][filter_name] = {}

        for header in dual_catalogue[filter_name].keys():
            matched_catalogue[object][filter_name][header] = dual_catalogue[filter_name][header][i]

with open('%s/dual_matched_catalogue_NIRCAMONLY.pickle' %dict_directory, 'wb') as handle:
    pkl.dump(matched_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

#############################
# Checking the AGN criteria #
#############################

# matching these with the MIRI catalog
with open('%s/table_matched_catalogue_f444w.pickle' %dict_directory, 'rb') as handle: # change this to the right path
    MIRI_matched_catalogue = pkl.load(handle)

all_catalogue = {}

ra_threshold_arcsec = 0.6 #arcsec
dec_threshold_arcsec = 0.6 #arcsec
threshold_deg = SkyCoord(ra_threshold_arcsec*u.arcsec, dec_threshold_arcsec*u.arcsec, frame='icrs')

for MIRI_object in tqdm(list(MIRI_matched_catalogue.keys())):
    ra_MIRI = float(MIRI_matched_catalogue[MIRI_object]['f444w']['ALPHA_SKY'])
    dec_MIRI = float(MIRI_matched_catalogue[MIRI_object]['f444w']['DELTA_SKY'])

    for NIRCam_object in matched_catalogue.keys():
        ra_f444w = float(matched_catalogue[NIRCam_object]['f444w']['ALPHA_SKY'])
        dec_f444w = float(matched_catalogue[NIRCam_object]['f444w']['DELTA_SKY'])
        class_star_f444w = float(matched_catalogue[NIRCam_object]['f444w']['CLASS_STAR'])

        if class_star_f444w < 1.5: #not excluding based on this yet
            ra_diff = abs(ra_MIRI-ra_f444w)
            dec_diff = abs(dec_MIRI-dec_f444w)

            if ra_diff < threshold_deg.ra.value and dec_diff < threshold_deg.dec.value:
                #print('matched')
                if MIRI_object in all_catalogue.keys(): print('twice')
                else:
                    all_catalogue[MIRI_object] = MIRI_matched_catalogue[MIRI_object]
                    all_catalogue[MIRI_object]['f277w'] = matched_catalogue[NIRCam_object]['f277w']
                    all_catalogue[MIRI_object]['f356w'] = matched_catalogue[NIRCam_object]['f356w']
                    all_catalogue[MIRI_object]['f444w_2'] = matched_catalogue[NIRCam_object]['f444w']

with open('%s/all_catalogue.pickle' %dict_directory, 'wb') as handle:
    pkl.dump(all_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

for filter in filter_list:
    print(all_catalogue[2][filter].keys)
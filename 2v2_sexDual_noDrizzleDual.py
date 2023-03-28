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
from reproject import reproject_adaptive
import os
import astropy.units as u

main_dir = '/home/stella/1st-year-project'
data_dir = '%s/JWST' %main_dir
processedData_dir = '%s/processedData' %main_dir
dict_dir = '%s/dict' %main_dir
sex_dir = '%s/sex' %main_dir
dualsex_dir = '%s/dualsex' %main_dir
default_config = '%s/config.sex' %sex_dir

if not os.path.exists(dict_dir): os. makedirs(dict_dir)
if not os.path.exists(processedData_dir): os. makedirs(processedData_dir)
if not os.path.exists(sex_dir): os. makedirs(sex_dir)
if not os.path.exists(dualsex_dir): os. makedirs(dualsex_dir)

#dictionaries of miri filter properties
#zp_dict        = {'f560w':25.265227862373404, 'f770w':25.265227862373404,  'f1000w':25.265227862373404}
#pixelsize_dict = {'f560w':0.11,               'f770w':0.11,                'f1000w':0.11}
#FWHM_dict      = {'f560w':0.22,               'f770w':0.25,                'f1000w':0.32}

#zp_dict        = {'f090w':26.64945925936983, 'f444w':27.54550869416636,  'f1000w':25.265227862373404}
#pixelsize_dict = {'f090w':0.0009453,               'f444w':0.00397,                'f1000w':0.11}
zp_dict = {}
pixelsize_dict = {}
FWHM_dict      = {'f090w':0.030,               'f444w':0.140, "f277w":0.088, "f356w":0.114,                'f1000w':0.32} #found by searching psf fwhm nircam

#detection limits
detectMult = 1
MINAREA_dict = {'f090w':3,            'f444w':3,            'f1000w':3}
THRESH_dict  = {'f090w':1*detectMult, 'f444w':1*detectMult, 'f1000w':1*detectMult}

##########
# CONFIG #
##########
filter_list = ['f277w', "f356w", 'f444w']  #'f150w', 'f200w', 'f277w', 'f356w', 
detection_filter_id = 2
measurement_filter_ids = [0, 1]

#####################
# SCI and WHT files #
#####################
filter_address = {}
for observation in glob.glob('%s/*/*' %data_dir):
    if 'i2d' in observation:
        filter = observation.split('_')[-2].split('-')[-1].lower()
        filter_address[filter] = observation
        
        #Find parameters from fits header
        header = fits.open(observation)[1].header
        
        f = header["PIXAR_SR"] * header["PHOTMJSR"] * 10**6
        mag = -2.5 * np.log10(f/3631)

        pixelsize = np.sqrt(header["PIXAR_A2"])

        zp_dict[filter] = mag
        pixelsize_dict[filter] = pixelsize
        
        if not filter in MINAREA_dict:
            MINAREA_dict[filter] = 3
        if not filter in THRESH_dict:
            THRESH_dict[filter] = 1*detectMult
    
for filter in filter_list:         
    assert (filter in FWHM_dict), "missing FWHM values"

sci_array = []
wht_array = []
for filter_name in filter_list:
    fits_file = fits.open(filter_address[filter_name])
    fits.writeto('%s/Original_SCI_%s.fits' %(processedData_dir, filter_name), fits_file['SCI'].data, fits_file['SCI'].header, overwrite=True)
    fits.writeto('%s/Original_WHT_%s.fits' %(processedData_dir, filter_name), fits_file['WHT'].data, fits_file['WHT'].header, overwrite=True)

    sci_array.append('%s/Original_SCI_%s.fits' %(processedData_dir, filter_name))
    wht_array.append('%s/Original_WHT_%s.fits' %(processedData_dir, filter_name))

#############
# DETECTION #
#############
#going through the filters and detect 3 point sources to find the pixel offset between the images
for filter_id in range(len(filter_list)):

    f = open(default_config,'r')
    lines = f.readlines()
    f.close()

    lines[6]  = 'CATALOG_NAME     %s/sex_%s.cat           # name of the output catalog\n' %(sex_dir, filter_list[filter_id])
    lines[14] = 'DETECT_MINAREA   %s                       # min. # of pixels above threshold\n' %MINAREA_dict[filter_list[filter_id]]
    lines[15] = 'DETECT_THRESH    %s                       # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n' %THRESH_dict[filter_list[filter_id]]
    lines[16] = 'ANALYSIS_THRESH  %s                       # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n' %THRESH_dict[filter_list[filter_id]]
    lines[33] = 'WEIGHT_IMAGE     %s                       # <detection>,<measurement> weight-image(s)\n' %wht_array[filter_id]
    lines[34] = 'WEIGHT_TYPE      MAP_WEIGHT               # weighting scheme (for single image, or\n'
    lines[49] = 'PHOT_APERTURES   %s                       # MAG_APER aperture diameter(s) in pixels\n' %(1.0/pixelsize_dict[filter_list[filter_id]])
    lines[59] = 'MAG_ZEROPOINT    %s                       # magnitude zero-point\n' %zp_dict[filter_list[filter_id]]
    lines[63] = 'PIXEL_SCALE      %s                       # size of pixel in arcsec (0=use FITS WCS info)\n' %pixelsize_dict[filter_list[filter_id]]
    lines[67] = 'SEEING_FWHM      %s                       # stellar FWHM in arcsec\n' %FWHM_dict[filter_list[filter_id]]
    lines[84] = 'CHECKIMAGE_NAME  %s/sex_%s.fits          # Filename for the check-image\n' %(sex_dir, filter_list[filter_id])

    sex_address = '%s/%s_new.sex' %(sex_dir, filter_list[filter_id])
    f = open(sex_address,'w')
    f.writelines(lines)
    f.close()

    #sex_command = 'sex'
    #p = subprocess.Popen(sex_command, stdout=subprocess.PIPE, shell=True)
    #print(p.communicate())

    sex_command = 'sex %s -c %s -WEIGHT_IMAGE %s' %(sci_array[filter_id], sex_address, wht_array[filter_id])
    p = subprocess.Popen(sex_command, stdout=subprocess.PIPE, shell=True)
    print(p.communicate())

##################
# cross-matching #
##################
#making the non-cross-matched catalogue
sex_catalogue = {}
for filter_name in filter_list:
    sex_catalogue[filter_name] = {}

    cat_address = '%s/sex_%s.cat' %(sex_dir, filter_name)
    f = open(cat_address,'r')
    lines = f.readlines()
    f.close()

    header_list = []
    count = 0
    for line in lines[:40]:
        if line.split()[0] == '#':
            header_list.append(line.split()[2])
            sex_catalogue[filter_name][line.split()[2]] = []
            first_cat_line = count + 1
        count = count + 1

    for line in lines[first_cat_line:]:
        values = line.split()
        for i in range(len(header_list)):
            if i == 0:
                sex_catalogue[filter_name][header_list[i]].append(int(values[i]))
            else:
                sex_catalogue[filter_name][header_list[i]].append(values[i])

with open('%s/sex_catalogue.pickle' %dict_dir, 'wb') as handle:
    pkl.dump(sex_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

#cross-matching the stars and finding the offset between the filters
# sexed_ids = [3, 4] #226 #these are the stars
reference_filter = 'f444w'
from copy import deepcopy
nonref_filteres = deepcopy(filter_list)
nonref_filteres.remove(reference_filter)
print(nonref_filteres)

sexed_ids = sex_catalogue[reference_filter]['NUMBER'] #226

ra_threshold_arcsec = 0.3 #arcsec
dec_threshold_arcsec = 0.3 #arcsec
threshold_deg = SkyCoord(ra_threshold_arcsec*u.arcsec, dec_threshold_arcsec*u.arcsec, frame='icrs')

matched_catalogue = {}
for sexed in tqdm(sexed_ids):

    matched_catalogue[sexed] = {}
    matched_catalogue[sexed][reference_filter] = {}
    for header in sex_catalogue[reference_filter].keys():
        matched_catalogue[sexed][reference_filter][header] = sex_catalogue[reference_filter][header][sexed-1]

    id_reference  = sex_catalogue[reference_filter]['NUMBER'][sexed-1]
    ra_reference  = float(sex_catalogue[reference_filter]['ALPHA_SKY'][sexed-1])
    dec_reference = float(sex_catalogue[reference_filter]['DELTA_SKY'][sexed-1])

    for filter_name in nonref_filteres:
        for object in range(len(sex_catalogue[filter_name]['NUMBER'])):

            class_star = float(sex_catalogue[filter_name]['CLASS_STAR'][object-1])
            ra_nonref  = float(sex_catalogue[filter_name]['ALPHA_SKY'][object-1])
            dec_nonref = float(sex_catalogue[filter_name]['DELTA_SKY'][object-1])

            ra_diff = abs(ra_reference-ra_nonref)
            dec_diff = abs(dec_reference-dec_nonref)

            if ra_diff < threshold_deg.ra.value and dec_diff < threshold_deg.dec.value:
                if class_star > -0.1:

                    if filter_name in matched_catalogue[sexed].keys(): print('twice!!!')
                    matched_catalogue[sexed][filter_name] = {}
                    for header in sex_catalogue[filter_name].keys():
                        matched_catalogue[sexed][filter_name][header] = sex_catalogue[filter_name][header][object-1]

#deleting everything that hasn't been matched in all the filters
delete_ids = []
for object in matched_catalogue:

    if len(matched_catalogue[object]) < 3: delete_ids.append(object)

for id in delete_ids:
    del(matched_catalogue[id])

#finding a linear pixel transformation that puts the stars on top of eachother
def offset_finder(matched_catalogue, filter_1, filter_2):

    x_offset_array = []
    y_offset_array = []
    for star in matched_catalogue.keys():
        print(star)
        x1 = float(matched_catalogue[star][filter_1]['X_IMAGE'])
        x2 = float(matched_catalogue[star][filter_2]['X_IMAGE'])
        x_offset = x2-x1
        x_offset_array.append(x_offset)

        y1 = float(matched_catalogue[star][filter_1]['Y_IMAGE'])
        y2 = float(matched_catalogue[star][filter_2]['Y_IMAGE'])
        y_offset = y2-y1
        y_offset_array.append(y_offset)

        x_offset = np.mean(x_offset_array)
        y_offset = np.mean(y_offset_array)

    return x_offset, y_offset

print(offset_finder(matched_catalogue, 'f444w', 'f277w'))
print(offset_finder(matched_catalogue, 'f444w', 'f356w'))

#generating new SCI and WHT files with adjusted pixel positions
targetShape_y, targetShape_x = np.shape(fits.open(sci_array[0])[0].data)
for filter_id in measurement_filter_ids:

    sci_fits = fits.open(sci_array[filter_id])
    wht_fits = fits.open(wht_array[filter_id])

    x,y = offset_finder(matched_catalogue, 'f444w', filter_list[filter_id])
    image = sci_fits[0].data[int(round(y)):, int(round(x)):]
    whtMap = wht_fits[0].data[int(round(y)):, int(round(x)):]

    #making the header
    sci_fits[0].header['CRPIX1'] = sci_fits[0].header['CRPIX1'] - x
    sci_fits[0].header['CRPIX2'] = sci_fits[0].header['CRPIX2'] - y

    #adding zero out of range
    target_image = np.zeros((targetShape_y, targetShape_x))
    target_whtMap = np.zeros((targetShape_y, targetShape_x))
    shape_y, shape_x = np.shape(image)
    for y in range(shape_y):
        for x in range(shape_x):
            if y < targetShape_y and x < targetShape_x:
                target_image[y,x] = image[y,x]
                target_whtMap[y,x] = whtMap[y,x]

    sci_fits[0].data = target_image
    wht_fits[0].data = target_whtMap

    fits.writeto('%s/reframed_SCI_%s.fits' %(processedData_dir, filter_list[filter_id]), sci_fits[0].data, sci_fits[0].header, overwrite=True)
    fits.writeto('%s/reframed_WHT_%s.fits' %(processedData_dir, filter_list[filter_id]), wht_fits[0].data, wht_fits[0].header, overwrite=True)

cp_command = 'cp %s/Original_SCI_f444w.fits %s/reframed_SCI_f444w.fits' %(processedData_dir, processedData_dir)
p = subprocess.Popen(cp_command, stdout=subprocess.PIPE, shell=True)
print(p.communicate())

cp_command = 'cp %s/Original_WHT_f444w.fits %s/reframed_WHT_f444w.fits' %(processedData_dir, processedData_dir)
p = subprocess.Popen(cp_command, stdout=subprocess.PIPE, shell=True)
print(p.communicate())

sci_array = []
wht_array = []
for filter_name in filter_list:
    sci_array.append('%s/reframed_SCI_%s.fits' %(processedData_dir, filter_name))
    wht_array.append('%s/reframed_WHT_%s.fits' %(processedData_dir, filter_name))

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
lines[33] = 'WEIGHT_IMAGE     %s                       # <detection>,<measurement> weight-image(s)\n' %wht_array[detection_filter_id]
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
    lines[33] = 'WEIGHT_IMAGE     %s, %s                   # <detection>,<measurement> weight-image(s)\n' %(wht_array[detection_filter_id], wht_array[measurement_filter_id])
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

with open('%s/dual_catalogue.pickle' %dict_dir, 'wb') as handle:
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

with open('%s/dual_matched_catalogue.pickle' %dict_dir, 'wb') as handle:
    pkl.dump(matched_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

############
# CLEANING #
############
clean_catalogue = {}
for object in objects_array:

    if float(matched_catalogue[object][filter_list[0]]['MAGERR_APER']) < 0.2:
        if float(matched_catalogue[object][filter_list[1]]['MAGERR_APER']) < 0.2:
            if float(matched_catalogue[object][filter_list[2]]['MAGERR_APER']) < 0.2:
                clean_catalogue[object] = matched_catalogue[object]

with open('%s/dual_clean_catalogue.pickle' %dict_dir, 'wb') as handle:
    pkl.dump(clean_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

#########
# STARS #
#########
star_catalogue = {}
for object in clean_catalogue.keys():

    if float(clean_catalogue[object][filter_list[0]]['CLASS_STAR']) >= 0.7:
        if float(clean_catalogue[object][filter_list[1]]['CLASS_STAR']) >= 0.7:
            if float(clean_catalogue[object][filter_list[2]]['CLASS_STAR']) >= 0.7:
                star_catalogue[object] = clean_catalogue[object]

with open('%s/dual_star_catalogue.pickle' %dict_dir, 'wb') as handle:
    pkl.dump(star_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

#########
# EDGES #
#########
# function that decides if a given cetnroid is in the edge regions
def mask(x_check, y_check, tolorence, image):

    # 0, if out of range
    shape_y, shape_x = np.shape(image)
    if x_check <= 0+tolorence or x_check >= shape_x-tolorence: return 0
    if y_check <= 0+tolorence or y_check >= shape_y-tolorence: return 0

    # 0, if 0
    if image[y_check, x_check] == 0.0: return 0

    # 0, if close to a zero patch
    for t in range(tolorence-2):
        if image[y_check, x_check-t] == 0.0 and image[y_check, x_check-t-1] == 0.0 and image[y_check, x_check-t-2] == 0.0:
            return 0
        if image[y_check, x_check+t] == 0.0 and image[y_check, x_check+t+1] == 0.0 and image[y_check, x_check+t+2] == 0.0:
            return 0

        if image[y_check-t, x_check] == 0.0 and image[y_check-t-1, x_check] == 0.0 and image[y_check-t-2, x_check] == 0.0:
            return 0
        if image[y_check+t, x_check] == 0.0 and image[y_check+t+1, x_check] == 0.0 and image[y_check+t+2, x_check] == 0.0:
            return 0

    # else
    return 1
#mask(1143, 1084, 20, reference_image)

# function that uses the above function to decide which clean_catalogue enteries to discard
def edger(dictionary, tolorence, filter_id):

    reference_filter = filter_list[filter_id]

    reference_image = fits.open(sci_array[filter_id])[0].data
    print('\n + using the image in %s as the reference' %sci_array[filter_id])
    print(' + using the filter %s as the reference' %reference_filter)
    print('\nInput number: %s' %len(dictionary))

    discarded_count = 0
    output_dictionary = {}
    for object_id in dictionary.keys():
        x_centroid = int(float(dictionary[object_id][reference_filter]['X_IMAGE']))
        y_centroid = int(float(dictionary[object_id][reference_filter]['Y_IMAGE']))

        decision = mask(x_centroid, y_centroid, tolorence, reference_image)
        if decision == 1:
            output_dictionary[object_id] = dictionary[object_id]
        if decision == 0:
            discarded_count = discarded_count + 1

    print('Discarded: %s' %discarded_count)
    return output_dictionary

# function that puts the edger function from above in a loop to be applied to more than one filter
def edger_multiBand(dictionary, tolorence, filter_ids):

    for filter_id in filter_ids:
        edged_catalogue = edger(dictionary, tolorence, filter_id)
        dictionary = edged_catalogue

    return dictionary

# producing the edged catalog
edged_catalogue = edger_multiBand(clean_catalogue, 35, [0,1,2])
with open('%s/dual_edged_catalogue.pickle' %dict_dir, 'wb') as handle:
    pkl.dump(edged_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

edged_star_catalogue = edger_multiBand(star_catalogue, 30, [0,1,2])
with open('%s/dual_edged_star_catalogue.pickle' %dict_dir, 'wb') as handle:
    pkl.dump(edged_star_catalogue, handle, protocol=pkl.HIGHEST_PROTOCOL)

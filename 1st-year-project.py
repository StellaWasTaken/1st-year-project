import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt

dict_directory = '/home/stella/1st-year-project/dict'

filter_list = ['f444w_2','f356w','f277w', "f770w", "f1000w", "f1500w", "f1800w"]

#dictionary for translation id's
trans_dir = {"f444w_2": 377, "f356w": 376, "f277w": 375, "f770w": 396, "f1000w": 397, "f1500w": 400, "f1800w": 401} #f770w uses f0770w data

with open('\\\\wsl.localhost\\Ubuntu\\home\\stella\\1st-year-project\\dict\\all_catalogue.pickle', 'rb') as handle:
    all_catalogue = pkl.load(handle)

print(all_catalogue[2].keys())

FluxData = []
FluxErr = []
for i in range(len(filter_list)):
    FluxData.append([])
    FluxErr.append([])

for key in all_catalogue.keys():
    for i in range(len(filter_list)):
        filter = filter_list[i]
        flux = float(all_catalogue[key][filter]["MAG_APER"]) #"FLUX_AUTO"
        err = float(all_catalogue[key][filter]["MAGERR_APER"]) #"FLUXERR_AUTO"
        FluxData[i].append(flux)
        FluxErr[i].append(err)
        print(all_catalogue[key][filter])
        if "spec_z" in all_catalogue[key][filter]:
            print()

#Header
main_txt = "# id"
line2 = "# id"
trans_txt = ""
for filter in filter_list:
    trans_index = str(trans_dir[filter])
    main_txt += " f_" + filter + " e_" + filter
    line2 += " F" + trans_index + " E" + trans_index
    #translation file
    if trans_txt != "":
        trans_txt += "\n"
    trans_txt += "f_" + filter + " F" + trans_index + "\ne_" + filter + " E" + trans_index
main_txt += "\n" + line2
#Data
for i in range(len(FluxData[0])): 
    main_txt += "\n" + str(i+1).rjust(6)
    for j in range(len(filter_list)):
        main_txt += f"  {FluxData[j][i]}  {FluxErr[j][i]}"

cat_dir = "C:\\Users\\Bruger\\Studie\\Data_og_projekt\\Projekt"
f = open(cat_dir + "\\file.cat", "w")
f.writelines(main_txt)
f.close()

f = open(cat_dir + "\\translation.cat", "w")
f.writelines(trans_txt)
f.close()

FluxData = np.array(FluxData)
FluxErr = np.array(FluxErr)
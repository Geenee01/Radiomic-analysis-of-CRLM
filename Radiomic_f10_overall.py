## make csvs with 10 selected features by using mrmr selecetion for hepatic recurrence

import numpy as np
import pandas as pd
import os

# defining Base path, input path for mrmr seleceted features and output paths
Base_path = "path/to/saved-data"
Base_path_mrmr = "path/to/mrmr-data"
Base_path_o ="path/to save/output/files"

# loading and reading necessary files and dictionaries
Clinical_df = pd.read_excel(Base_path + "clinical.xlsx")
totalL1_vol_dic = np.load(Base_path + 'totalL1voldic.npy', allow_pickle=True).item()



### UWA ###

# load radiomic features data
UWA_df = pd.read_csv(Base_path + "UWA_aggregation.csv", index_col=0)
#print(UWA_df)
# load mrmr data
UWA_mrmr_indices = pd.read_csv(Base_path_mrmr + "UWA_10features_mRMR.csv")

# making list of column indices
UWA_mrmr_indices = UWA_mrmr_indices['x'].tolist()

# change "." to "_", as the column names in the radiomic data have "-" instead of "." as per 'R' language
UWA_mrmr_indices = ['_'.join(feature.split('.')) for feature in UWA_mrmr_indices]

# create dataframe with top 10 feature values
#UWA_df.columns = UWA_df.columns.astype(int)
UWA_mRMRe = UWA_df[ UWA_mrmr_indices].copy()
print("UWA features:",UWA_mRMRe.columns)

# adding 'Time' and 'Event' columns accoding to index column

UWA_mRMRe["Time"]=Clinical_df["overall surv months"].astype(float)
UWA_mRMRe["Event"]=Clinical_df["fu status"].astype(int)
print(UWA_mRMRe)

# Save the selected data to a new CSV
UWA_mRMRe.to_csv(Base_path_o + 'UWA_f10.csv', index=False)


### WA ###

# load data
WA_df = pd.read_csv(Base_path + "WA_aggregation.csv", index_col=0)
WA_mrmr_indices = pd.read_csv(Base_path_mrmr + "WA_10features_mRMR.csv")

# making list of column indices
WA_mrmr_indices = WA_mrmr_indices['x'].tolist()

# change "." to "_", as the column names in the radiomic data have "-" instead of "." as per 'R' language
WA_mrmr_indices = ['_'.join(feature.split('.')) for feature in WA_mrmr_indices]

# create dataframe with top 10 feature values
WA_mRMRe = WA_df[WA_mrmr_indices].copy()
print("WA features:",WA_mRMRe.columns)

# adding 'Time' and 'Event' columns by mapping the values from clinical data accordning to patient ID
WA_mRMRe["Time"]=WA_mRMRe.index.to_series().map(Clinical_df.set_index('Patient-ID')['overall surv months'])
WA_mRMRe["Event"]=WA_mRMRe.index.to_series().map(Clinical_df.set_index('Patient-ID')['fu status'].astype(int))

# adding 'Time' and 'Event' columns
WA_mRMRe["Time"]=Clinical_df["overall surv months"]
WA_mRMRe["Event"]=Clinical_df["fu status"].astype(int)
#print(WA_mRMRe)

# Save the selected data to a new CSV
WA_mRMRe.to_csv(Base_path_o + 'WA_f10.csv', index=False)



### largest3 ###

# load data
largest3_df = pd.read_csv(Base_path + "largest3_aggregation.csv", index_col=0)
largest3_mrmr_indices = pd.read_csv(Base_path_mrmr + "largest3_10features_mRMR.csv")

# making list of column indices
largest3_mrmr_indices = largest3_mrmr_indices['x'].tolist()

# change "." to "_", as the column names in the radiomic data have "-" instead of "." as per 'R' language
largest3_mrmr_indices = ['_'.join(feature.split('.')) for feature in largest3_mrmr_indices]

# create dataframe with top 10 feature values
#UWA_df.columns = UWA_df.columns.astype(int)
largest3_mRMRe = largest3_df[largest3_mrmr_indices].copy()
print("largest3 features:",largest3_mRMRe.columns)

# adding 'Time' and 'Event' columns
largest3_mRMRe["Time"]=Clinical_df["overall surv months"]
largest3_mRMRe["Event"]=Clinical_df["fu status"].astype(int)

# Save the selected data to a new CSV
largest3_mRMRe.to_csv(Base_path_o + 'largest3_f10.csv', index=False)


### largest1, largest1 with number of tumors and largest 1 with total volume of tumors ### 

# load data
largest1_df = pd.read_csv(Base_path + "largest1_radiomic.csv", index_col=0)
largest1_mrmr_indices = pd.read_csv(Base_path_mrmr + "largest1_10features_mRMR.csv")

# making list of column indices
largest1_mrmr_indices = largest1_mrmr_indices['x'].tolist()

# change "." to "_", as the column names in the radiomic data have "." instead of "_" as per 'R' language
largest1_mrmr_indices = ['_'.join(feature.split('.')) for feature in largest1_mrmr_indices]

# create dataframe with top 10 feature values
largest1_mRMRe = largest1_df[largest1_mrmr_indices].copy()
print("largest1 features:",largest1_mRMRe.columns)

# adding 'Time' and 'Event' columns
largest1_mRMRe["Time"]=Clinical_df["overall surv months"]
largest1_mRMRe["Event"]=Clinical_df["fu status"].astype(int)

# adding 'tumor_num' column
largest1_num = largest1_mRMRe.copy()
#print(largest1_num)

# a folder where all the CT scan images volume and segments are stored
data_folders = "path/to/data"

# to create list of number of tumors in each patient file
tumor_num_list=[]

# looping through all the folders and then in files
for folder in os.listdir(data_folders):
    folder_path = os.path.join(data_folders, folder)
    if os.path.isdir(folder_path):
        tumor_counts_in_folder = 0        
        for file in os.listdir(folder_path):          
            if "Tumor" in file and ".mhd" in file:
                tumor_counts_in_folder+=1                
    tumor_num_list.append(tumor_counts_in_folder) 

# calculate median of the tumor counts
median_tumor_num = np.median(tumor_num_list)
print("Median tumor count: ", median_tumor_num)      
                            
# adding 'tumor_num' column
largest1_num["Tumor num"] = tumor_num_list

# Save the selected data to a new CSV
largest1_num.to_csv(Base_path_o + 'largest1_f10_num.csv', index=False)

#Add another column for 'total volume' and save in seperate csv file
largest1_vol = largest1_mRMRe.copy()
largest1_vol['total_volume'] = largest1_vol.index.to_series().map(totalL1_vol_dic)

# Save the selected data to a new CSV
largest1_vol.to_csv(Base_path_o + 'largest1_f10_vol.csv', index=False)


### smallest1 ###

# load data
smallest1_df = pd.read_csv(Base_path + "smallest1_radiomic.csv", index_col=0)
smallest1_mrmr_indices = pd.read_csv(Base_path_mrmr + "smallest1_10features_mRMR.csv")

# making list of column indices
smallest1_mrmr_indices = smallest1_mrmr_indices['x'].tolist()

# change "." to "_", as the column names in the radiomic data have "-" instead of "." as per 'R' language
smallest1_mrmr_indices = ['_'.join(feature.split('.')) for feature in smallest1_mrmr_indices]

# create dataframe with top 10 feature values
smallest1_mRMRe = smallest1_df[smallest1_mrmr_indices].copy()
print("smallest1 features:",smallest1_mRMRe.columns)

# adding 'Time' and 'Event' columns
smallest1_mRMRe["Time"]=Clinical_df["overall surv months"]
smallest1_mRMRe["Event"]=Clinical_df["fu status"].astype(int)

# Save the selected data to a new CSV
smallest1_mRMRe.to_csv(Base_path_o + 'smallest1_f10.csv', index=False)
print("done")

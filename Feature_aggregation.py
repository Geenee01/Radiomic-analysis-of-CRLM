### Radiomic features aggregation ###

import numpy as np
import pandas as pd
import os
from sklearn.preprocessing import StandardScaler

# base path to folder
Base_path = "C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/"
 
# load radiomic data
radiomics_data = pd.read_csv(Base_path + "extracted_radiomic_features_0-300.csv")

# load clinical data
clinical_data_df = pd.read_excel(Base_path + "CRLM Clinical data.xlsx")

# create list of columns starting with 'original' (copied from brain_metastases)
value_feature_names = [col for col in radiomics_data.columns if (col.startswith('original'))] + [col for col in radiomics_data.columns if (col.startswith('wavelet'))]

# Creating new df "radiomics_d" to store patient ID and tumor Id with features' name
radiomics_d = radiomics_data[['patient ID', 'Tumor'] + value_feature_names]

# clean column names, removing 'original'
radiomics_d.columns = ['patient ID', 'Tumor'] + ['_'.join(col.split('-')) for col in value_feature_names]

# keep names of radiomic features only
clean_col_names = radiomics_d.columns.tolist()[3:]

# Common among all methods: sum of shape features and average of other radiomic features
# save "shape" features in a list
shape_features = [col for col in radiomics_d.columns if 'shape' in col]

# save "non-shape" features in a list which also includes "patient ID" and "Tumor"
non_shape_features = [col for col in radiomics_d.columns if 'shape' not in col]

# initializing standard scaler for normalization
scaler = StandardScaler()


### UWA ###

# sum of "shape" features
columns_to_sum = ["patient ID"] + shape_features
sum_df = radiomics_d[columns_to_sum].copy()
sum_df = sum_df.groupby('patient ID').sum()

# Mean of other features:
# 'patient ID' and 'Tumor' do not have 'shape' in their names, so they are already included
mean_df = radiomics_d[non_shape_features].copy()

# Exclude non-numeric column "Tumor"from mean_df
mean_df = mean_df.drop(['Tumor'], axis=1)

# Group by 'patient ID' and calculate the mean
mean_df = mean_df.groupby('patient ID').mean()

# Reset the index to make 'patient ID' a regular column
mean_df.reset_index(inplace=True)

# combine sum and mean dfs
UWA_df = pd.merge(sum_df, mean_df, on = "patient ID")

# reset the index to make 'patient ID' a index column
UWA_df.set_index('patient ID', inplace=True)

# normalization
scaler = StandardScaler()
UWA_df.iloc[:,1 :] = scaler.fit_transform(UWA_df.iloc[:, 1:])

# change the index column of UWA_df as clinical data's 'Patient-ID' column values, if UWA_df index coumns values are same as clinical df's De-identify Scout Name values
UWA_df.index = UWA_df.index.map(clinical_data_df.set_index('De-identify Scout Name')['Patient-ID'])

# reset the index
UWA_df.reset_index(inplace=True)

# drop the 'patient ID' column (for mRMRe)
UWA_df = UWA_df.drop('patient ID', axis=1)

# save radiomic features for feature selection
UWA_radiomic_Ofile = Base_path + 'UWA_radiomic_aggregation.csv'
UWA_df.to_csv(UWA_radiomic_Ofile, index=True)

## restore dictonary##
# restore column names
#inv_col_names = {v: k for k, v in col_names.items()}


# save the dictionary
#np.save(Base_path + 'colnames.npy', inv_col_names, allow_pickle=True)
#print(inv_col_names_with_start_index_1)


### WA ###

# radiomic_d copy to new df
rad_df = radiomics_d.copy()
# calculating total volume per patient
total_volume = rad_df[['patient ID', 'original_shape_MeshVolume']].copy()
total_volume = total_volume.groupby('patient ID').sum()

# save labels mapped to ID
total_vol_dic = pd.Series(total_volume['original_shape_MeshVolume'].values, index=total_volume.index).to_dict()
np.save(Base_path + 'totalvoldic.npy', total_vol_dic)

# map labels to radiomic data
rad_df['total volume'] = rad_df['patient ID'].map(total_vol_dic)

# calculate weight
rad_df['weight'] = rad_df['original_shape_MeshVolume'] / rad_df['total volume']

# store sum of "shape" features
columns_to_sum = ["patient ID"] + shape_features
sum_df = rad_df[columns_to_sum].copy()
sum_df = sum_df.groupby('patient ID').sum()

# 'patient ID' and 'Tumor' do not have 'shape' in their names, so they are already included
w_ave_df = rad_df[(non_shape_features) + ['weight']].copy()

# Exclude column "Tumor" from w_ave_df
w_ave_df = w_ave_df.drop(['Tumor'], axis=1)

# Setting "patient ID" as an index
w_ave_df = w_ave_df.set_index('patient ID')

# multiply all values by weight
w_ave_df = w_ave_df.iloc[:, 0:837 ].multiply(w_ave_df['weight'], axis='index')

# sum the weighted values to get weighted average
w_ave_df= w_ave_df.groupby("patient ID").sum()

# Reset the index to make 'patient ID' a regular column
w_ave_df.reset_index(inplace=True)

# combine sum and mean df
WA_df = pd.merge(sum_df, w_ave_df, on = "patient ID")

# set 'patient ID' as index
WA_df.set_index('patient ID', inplace=True)

# standardize
#WA_df.columns = WA_df.columns.astype(int)
scaler = StandardScaler()
WA_df.iloc[:, 1:] = scaler.fit_transform(WA_df.iloc[:, 1:])

# change the index column of UWA_df as clinical data's 'Patient-ID' column values, if WA_df index columns values are same as clinical df's De-identify Scout Name values
WA_df.index = WA_df.index.map(clinical_data_df.set_index('De-identify Scout Name')['Patient-ID'])

# reset the index
WA_df.reset_index(inplace=True)

# drop the 'patient ID' column (for mRMRe)
WA_df = WA_df.drop('patient ID', axis=1)

# save radiomic features for feature selection
WA_radiomic_Ofile = Base_path + 'WA_radiomic_aggregation.csv'
WA_df.to_csv(WA_radiomic_Ofile, index=True)


### WA of largest 3 ###

large_3 = radiomics_d.copy()

large_3 = large_3.sort_values(by=['patient ID', 'original_shape_MeshVolume'], ascending=False).groupby('patient ID').head(3)

# store sum of "shape" features
columns_to_sum = ["patient ID"] + shape_features
sum_df = large_3[columns_to_sum].copy()
sum_df = sum_df.groupby('patient ID').sum()

# calculating weighted average of 3 largest
# calculate total volume of largest 3
large_3_total_vol = large_3[['patient ID', 'original_shape_MeshVolume']].copy()
large_3_total_vol = large_3.groupby('patient ID').sum()
 
# save labels mapped to ID
large_3_t_vol_dic = pd.Series(large_3_total_vol['original_shape_MeshVolume'].values, index=large_3_total_vol.index).to_dict()

# maps labels to radiomic data
large_3['Largest 3 total volume'] = large_3['patient ID'].map(large_3_t_vol_dic)

# calculate weight
large_3['largest 3 weight'] = large_3['original_shape_MeshVolume'] / large_3['Largest 3 total volume']

# 'patient ID' and 'Tumor' do not have 'shape' in their names, so they are already included
large_3_wa_df = large_3[(non_shape_features) + ['largest 3 weight']].copy()

# Exclude column "Tumor" from w_ave_df
large_3_wa_df = large_3_wa_df.drop(['Tumor'], axis=1)

# Setting "patient ID" as an index
large_3_wa_df = large_3_wa_df.set_index('patient ID')

# multiply all values by weight
large_3_wa_df = large_3_wa_df.iloc[:, 0:837 ].multiply(large_3_wa_df['largest 3 weight'], axis='index')

# sum the weighted values to get weighted average
large_3_wa_df= large_3_wa_df.groupby("patient ID").sum()

# Reset the index to make 'patient ID' a regular column
large_3_wa_df.reset_index(inplace=True)

# combine sum and mean df
largest_3_df = pd.merge(sum_df, large_3_wa_df, on = 'patient ID')

# set 'patient ID' as index
largest_3_df.set_index('patient ID', inplace=True)

# standardize
#largest_3_df.columns = largest_3_df.columns.astype(int)
scaler = StandardScaler()
largest_3_df.iloc[:, 1:] = scaler.fit_transform(largest_3_df.iloc[:,1:])

# change the index column of UWA_df as clinical data's 'Patient-ID' column values, if WA_df index columns values are same as clinical df's De-identify Scout Name values
largest_3_df.index = largest_3_df.index.map(clinical_data_df.set_index('De-identify Scout Name')['Patient-ID'])

# reset the index
largest_3_df.reset_index(inplace=True)

# drop the 'patient ID' column (for mRMRe)
largest_3_df = largest_3_df.drop('patient ID', axis=1)

# save radiomic features for feature selection
largest3_radiomic_Ofile = Base_path + 'largest3_radiomic_aggregation.csv'
largest_3_df.to_csv(largest3_radiomic_Ofile, index=True)



### Largest 1 ###

# copying the data
largest_1_df = radiomics_d.copy()

# sorting the largest tumor by mesh volume
largest_1_df = largest_1_df.sort_values(by=['patient ID', 'original_shape_MeshVolume'], ascending=False).groupby('patient ID').head(1)

print(largest_1_df['original_shape_MeshVolume'].describe())
print(np.percentile(largest_1_df['original_shape_MeshVolume'], 33))
print(np.percentile(largest_1_df['original_shape_MeshVolume'], 50))
print(np.percentile(largest_1_df['original_shape_MeshVolume'], 66))

# making df ascending for patient IDs 
largest_1_df = largest_1_df.sort_values(by="patient ID", ascending=True)

# reset the index
largest_1_df.reset_index(inplace=True)

# drop the 'patient ID' column (for mRMRe)
new_largest_1_df = largest_1_df.drop(["index","Tumor"], axis=1)

# set 'patient ID' as index
new_largest_1_df.set_index('patient ID', inplace=True)

# standardize
#new_largest_1_df.columns = new_largest_1_df.columns.astype(int)
scaler = StandardScaler()
new_largest_1_df.iloc[:, 1:] = scaler.fit_transform(new_largest_1_df.iloc[:, 1:])

# change the index column of UWA_df as clinical data's 'Patient-ID' column values, if WA_df index columns values are same as clinical df's De-identify Scout Name values
new_largest_1_df.index = new_largest_1_df.index.map(clinical_data_df.set_index('De-identify Scout Name')['Patient-ID'])

# reset the index
new_largest_1_df.reset_index(inplace=True)

# drop the 'patient ID' column (for mRMRe)
new_largest_1_df = new_largest_1_df.drop('patient ID', axis=1)

# save 'total_volume' column for 'largest 1 and total volume of tumors)
total_volume = radiomics_d[['patient ID', 'original_shape_MeshVolume']].copy()
total_volume = total_volume.groupby('patient ID').sum()

# save 'total_volume' dictionary
totalL1_vol_dic = pd.Series(total_volume['original_shape_MeshVolume'].values, new_largest_1_df.index).to_dict()
np.save(Base_path + 'totalL1voldic.npy', totalL1_vol_dic)


# path to save csv file
largest1_ts_radiomic_Ofile = Base_path + 'largest1_radiomic.csv'

# save to csv
new_largest_1_df.to_csv(largest1_ts_radiomic_Ofile, index=True)


### smallest 1 ###

smallest_1_df = radiomics_d.copy()
#print(smallest_1_df)
# sorting the smallest tumor by mesh volume
smallest_1_df = smallest_1_df.sort_values(by=['patient ID', 'original_shape_MeshVolume'], ascending=True).groupby('patient ID').head(1)

# making df ascending for patient IDs
smallest_1_df = smallest_1_df.sort_values(by=['patient ID'], ascending=True)

# reset the inde
smallest_1_df.reset_index(inplace=True)

# drop the 'patient ID' column (for mRMRe)
smallest_1_df = smallest_1_df.drop(['index', 'Tumor'], axis=1)

# set 'patient ID' as index
smallest_1_df.set_index('patient ID', inplace=True)

# standardize
scaler = StandardScaler()
smallest_1_df.iloc[:, 1:] = scaler.fit_transform(smallest_1_df.iloc[:, 1:])

# change the index column of UWA_df as clinical data's 'Patient-ID' column values, if WA_df index columns values are same as clinical df's De-identify Scout Name values
smallest_1_df.index = smallest_1_df.index.map(clinical_data_df.set_index('De-identify Scout Name')['Patient-ID'])

# reset the index
smallest_1_df.reset_index(inplace=True)

# drop the 'patient ID' column (for mRMRe)
smallest_1_df = smallest_1_df.drop('patient ID', axis=1)

# path to save csv file
smallest1_radiomic_Ofile = Base_path + 'smallest1_radiomic.csv'

# save to csv
smallest_1_df.to_csv(smallest1_radiomic_Ofile, index=True) 
print("done")

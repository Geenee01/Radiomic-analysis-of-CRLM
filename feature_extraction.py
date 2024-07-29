# Importing libraries
import os
import logging
import SimpleITK as sitk
import radiomics
import numpy as np
import pandas as pd
from radiomics import featureextractor


## Logging is straight from helloRadiomics.ipynb from pyradiomics ##
# Get the PyRadiomics logger (default log-level = INFO)
logger = radiomics.logger
logger.setLevel(logging.DEBUG)  # set level to DEBUG to include debug log messages in log file

# Write out all log entries to a file
handler = logging.FileHandler(filename='testLog.txt', mode='w')
formatter = logging.Formatter('%(levelname)s:%(name)s: %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

Base_path = "path/to/data"
Nii_path = Base_path 

# Folder path
data_folders = Base_path + "MCRC"

# Initialize the PyRadiomics feature extractor
params = "path/to/Params.yml"
extractor = featureextractor.RadiomicsFeatureExtractor(params)

# An empty DataFrame 
radiomic_data = pd.DataFrame()

# Iterate through the folders in the data directory
for folder in os.listdir(data_folders):
    folder_path = os.path.join(data_folders, folder)
    if os.path.isdir(folder_path):

        # An empty list to store extracted features from this folder
        folder_features = []
       
        ct_scan_image = None
        tumor_segmentation = None
        tumor_list=[]

        # To find the total CT scan image and tumor segment files in the folder
        for filename in os.listdir(folder_path):
            
            if "volume.mhd" in filename:

                # Read the '.mhd' files
                ct_scan_image = sitk.ReadImage(os.path.join(folder_path, filename))

            elif "Tumor" in filename and ".mhd" in filename:
                tumor_list.append(filename)
        
         # Iterating through each tumor segment to extract radiomic features    
        for tumor in tumor_list:
                
                # Read the '.mhd' files
                tumor_seg = sitk.ReadImage(os.path.join(folder_path, tumor))
                tumor_segmentation = (tumor_seg != -1000)

                # Extracting  radiomic features
                feature_vector = extractor.execute(ct_scan_image, tumor_segmentation)                      
                folder_features.append(feature_vector)
     

        # Create a DataFrame from the folder_features list
        # Adding patinet_ID and Tumor file in the dataframe
        if folder_features:       
            folder_data = pd.DataFrame(folder_features)
            folder_data.insert(0, "patient ID", folder)
            folder_data.insert(1, "Tumor", tumor_list)
            
            # Remove ".nii.gz" extension from Tumor filenames
            folder_data['Tumor'] = folder_data['Tumor'].str.replace('.mhd', '')
            radiomic_data = pd.concat([radiomic_data, folder_data], ignore_index=True)

# Save the radiomic features to a CSV file
output_csv_path = Base_path + "radiomic.csv"
radiomic_data.to_csv(output_csv_path, index=False)


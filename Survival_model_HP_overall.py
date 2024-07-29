### Creating and comparing 3 models ###

'''
# 3 Survival models:
(1) Cox Proposional Hazard
(2) Cox Proposional Hazard Lasso Regression
(3) Random Survival Forest

# 6 aggregation methods:
(1) Unweighted Average
(2) Weighted Average
(3) Weighted Average of Largest 3 tumors
(4) Largest 1 tumor
(5) Largest 1 tumor and number of metastases(tumor) in liver
(7) Largest 1 tumor and total volume of tumors in liver
(6) Smallest 1 tumor

'''

# import libraries
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.utils import resample
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest
from sklearn.model_selection import train_test_split, GridSearchCV
import seaborn as sns
import matplotlib.pyplot as plt

# Base path to selceted top10 feature csvs folder
Base_path = "path/to/features"

# csv files with 10 selected radiomic features and 'time' and 'event' data for whole dataset
# have different files for largest 1 with numbers of tumors and laregst 1 with total volume of tumors
UWA_top10 = Base_path + "UWA_f10.csv"
WA_top10 = Base_path + "WA_f10.csv"
largest3_top10 = Base_path + "largest3_f10.csv"
largest1_top10 = Base_path + "largest1_f10_num.csv"
largest1_top10_vol = Base_path + "largest1_f10_vol.csv"
smallest1_top10 = Base_path + "smallest1_f10.csv"


###### Models ######

### CPH Model ###
def CPH_Bootstrap(file_path, num=False, vol=False):
    '''
    This function will compute CPH with bootstrapping method.

    :param file_path: (str) filename of selected features
    :param num: (bool) set True to include the number of mets (col 13)
    :return: (str) C-index (95% confidence interval)
    '''

    df = pd.read_csv(file_path)
   
    
    # Configure bootstrap 
    num_iterations = 100
    num_size = int(len(df) * 0.50) # sampling 50% (98 samples) of data

    # Calculate population of statistics
    metrics = []
    for i in range(num_iterations):
        # Prepare the sample
        # If indicated, include the number of tumors or total volume of tumors (col 13)
        if num:
            sample = resample(df.iloc[:, np.r_[:13]], n_samples=num_size, random_state=i)           
        elif vol:
            sample = resample(df.iloc[:, np.r_[:13]], n_samples=num_size, random_state=i)
        else:
            sample = resample(df.iloc[:, np.r_[:12]], n_samples=num_size, random_state=i)
        
        # Calculate c-index and append to the list
        cph = CoxPHFitter(penalizer=0.001).fit(sample, 'Time', 'Event')  #penalizer=0.01
        score = concordance_index(sample['Time'], -cph.predict_partial_hazard(sample), sample['Event'])
        metrics.append(score)
        
    # Calculate confidence interval 
    alpha = 0.95
    p = ((1.0 - alpha) / 2.0) * 100
    lower = max(0.0, np.percentile(metrics, p)) #2.5th percentile
    p = (alpha + ((1.0 - alpha) / 2.0)) * 100
    upper = min(1.0, np.percentile(metrics, p)) #97.5th percentile
    med = np.percentile(metrics, 50)

    # Identify the aggregation method name
    if num:
        name = file_path.split('/')[-1].split('_')[0] + ' and number of Tumors'
    elif vol:
        name = file_path.split('/')[-1].split('_')[0] + ' and total volume of Tumors'
    else:
        name = file_path.split('/')[-1].split('_')[0]


    return metrics, print(name, 'CPH model', '%.3f (%.3f-%.3f)' % (med, lower, upper))   

       

### CPH Lasso regression modle ###

def Lasso_CPH_Bootstrap(file_path, num= False, vol=False):
    '''
    This function will compute CPH-Lasso regression with bootstrapping method. And will do Gridsearch hyperparameter tuning with 5 fold cross validation.

    :param file_path: (str) filename of selected features
    :param num: (bool) set True to include the number of mets (col 13)
    :return: (str) C-index (95% confidence interval)
    '''
    #seed_value = 98
    #np.random.seed(seed_value)
    df = pd.read_csv(file_path)
    
    # Configure bootstrap 
    num_iterations = 100
    num_size = int(len(df) * 0.50) # sampling 50% (98 samples) of data

    # Calculate population of statistics
    metrics = []
    for i in range(num_iterations):
        # Prepare the sample

        # If indicated, include the number of tumors (col 13)
        if num:
            sample = resample(df.iloc[:, np.r_[:13]], n_samples=num_size, random_state=i)
            X = sample.iloc[:, np.r_[:10, 12]].copy() #indices are different in 'sample' dataframe

        elif vol:
            sample = resample(df.iloc[:, np.r_[:13]], n_samples=num_size, random_state=i)
            X = sample.iloc[:, np.r_[:10, 12]].copy() #indices are different in 'sample' dataframe
                
        else:
            sample = resample(df.iloc[:, np.r_[:12]], n_samples=num_size, random_state=i)
            X = sample.iloc[:, np.r_[:, :10]].copy() #indices are different in 'sample' dataframe 
            
        X = X.to_numpy()
        y = sample[['Event', 'Time']].copy()
        y['Event'] = y['Event'].astype('bool')
        y = y.to_records(index=False)
                                                        
        grid_val = {
            'alphas': [[0.045], [0.050],[0.055],[0.060], [0.065],[0.070], [0.075],[0.080], [0.085],[0.090], [0.095]]
        }

        coxnet = CoxnetSurvivalAnalysis(l1_ratio=1)
        grid_search = GridSearchCV(estimator=coxnet, param_grid =grid_val, cv=5) 
        grid_search.fit(X, y)

        best_estimator = grid_search.best_estimator_
        score = best_estimator.score(X, y)
        metrics.append(score)
        
    # Calculate confidence interval 
    alpha = 0.95
    p = ((1.0 - alpha) / 2.0) * 100
    lower = max(0.0, np.percentile(metrics, p))
    p = (alpha + ((1.0 - alpha) / 2.0)) * 100
    upper = min(1.0, np.percentile(metrics, p))
    med = np.percentile(metrics, 50)

    # Identify the aggregation method name
    if num:
        name = file_path.split('/')[-1].split('_')[0] + ' and number of Tumors'
    elif vol:
        name = file_path.split('/')[-1].split('_')[0] + ' and total volume of Tumors'
    else:
        name = file_path.split('/')[-1].split('_')[0]

    return print(name, 'CPH-Lasso regression model', '%.3f (%.3f-%.3f)' % (med, lower, upper))


### RSF model ###

def RSF_Bootstrap(file_path, num=False, vol=False):
    '''
    This function will compute RSF with bootstrapping method. And will do Gridsearch hyperparameter tuning with 5 fold cross validation.

    :param file_path: (str) filename of selected features
    :param num: (bool) set True to include the number of mets (col 13)
    :return: (str) C-index (95% confidence interval)
    '''

    df = pd.read_csv(file_path)

    # configure bootstrap 
    num_iterations = 100
    num_size = int(len(df) * 0.50) # sampling 50% (98 samples) of data

	# parameters for RSF
    #NUMESTIMATORS = 100
    TESTSIZE = 0.20
    random_state = 20

    # Calculate population of statistics
    metrics = []
    for i in range(num_iterations):
        # Prepare the sample

        # If indicated, include the number of tumors (col 13)
        if num:
            sample = resample(df.iloc[:, np.r_[:13]], n_samples=num_size, random_state=i)
            X = sample.iloc[:, np.r_[:10, 12]].copy() #indices are different in 'sample' dataframe
        elif vol:
            sample = resample(df.iloc[:, np.r_[:13]], n_samples=num_size, random_state=i)
            X = sample.iloc[:, np.r_[:10, 12]].copy() #indices are different in 'sample' dataframe 
        else:
            sample = resample(df.iloc[:, np.r_[:12]], n_samples=num_size, random_state=i)
            X = sample.iloc[:, np.r_[:10]].copy() #indices are different in 'sample' dataframe 

        X = X.to_numpy().astype('float64')
        y = sample[['Event', 'Time']].copy()
        y['Event'] = y['Event'].astype('bool')
        y['Time'] = y['Time'].astype('float64')
        y = y.to_records(index=False)
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TESTSIZE, random_state=random_state)
        grid_val_rsf = {
            'n_estimators': [100],
            'min_samples_split': [10, 15, 20, 25],
            'min_samples_leaf': [5, 8, 10, 15],
            'max_features': ['sqrt']
        }

        rsf = RandomSurvivalForest(n_jobs=-1, random_state=random_state)
        grid_search = GridSearchCV(estimator=rsf, param_grid=grid_val_rsf, cv=5)
        grid_search.fit(X_train, y_train)
        
        best_estimator = grid_search.best_estimator_
        score = best_estimator.score(X_test, y_test)
        metrics.append(score)

    # Calculate confidence interval 
    alpha = 0.95
    p = ((1.0 - alpha) / 2.0) * 100
    lower = max(0.0, np.percentile(metrics, p))
    p = (alpha + ((1.0 - alpha) / 2.0)) * 100
    upper = min(1.0, np.percentile(metrics, p))
    med = np.percentile(metrics, 50)

    # Identify the aggregation method name
    if num:
        name = file_path.split('/')[-1].split('_')[0] + ' and Num of Tumors'
    elif vol:
        name = file_path.split('/')[-1].split('_')[0] + ' and total volume of Tumors'   
    else:                                                                         
        name = file_path.split('/')[-1].split('_')[0]

    return print(name, 'RSF model', '%.3f (%.3f-%.3f)' % (med, lower, upper))


### Applying and comparing models###


# Cox Proportional Hazards

CPH_Bootstrap(UWA_top10)
CPH_Bootstrap(WA_top10)
CPH_Bootstrap(largest3_top10)
CPH_Bootstrap(largest1_top10)
CPH_Bootstrap(largest1_top10, num=True)
CPH_Bootstrap(largest1_top10_vol, vol=True)
CPH_Bootstrap(smallest1_top10)


# Cox Proportional Hazards with Lasso Regularization

Lasso_CPH_Bootstrap(UWA_top10)
Lasso_CPH_Bootstrap(WA_top10)
Lasso_CPH_Bootstrap(largest3_top10)
Lasso_CPH_Bootstrap(largest1_top10) # need to remove alphas from [0.005] to [0.045]
Lasso_CPH_Bootstrap(largest1_top10, num=True)
Lasso_CPH_Bootstrap(largest1_top10_vol, vol=True)
Lasso_CPH_Bootstrap(smallest1_top10)


# RSF: Random Survival Forest
RSF_Bootstrap(UWA_top10)
RSF_Bootstrap(WA_top10)
RSF_Bootstrap(largest3_top10)
RSF_Bootstrap(largest1_top10)
RSF_Bootstrap(largest1_top10, num=True)
RSF_Bootstrap(largest1_top10_vol, vol=True)
RSF_Bootstrap(smallest1_top10)

# visualize the results for comparison
metrics_UWA, cph = CPH_Bootstrap(UWA_top10)
metrics_WA, cph = CPH_Bootstrap(WA_top10)
metrics_largest3, cph = CPH_Bootstrap(largest3_top10)
metrics_largest1, cph = CPH_Bootstrap(largest1_top10)
metrics_largest1_num, cph = CPH_Bootstrap(largest1_top10, num=True)
metrics_largest1_vol, cph = CPH_Bootstrap(largest1_top10_vol, vol=True)
metrics_smallest1, cph = CPH_Bootstrap(smallest1_top10)

# Create a dataframe with all the metrics
data = {'UWA': metrics_UWA, 'WA': metrics_WA, 'Largest3': metrics_largest3, 'Largest1': metrics_largest1, 'Largest1_num': metrics_largest1_num, 'Largest1_vol': metrics_largest1_vol, 'Smallest1': metrics_smallest1}

df = pd.DataFrame(data)
df.to_csv(Base_path + 'CPH_metrics.csv', index=False)

data = pd.read_csv(Base_path + 'CPH_metrics.csv')
df_for_plt ={ 'Unweighted average': data['UWA'], 'Weighted average': data['WA'], 'Largest3': data['Largest3'], 'Largest1_vol': data['Largest1_vol'], 'Largest1_num': data['Largest1_num'], 'Largest1': data['Largest1'],'Smallest1': data['Smallest1']}

#Plot violin plot
#plt.axvline(x=uni_dict[dataName],linestyle='--',color='black')
sns.violinplot(data=df_for_plt,orient='h',palette='dark')
plt.xlabel('Concordance Index (C-Index)')
plt.ylabel('Method')
plt.title('CPH Model for Overall Survival with single tumors(197 patients)')
plt.savefig('testing.png',dpi=200,bbox_inches='tight')
plt.show()


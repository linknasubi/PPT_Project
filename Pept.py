import pandas as pd 
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA, KernelPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.manifold import TSNE
from sklearn.feature_selection import SelectKBest, chi2, f_classif


dataSet = ['Dataset/df_cppsite.csv', 'Dataset/df_pc2pred.csv', 'Dataset/df_nc2pred.csv',
          'Dataset/df_pselect1.csv', 'Dataset/df_nselect1.csv', 'Dataset/df_pselect2.csv',
          'Dataset/df_nselect2.csv', 'Dataset/seq_cppsite.csv', 'Dataset/seq_pc2pred.csv',
          'Dataset/seq_nc2pred.csv', 'Dataset/seq_pselect1.csv', 'Dataset/seq_nselect1.csv',
          'Dataset/seq_pselect2.csv', 'Dataset/seq_nselect2.csv']





def plotMaker(data: pd.DataFrame, label_data: str) -> None:
    lista1 = []
    
    for _ in range (len(data[label_data])):
        lista1.append(_)

    
    plt.figure(1)
    data[label_data].plot(kind='bar')
    plt.show()
    

plotMaker(dic_data['df_cppsite'], 'NHA')


def sampleMaker(data: pd.DataFrame) -> int:
    
    std_dev_list = []
    scaler = MinMaxScaler()
    
    for _ in range (1, data.shape[1]-1):
        
        std_aux_list = data[data.columns[_]].std()
        
        
        aux_list = (data[data.columns[_]] - std_aux_list/(2*std_aux_list))
        
#        aux_list = ((data[data.columns[_]] - data[data.columns[_]].max())/(data[data.columns[_]].max() - data[data.columns[_]].min()))
        
        mu, std = norm.fit(data[data.columns[_]])
        
        std_dev_list.append(aux_list.std(axis=0))
        
        
        plt.figure(1)
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 200)
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, 'k', linewidth=2)
        title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
        plt.title(title)
#        aux_list.plot(kind='bar')
        plt.show()
        
    
    
    print(std_dev_list)

"""[0.12262558166225068, 0.12292542572090273, 0.15021475920771704, 
0.12577227285815143, 0.1340860140623043, 0.1283058902867711, 0.
12364030536140189, 0.15835296459736659, 0.2170518093063288, 0.14310368662148973]"""


"""

[897.0519918269031, 414.29310379864535, 10.57415767376415, 12.828771831531448, 16.76075175778803,
 2.566117805735418, 32.27011969932588, 0.09765099483504283, 0.21705180930632875, 0.1431036866214896]
            """
sampleMaker(dic_data['df_cppsite'])


def prepareData(data_Set: list) -> dict:
    
    
    dic_data = {}
    form_dic = {}
    
    """
    
        A iteração abaixo serve para associar cada path do CSV a uma chave correspondente ao próprio nome do CSV.
    
    
    """
    
    for _ in range(len(dataSet)):
        dic_data[dataSet[_][8:-4]] = pd.read_csv(dataSet[_])
        
    
    """
    
        A iteração abaixo faz todo tratamento para cada Data Set conforme demonstrado no Jupyter


    """
    
    for _ in range(len(dataSet)):
        
        dic = dataSet[_]
        if dic[8:10] == 'df':
        
                 
            
            if dic[11:13] == 'pc' or dic[11:13] == 'ps':
                
                print(dic[8:-4])
                
                dic_data[dic[8:-4]] = dic_data[dic[8:-4]][dic_data[dic[8:-4]]['total[AA]'] < 31]
    
                
                try:
                
                    dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop(columns=["total[AA]","Unnamed: 0"])
                
                except KeyError:
                    
                    dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop(columns="total[AA]")
                
                
                dic_data[dic[8:-4]].insert(column='Class', value=np.ones(dic_data[dic[8:-4]].shape[0]).T, 
                 allow_duplicates=True, loc=10)
                
                dic_data[dic[8:-4]].insert(column='Seq', value= dic_data['seq'+ dic[10:-4]]['Seq'],
                    allow_duplicates=True, loc=0)
                
                dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop_duplicates(subset='Seq')
            
            elif dic[11:13] == 'nc' or dic[11:13] == 'ns':
                
                
                dic_data[dic[8:-4]] = dic_data[dic[8:-4]][dic_data[dic[8:-4]]['total[AA]'] < 31]
                
                try:
                
                    dic_data[dic[8:-4]] = dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop(columns=["total[AA]","Unnamed: 0"])
                
                except KeyError:
                    
                    dic_data[dic[8:-4]] = dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop(columns="total[AA]")
                
                dic_data[dic[8:-4]].insert(column='Class', value=np.zeros(dic_data[dic[8:-4]].shape[0]).T, 
                 allow_duplicates=True, loc=10)
                
                dic_data[dic[8:-4]].insert(column='Seq', value= dic_data['seq'+ dic[10:-4]]['Seq'],
                    allow_duplicates=True, loc=0)
                
                dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop_duplicates(subset='Seq')
            
            else:

                dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop(columns="total[AA]")
                dic_data[dic[8:-4]].insert(column='Class', value=np.ones(dic_data[dic[8:-4]].shape[0]).T, 
                    allow_duplicates=True, loc=10)
                dic_data[dic[8:-4]].insert(column='Seq', value= dic_data['seq'+ dic[10:-4]]['Seq'],
                    allow_duplicates=True, loc=0)
                dic_data[dic[8:-4]] = dic_data[dic[8:-4]].drop_duplicates(subset='Seq')
            
            
        form_dic[dic[8:-4]] = dic_data[dic[8:-4]].shape
            

            
    return dic_data, form_dic

"""
    
    É criado dois dicionários, sendo o primeiro referente ao DataFrame, e o segundo
    referente ao formato dos DataFrames.

"""


dic_all = prepareData(dataSet)
dic_data = dic_all[0]
form_data = dic_all[1]


neg_pep1 = pd.concat([dic_data['df_nc2pred'], dic_data['df_nselect2']])
##aux = pd.concat([dic_data['df_cppsite'], dic_data['df_pc2pred']]).drop_duplicates(subset='Seq')
#pos_pep1 = pd.concat([dic_data['df_pc2pred'], dic_data['df_pselect2']]) #focar no sample

#neg_pep1 = pd.concat([df_nc2pred,df_nselect2])
aux = pd.concat([dic_data['df_cppsite'], dic_data['df_pc2pred']]).drop_duplicates(subset='Seq')
pos_pep1 = pd.concat([dic_data['df_pselect2'],aux.sample(n=369,random_state=1)])




data_train1 = pd.concat([pos_pep1,neg_pep1]) 
    
    
all_pep1 = pd.concat([dic_data['df_cppsite'],dic_data['df_pc2pred'], dic_data['df_pselect2'], dic_data['df_nc2pred'],
                      dic_data['df_nselect2']])

    
    
    
    

    
X_train1 = data_train1.drop(columns=["Class","Seq"]).values
y_train1 = data_train1["Class"].values

X_all1 = all_pep1.drop(columns=["Class","Seq"]).values
y_all1 = all_pep1["Class"].values




# Data Teste II

data_test2 = pd.concat([dic_data['df_pselect1'], dic_data['df_nselect1']])

X_test2 = data_test2.drop(columns=["Class","Seq"]).values
y_test2 = data_test2["Class"].values


# ANOVA teste

f1 = SelectKBest(f_classif, k=2)
f2 = SelectKBest(f_classif, k=2)
f3 = SelectKBest(f_classif, k=2)

f1.fit_transform(X_all1, y_all1)
f2.fit_transform(X_train1, y_train1)
f3.fit_transform(X_test2, y_test2)

y1 = f1.scores_
y2 = f2.scores_
y3 = f3.scores_

x1=np.array([['MolWt'], ['TPSA'], ['LogP'], ['NHA'], ['NHD'], ['NAR'], ['NRB'], ['fcSP3'], ['f[ARG]'],['f[LYS]']]).ravel()

figA1 = plt.figure(figsize=(20,5))

fig1 = plt.figure(figsize=(20,20))




plt.subplot(4,3,1)

b = all_pep1.corr()['NHA'][:-1].sort_values()
all_pep1.corr()['Class'][:-1].sort_values().plot(kind='bar',rot=0)
plt.title("Correlation All Peptides")
plt.savefig("Correlation_All_Peptides",  dpi=300)


plt.subplot(4,3,2)
data_train1.corr()['Class'][:-1].sort_values().plot(kind='bar',rot=0)
plt.title("Correlation Data Training / Test I")
plt.savefig("Correlation Data Training_Test I",  dpi=300)

plt.subplot(4,3,3)
data_test2.corr()['Class'][:-1].sort_values().plot(kind='bar',rot=0)
plt.title("Correlation Data Test II")
plt.savefig("Correlation Data Test II",  dpi=300)








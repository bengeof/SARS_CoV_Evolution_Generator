
#Creating the Distance Matrix for clustering 

import os

os.system('python helper.py -i evolved_variants.fasta --metric Levenshtein --reducer UMAP')


#Clustering

import pandas as pd


dfd = pd.read_csv('DM-checkpoint.csv')

from sklearn.cluster import AgglomerativeClustering

model = AgglomerativeClustering(n_clusters=25, linkage='ward').fit(data)
labels = model.fit_predict(data)

print(labels)


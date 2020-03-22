from sklearn import metrics
import crecknn
import math
import numpy as np

# Data load
D = np.load('data_sample/MPEG7-CFD.npy')

# Definition of true labels to compute F-measure
labels_true = [math.floor(i/20) for i in range(len(D))]

# Cluster initialization
cluster = crecknn.Cluster(k=15, ck=4)
cluster.fit(D)

# Cluster evaluation
rand = metrics.adjusted_rand_score(labels_true, cluster.labels_)
nmi_min = metrics.normalized_mutual_info_score(labels_true, cluster.labels_, 'min')
ami = metrics.adjusted_mutual_info_score(labels_true, cluster.labels_, 'min')
ocv = metrics.homogeneity_completeness_v_measure(labels_true, cluster.labels_)
v = metrics.v_measure_score(labels_true, cluster.labels_)
print("adjusted rand score -> ",rand)
print("Normalized mutual info score -> ", nmi_min)
print("Adjusted mutual info score -> ", ami)
print("Homogeneity / Completeness / V-Measure -> ", ocv)
print("V-Measure -> ", v)
print('F-Measure -> ', cluster.calculate_f_measure(cluster.labels_, labels_true))
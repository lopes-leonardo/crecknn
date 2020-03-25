# crecknn (Clustering through Reciprocal kNN graph and connected components)

## Abstract

The huge increase in the amount of multimedia data available and the pressing need for organizing them in different categories, especially in scenarios where there are no labels available, makes data clustering an essential task in different scenarios.
In this work, we present a novel clustering method based on an unsupervised manifold learning algorithm, in which a more effective similarity measure is computed by the manifold learning and used for clustering purposes.
The proposed approach is applied to anomaly detection in videos and used in combination with different background segmentation methods to improve their effectiveness.
An experimental evaluation is conducted on three different image datasets and one video dataset.
The obtained results indicate superior accuracy in most clustering tasks when compared to the baselines.
Results also demonstrate that the clustering step can improve the results of background subtraction approaches in the majority of cases.

## Overview

![Process Workflow](/img/process_image.png)

The **crecknn** works based on three stages:

1. Ranked lists are extracted from the input data, by using the Euclidean distance
2. Those ranked lists are applied to a [manifold learning algorithm](https://www.sciencedirect.com/science/article/abs/pii/S0031320317301978) in order to retrieve a new distance measure and an enhanced set of ranked lists from the data
2. Based on the new enhanced ranked lists retrieved from the stage 1, the algorithm is repeated with a low neighborhood size in order to extract high-reliable initial cluster configuration from the dataset
3. The initial clustering configuration is agglomerated based on the enhanced distance measure, also obtained on stage 1


## Instalation

Clone the repository and enter the project folder:

```
cd crecknn
```

Create and activate a python virtual enviroment:

```
virtualenv venv
source venv/bin/activate
```

Install the project dependencies:

```
pip install -r requirements.txt
```

Compile the source files:

```
make
```

After the described steps, you can run the example code in order to check the instalation:

```
python example.py
```

## Parameters and utilization

The **crecknn** contains three main parameters, one being optional:

Parameter | Decription | Default Value
--------- | ---------- | -------------
k | The size of reciprocal neighborhood exploited by the manifold learning algorithm. | 15
ck | The low size reciprocal neighborhood applied on stage 2, in order to retrieve the first cluster configuration from the dataset | 3
min_elem | The minimum number of elements of each cluster after the agglomeration on stage 3. This parameter defines the number of cluster fusions executed by the algorithm | Being optional, the default value will be equal to k

To utilize the **crecknn**, first import the package on your python code:

```python
import crecknn
```
(Please note that the python code must be on the same folder where the `make` command was executed)

Then, initialize the cluster object:

```python
cluster = crecknn.Cluster(k=15, ck=4, min_elem=10)
```

Finally, insert the input data:

```python
cluster.fit(data)
```

Notice that `data` must be array-like and can have shape(n, m), being interpreted as a feature matrix, or shape(n, n), being interpreted as a distance matrix.

After the clusterization, the labels can be retrieve by:

```python
cluster.labels_
```

## Acknowledgements

The authors are grateful to the São Paulo Research Foundation - FAPESP (#2013/07375-0, #2014/12236-1, #2017/25908-6, #2018/15597-6, #2018/21934-5, #2019/07825-1, and #2019/02205-5), the Brazilian National Council for Scientific and Technological Development - CNPq (#308194/2017-9, #307066/2017-7, and #427968/2018-6), and Petrobras (#2017/00285-6).

## Reference

Lopes, L.; Valem, L.; Pedronette, D.; Guilherme, I.; Papa, J.; Santana, M. and Colombo, D. (2020). **Manifold Learning-based Clustering Approach Applied to Anomaly Detection in Surveillance Videos**.In Proceedings of the 15th International Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications - Volume 5: VISAPP, ISBN 978-989-758-402-2, ISSN 2184-4321, pages 404-412. DOI: 10.5220/0008974604040412

For more information, experimental procedures and achieved results refer to the final [publication](https://www.scitepress.org/PublicationsDetail.aspx?ID=XNW/adjpjh0=&t=1).

If you use our algorithm, please refer to:

```
@conference{visapp20,
author={Leonardo Tadeu Lopes. and Lucas Pascotti Valem. and Daniel Carlos Guimarães Pedronette. and Ivan Rizzo Guilherme. and João Paulo Papa. and Marcos Cleison Silva Santana. and Danilo Colombo.},
title={Manifold Learning-based Clustering Approach Applied to Anomaly Detection in Surveillance Videos},
booktitle={Proceedings of the 15th International Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications - Volume 5: VISAPP,},
year={2020},
pages={404-412},
publisher={SciTePress},
organization={INSTICC},
doi={10.5220/0008974604040412},
isbn={978-989-758-402-2},
}
``` 


Feel free to [contact us](mailto:leonardo.lopes@unesp.br) for any questions or sugestions.
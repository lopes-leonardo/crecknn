# distutils: language = c++

from Data cimport Distance, Item, Database
from libcpp.map cimport map
from libcpp.vector cimport vector
from sklearn import metrics
from sklearn.neighbors import BallTree
import itertools
import math
import numpy as np
import time


cdef class Cluster:
    cdef Database database
    cdef Item* item
    cdef int k, t, ck, initialized, min_elem, classSize
    cdef public vector[int] labels_
    cdef vector[vector[int]] D
    cdef set true_pairs

    def __init__(self, 
                  k=15,
                  ck=3,
                  min_elem=None):
        self.initialized = 0
        self.k = k
        self.ck = ck
        self.min_elem = k if min_elem is None else min_elem

    def fit(self,
            data):
        np_data = np.array(data)
        x, y = np_data.shape
        if x == y:
            self.__init_with_distance_matrix(np_data)
        else:
            self.__init_with_features(np_data)
        self.database.clusterize(1, self.k, self.ck)
        self.database.postProcessing(self.min_elem)
        self.labels_ = self.database.getCCMap()


    def __init_with_features(self, points):
        self.database =  Database()
        dataSize = len(points)
        balltree = BallTree(points, leaf_size=50)
        D = []
        for index, point in enumerate(points):
            rank, rank_map = balltree.query([point], k=dataSize)
            D.append(rank[0])
            item = new Item(index, rank[0], rank_map[0])
            self.database.addItem(item)
        self.database.passDistanceRanking()
        self.initialized = 1
        self.D = D

    def __init_with_distance_matrix(self, points):
        self.database =  Database()
        for index, ranking in enumerate(points):
            item = new Item(index, ranking)
            self.database.addItem(item)
        self.database.passDistanceRanking()
        self.initialized = 1
        self.D = points

    def __make_pairs(self, labels):
        elements = {}
        for index, cluster in enumerate(labels):
            if cluster in elements:
                elements[cluster].append(index)
            else:
                elements[cluster] = [index]

        pairs = set()
        for lista in elements.values():
            for par in itertools.combinations(lista,2):
                pairs.add(par)

        return pairs

    def calculate_f_measure(self, predict_labels, labels_true = None):
        if(self.true_pairs is None and labels_true is not None):
            self.true_pairs = self.__make_pairs(labels_true)
        elif(self.true_pairs is None and not self.labels_true.empty()):
            self.true_pairs = self.__make_pairs(self.labels_true)
        elif(self.true_pairs is None):
            print('Missing true labels for F_Measure calculation')
            return None

        predict_pairs = self.__make_pairs(predict_labels)
        a = predict_pairs.intersection(self.true_pairs)
        b = predict_pairs - self.true_pairs
        c = self.true_pairs - predict_pairs

        return float(2*len(a)) / float(2*len(a) + len(b) + len(c))
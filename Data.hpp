/* <Data.hpp>
 *
 *
 * This file contains main data structures for clusterization.
 *
 *
 * @author: Leonardo Tadeu Lopes
 */

#ifndef DATA_H
#define DATA_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <list> 

namespace crecknn {
    class Distance {
        public:
            Distance();
            Distance(int, double);
            int getItem();
            void setItem(int);
            double getValue();
            void setValue(double);
            bool operator <(const Distance& r);
            bool operator ==(const Distance& r);
            double value;
            int item;
        private:
    };

	class Item {
		public:
			Item();
			Item(int, std::vector<double>);
			Item(int, std::vector<double>, std::vector<int>);
			void setRanking(std::vector<double>);
			void setRanking(std::vector<Distance>);
			void setOrderedRanking(std::vector<double>*, std::vector<int>*);
			void setRankingAt(int, Distance);
			int getId();
			std::vector<Distance> getRanking();
			Distance getRankingAt(int);
			void orderRanking(int);
			int getRankingSize();
			int findRanking(int, int);
		private:
			int id;
			std::vector<Distance> ranking;
	};

    class Edge {
        public:
            Edge(int, int);
            int node1, node2;
        private:
    };

    class Node {
        public:
            Node(int, int, int);
            int parentId, imageId, level;
            std::vector<Edge> edges;
            bool hasEdge(int);
            void addEdge(Edge);
            void clearEdges();
        private:
    };

    class ConnectedComponent {
        public:
            ConnectedComponent(int);
            std::vector<Node*> nodes;
            int representant;
            std::vector<int> ids;
            bool isMember(int);
        private:
    };

    std::vector<Node> initializeNodes(int);
    std::vector<std::unordered_set<int>> initializeCCs(int, std::vector<int>*);
    bool findTkNeighbor(std::vector<std::vector<Distance>>*,int,int,int);
    bool isReciprocalNeighbor(std::vector<std::vector<Distance>>*,int,int,int);
    std::vector<std::vector<Distance>> getZeroMatrix(std::vector<std::vector<Distance>>, bool=false, bool=false);
    void printMatrix(std::vector<std::vector<Distance>>);
    void printNodes(std::vector<Node>*, int);
    std::vector<ConnectedComponent> getConnectedComponents(std::vector<Node>*);
    void incrementSimilarity(std::vector<std::vector<Distance>>*,int,int,int);
    void incrementByEdges(std::vector<std::vector<Distance>>*, std::vector<Node>*, int);
    void incrementByConnectedComponents(std::vector<std::vector<Distance>>*, std::vector<std::unordered_set<int>>, int);
    void sortMatrixAsc(std::vector<std::vector<Distance>>*, int);
    void sortMatrixDesc(std::vector<std::vector<Distance>>*, int);
    std::vector<std::vector<Distance>> convertToDistanceMatrix(std::vector<std::vector<Distance>>*);
    void printCCs(std::vector<std::unordered_set<int>>);
    void printCCMap(std::vector<int>);
    void createEdgesForTk(std::vector<Node>*,std::vector<std::vector<Distance>>*,int);
    int findSet(int, std::vector<int>*);
    void updateCCs(std::vector<Node>*,std::vector<std::unordered_set<int>>*,std::vector<int>*);
    void unionSet(std::vector<std::unordered_set<int>>* ccs, std::vector<int>*,int,int);
    void freeMatrix(int**, int);
    void initializePositionMap(int);
    void initializeRankedList(int);
    int findMatrixRanking(int, int, int);

    class Database {
        public:
            Database();
            void addItem(Item*);
            std::vector<Item*> getItems();
            Item* getItem(int);
            void setItem(int, Item*);
            void passDistanceRanking();
            void normalizeRanking();
            void clusterize(int, int, int);
            void runMethod(int, int, bool, bool=true, bool=true, bool=true);
            std::vector<std::vector<Distance>> getDistanceMatrix();
            int getDatabaseSize();
            bool isClusterized();
            void makeOutput();
            void joinClusters(int);
            std::vector<int> getCCMap();
            std::vector<int> getCC();
            void saveDistanceRanking();
            void loadDistanceRanking();
            void saveNormalizedRanking();
            void loadNormalizedRanking();
            void saveOriginalDistanceRanking();
            void loadOriginalDistanceRanking();
            void uniteClustersByElemNumber(int, bool=false);
            void postProcessing(int);
            float calculateClustersDistance(int, int);
            void fusionClusters(int, int, bool=false);
        private:     
            std::vector<Item*> items;
            std::vector<std::vector<Distance>> distanceMatrix, similarityMatrix, savedDistanceMatrix, originalDistanceMatrix, normalizedDistanceMatrix;
            int k, t, ck, l, clusterized;
            std::vector<int> ccMap, outputCCMap, outputCC;
            std::vector<std::unordered_set<int>> ccs, outputCCItens;
    };
}
#endif

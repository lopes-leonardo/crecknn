/* <Data.cpp>
 *
 * Objects method implementation file
 *
 * @Author: Leonardo Tadeu Lopes
 *
 *
 *****************************************************************************************************************
 *
 *
 * This file contains the methods implementation for rec-knn-cluster object classes
 */

#include "Data.hpp"

namespace crecknn {
    // GLOBAL VARIABLES ############################
    int** positionMap = nullptr;
    int** rankedList = nullptr;

    // DISTANCE ####################################################
    Distance::Distance() {}

    Distance::Distance(int item, double value) {
        this->item = item;
        this->value = value;
    }

    int Distance::getItem() {
        return this->item;
    }

    void Distance::setItem(int item) {
        this->item = item;
    }


    double Distance::getValue() {
        return this->value;
    }

    void Distance::setValue(double value) {
        this->value = value;
    }

    bool Distance::operator <(const Distance& r) {
        return value < r.value;
    }

    bool Distance::operator ==(const Distance& r) {
        return item == r.item;
    }



    // ITEM ########################################################
    Item::Item() {}

    Item::Item(int id, std::vector<double> list) {
        this->id = id;
        this->setRanking(list);
        this->orderRanking(list.size());
    }

    Item::Item(int id, std::vector<double> ranking, std::vector<int> rankingMap) {
        this->id = id;
        this->setOrderedRanking(&ranking, &rankingMap);
    }

    void Item::setOrderedRanking(std::vector<double>* ranking, std::vector<int>* map) {
        int dataSize = ranking->size();
        for(int x = 0; x < dataSize; x++) {
            Distance distance (map->at(x), ranking->at(x));
            this->ranking.push_back(distance);
        }
    }

    void Item::setRanking(std::vector<double> rawRanking) {
        int rawSize = rawRanking.size();
        for(int x = 0; x < rawSize; x++) {
            Distance distance (x, rawRanking.at(x));
            this->ranking.push_back(distance);
        }
    }

    void Item::setRanking(std::vector<Distance> ranking) {
        this->ranking = ranking;
    }

    Distance Item::getRankingAt(int position) {
        return this->ranking.at(position);
    }

    std::vector<Distance> Item::getRanking() {
        return this->ranking;
    }

    int Item::getId() {
        return this->id;
    }

    void Item::orderRanking(int size) {
        std::stable_sort(ranking.begin(), ranking.end(), [ ]( const Distance lhs, const Distance rhs )
        {
            return lhs.value < rhs.value;
        });
    }

    void Item::setRankingAt(int position, Distance newDistance) {
        this->ranking.at(position) = newDistance;
    }

    int Item::findRanking(int id, int maxLevel) {
        // If id equals Item->id, there's no need to initiate a search
        if(id == this->getId()) {
            return 1;
        }

        int x;
        for(x=1; x < maxLevel; x++)
        {
            Distance rank = ranking.at(x);
            if(rank.getItem() == id)
            {
                break;
            }
        }

        // Add 1 on x to return between 1 and maxLevel
        return x;
    }

    int Item::getRankingSize() {
        return ranking.size();
    }

    // EDGE ###################################################
    Edge::Edge(int a, int b) {
        this->node1 = a;
        this->node2 = b;
    }

    // NODE ###################################################
    Node::Node(int a, int b, int c) {
        this->imageId = a;
        this->parentId = b;
        this->level = c;
    }

    bool Node::hasEdge(int id)
    {
        int edgesSize = edges.size();
        if(edgesSize == 0)
        {
            return false;
        }
        for(int x = 0; x < edgesSize; x++)
        {
            Edge edge = edges.at(x);
            if(edge.node2 == id)
            {
                return true;
            }
        }
        return false;
    }

    void Node::addEdge(Edge newEdge) {
        this->edges.push_back(newEdge);
    }

    void Node::clearEdges() {
        this->edges.clear();
    }

    // CONNECTED COMPONENT ###################################################
    ConnectedComponent::ConnectedComponent(int id)
    {
        this->representant = id;
    }

    // SUPPORT FUNCTIONS
    std::vector<Node> initializeNodes(int size) {
        std::vector<Node> nodes;
        for(int x = 0; x < size; x++)
        {
            Node node (x, x, 1);
            nodes.push_back(node);
        }
        return nodes;
    }
    std::vector<std::unordered_set<int>> initializeCCs(int size, std::vector<int>* ccMap) {
        std::vector<std::unordered_set<int>> ccs;
        for(int i = 0; i < size; i++) {
            std::unordered_set<int> set;
            set.insert(i);
            ccs.push_back(set);
            ccMap->push_back(i);
        }
        return ccs;
    }
    bool findTkNeighbor(std::vector<std::vector<Distance>> *similarityMatrix, int source, int target, int tk) {
        std::vector<Distance> imageRanking = similarityMatrix->at(source);
        for(int x = 0; x < tk; x++) {
            Distance distance = imageRanking.at(x);
            if(distance.getItem() == target) {
                return true;
            }
        }
        return false;
    }
    bool isReciprocalNeighbor(std::vector<std::vector<Distance>> *similarityMatrix, int image1, int image2, int tk) {
        return findTkNeighbor(similarityMatrix, image1, image2, tk) && findTkNeighbor(similarityMatrix, image2, image1, tk);
    }
    std::vector<std::vector<Distance>> getZeroMatrix(std::vector<std::vector<Distance>> distanceMatrix, bool makePositionMap, bool makeRankedList) {
        std::vector<std::vector<Distance>> zeroMatrix;
        int size = distanceMatrix.size();
        for(int x = 0; x < size; x++) {
            std::vector<Distance> distances =  distanceMatrix.at(x);
            std::vector<Distance> zeroDistances;
            int size2 = distances.size();
            for(int y = 0; y < size2; y++) {
                Distance distance = distances.at(y);
                distance.setValue(0);
                zeroDistances.push_back(distance);
                if(makePositionMap) {
                    positionMap[x][distance.getItem()] = y;
                }
                if(makeRankedList) {
                    rankedList[x][y] = distance.getItem();
                }
            }
            zeroMatrix.push_back(zeroDistances);
        }
        return zeroMatrix;
    }
    void incrementSimilarity(std::vector<std::vector<Distance>> *similarityMatrix, int image1, int image2, int increment)
    {
        int position = positionMap[image1][image2];
        std::vector<Distance>* imageRanking = &similarityMatrix->at(image1);
        Distance* distance = &imageRanking->at(position);
        distance->value += increment;
    }
    void incrementByEdges(std::vector<std::vector<Distance>> *similarityMatrix, std::vector<Node> *nodes, int increment) {
        int nodesSize = nodes->size();
        for(int imageId = 0; imageId < nodesSize; imageId++)
        {
            Node imageNode = nodes->at(imageId);
            int edgesSize = imageNode.edges.size();
            for(int edge1Id = 0; edge1Id < edgesSize; edge1Id++)
            {
                Edge edge1 = imageNode.edges.at(edge1Id);
                int item1 = edge1.node2;
                for(int edge2Id = 0; edge2Id < edgesSize; edge2Id++) {
                    Edge edge2 = imageNode.edges.at(edge2Id);
                    int item2 = edge2.node2;
                    incrementSimilarity(similarityMatrix, item1, item2, increment);
                }
            }
        }
    }
    void incrementByConnectedComponents(std::vector<std::vector<Distance>> *similarityMatrix, std::vector<std::unordered_set<int>> ccs, int increment) {
        // Iterate over every item
        int size = similarityMatrix->size();
        std::unordered_set<int>::iterator it1, it2;
        for(int itemId = 0; itemId < size; itemId++) {
            std::unordered_set<int> cc = ccs.at(itemId);
            if(!cc.empty()) {
                // Iterate twice over cc items to increment all possibilities
                for(it1 = cc.begin(); it1 != cc.end(); ++it1) {
                    for(it2 = cc.begin(); it2 != cc.end(); ++it2) {
                        incrementSimilarity(similarityMatrix, *it1, *it2, increment);
                    }
                }
            }

        }
    }
    void sortMatrixAsc(std::vector<std::vector<Distance>> *similarityMatrix, int tk) {
        int matrixSize = similarityMatrix->size();
        for(int imageId = 0; imageId < matrixSize; imageId++) {
            std::vector<Distance> imageRanking = similarityMatrix->at(imageId);
            std::stable_sort(imageRanking.begin(), imageRanking.begin() + tk, [ ]( const Distance lhs, const Distance rhs )
            {
                return lhs.value < rhs.value;
            });
            similarityMatrix->at(imageId) =  imageRanking;
        }
    }
    void sortMatrixDesc(std::vector<std::vector<Distance>> *similarityMatrix, int tk) {
        int matrixSize = similarityMatrix->size();
        for(int imageId = 0; imageId < matrixSize; imageId++) {
            std::vector<Distance> imageRanking = similarityMatrix->at(imageId);
            std::stable_sort(imageRanking.begin(), imageRanking.begin() + tk, [ ]( const Distance lhs, const Distance rhs )
            {
                return lhs.value > rhs.value;
            });
            similarityMatrix->at(imageId) =  imageRanking;
        }
    }
    std::vector<std::vector<Distance>> convertToDistanceMatrix(std::vector<std::vector<Distance>>* similarityMatrix) {
        std::vector<std::vector<Distance>> distanceMatrix;
        int size = similarityMatrix->size();
        for(int i = 0; i < size; i++) {
            std::vector<Distance> itemRanking = similarityMatrix->at(i);
            std::vector<Distance> newRanking;
            for(int pos = 0; pos < size; pos++) {
                Distance distance = itemRanking.at(pos);
                distance.value = 1 / (1 + distance.value);
                // std::cout << "Distance from item " << i << " to item " << distance.getItem() << ": " << distance.value << "\n";
                newRanking.push_back(distance);
            }
            distanceMatrix.push_back(newRanking);
        }
        sortMatrixAsc(&distanceMatrix, size);
        return distanceMatrix;
    }
    void createEdgesForTk(std::vector<Node>* nodes, std::vector<std::vector<Distance>>* similarityMatrix, int tk){
        // Iterate item nodes to generate edges
        int size = nodes->size();
        for(int itemId = 0; itemId < size; itemId++) {
            Node* itemNode = &nodes->at(itemId);
            std::vector<Distance> itemRanking = similarityMatrix->at(itemId);

            // Get item distance at "tk-1" (tk runs from 1 to K) rank position and knnItem id
            Distance distance = itemRanking.at(tk-1);
            int knnItemId = distance.getItem();

            // Is item already connected to knnItem?
            // If true there is no need to check it again
            if(!itemNode->hasEdge(knnItemId)) {
                // Obtain knnItem individual Node
                if(isReciprocalNeighbor(similarityMatrix, itemId, knnItemId, tk)) {
                    // Create edge from item to knnItem
                    Edge newEdge (itemId, knnItemId);
                    itemNode->addEdge(newEdge);

                    // If item and knnItem are differents, create the reciprocal edge
                    if(itemId != knnItemId) {
                        Node* knnItemNode = &nodes->at(knnItemId);
                        Edge reciprocalEdge (knnItemId, itemId);
                        knnItemNode->addEdge(reciprocalEdge);
                    }
                }
            }
        }
    }
    int findSet(int imageId, std::vector<int>* ccMap) {
        int parentId = ccMap->at(imageId);
        if(parentId == imageId) {
            return imageId;
        }
        else
        {
            int newParent = findSet(parentId, ccMap);
            ccMap->at(imageId) = newParent;
            return newParent;
        }
    }
    void updateCCs(std::vector<Node>* nodes, std::vector<std::unordered_set<int>>* ccs, std::vector<int>* ccMap){
        // Iterate over nodes
        int size = nodes->size();
        for(int itemId = 0; itemId < size; itemId++) {
            Node* itemNode = &nodes->at(itemId);

            // Iterate over node edges
            int edgesSize =  itemNode->edges.size();
            for(int y = 0; y < edgesSize; y++) {
                Edge edge = itemNode->edges.at(y);

                // With edge, compares if items are in the same CC
                if(findSet(edge.node1, ccMap) != findSet(edge.node2, ccMap)) {
                    unionSet(ccs, ccMap, edge.node1, edge.node2);
                }
            }
        }
    }
    void unionSet(std::vector<std::unordered_set<int>>* ccs, std::vector<int>* ccMap, int item1, int item2) {
        int leader1 = ccMap->at(item1);
        int leader2 = ccMap->at(item2);
        std::unordered_set<int>* set1 =  &ccs->at(leader1);
        std::unordered_set<int>* set2 =  &ccs->at(leader2);
        std::unordered_set<int>::iterator it;
        for(it = set2->begin(); it != set2->end(); ++it) {
            int item =  *it;
            set1->insert(item);
            ccMap->at(item) = leader1;
        }
        set2->clear();
    }
    void freeMatrix(int** pointer, int size) {
        if(pointer != nullptr) {
            for(int i = 0; i < size; i++) {
                int* item = pointer[i];
                free(item);
            }
            free(pointer);
            pointer = nullptr;
        }
    }
    void initializePositionMap(int size) {
        if(positionMap != nullptr) {
            freeMatrix(positionMap, size);
        }
        positionMap = (int**)malloc(sizeof(int*) * size);
        for(int i = 0; i < size; i++) {
            positionMap[i] = (int*)malloc(sizeof(int) * size);
        }
    }
    void initializeRankedList(int size) {
        freeMatrix(rankedList, size);
        rankedList = (int**)malloc(sizeof(int*) * size);
        for(int i = 0; i < size; i++) {
            rankedList[i] = (int*)malloc(sizeof(int) * size);
        }
    }
    int findMatrixRanking(std::vector<std::vector<Distance>>* distanceMatrix, int sourceId, int id, int maxLevel) {
        // If id equals Item->id, there's no need to initiate a search
        if(id == sourceId) {
            return 1;
        }

        std::vector<Distance>* ranking = &distanceMatrix->at(sourceId);
        int x;
        for(x=1; x < maxLevel; x++)
        {
            Distance* rank = &ranking->at(x);
            if(rank->getItem() == id)
            {
                break;
            }
        }

        // Add 1 on x to return between 1 and maxLevel
        return x;
    }

    // DATABASE ###################################################
    Database::Database() {
        this->clusterized = 0;
    }

    void Database::addItem(Item* item) {
        this->items.push_back(item);
    }

    bool Database::isClusterized() {
        return (this->clusterized == 1);
    }

    std::vector<Item*> Database::getItems() {
        return this->items;
    }

    void Database::setItem(int position, Item* item) {
        items.at(position) = item;
    }

    Item* Database::getItem(int position) {
        return items.at(position);
    }

    void Database::passDistanceRanking() {
        this->distanceMatrix = this->getDistanceMatrix();
    }

    void Database::saveDistanceRanking() {
        this->savedDistanceMatrix = this->distanceMatrix;
    }

    void Database::loadDistanceRanking() {
        this->distanceMatrix = this->savedDistanceMatrix;
    }

    void Database::saveNormalizedRanking() {
        this->normalizedDistanceMatrix = this->distanceMatrix;
    }

    void Database::loadNormalizedRanking() {
        this->distanceMatrix = this->normalizedDistanceMatrix;
    }

    void Database::saveOriginalDistanceRanking() {
        this->originalDistanceMatrix = this->distanceMatrix;
    }

    void Database::loadOriginalDistanceRanking() {
        this->distanceMatrix = this->originalDistanceMatrix;
    }

    // Normalize on distanceMatrix
    void Database::normalizeRanking() {
        int itemsSize = this->getDatabaseSize(), normalizedValue, itemId, item2Id, position, normalizedY;
        int maxSize = itemsSize * itemsSize;
        for(int x = 0; x < itemsSize; x++) {
            itemId = x;
            std::vector<Distance>* itemRanking = &this->distanceMatrix.at(x);
            // Calculate the L fist positions normalized value
            for(int y = 0; y < itemsSize; y++) {
                Distance* ranking = &itemRanking->at(y);
                if(y < this->l) {
                    item2Id = ranking->getItem();
                    position = findMatrixRanking(&this->distanceMatrix, item2Id, itemId, itemsSize);
                    normalizedY = y + 1;
                    normalizedValue = normalizedY + position + std::max(normalizedY, position);
                } else {
                    normalizedValue = maxSize;
                }

                // Updates item ranking on position y
                ranking->value = normalizedValue;
            }
        }

        // Sort the L first itens by minor normalized value
        for(int x = 0; x < itemsSize; x++) {
            std::vector<Distance>* itemRanking = &this->distanceMatrix.at(x);
            std::stable_sort(itemRanking->begin(), itemRanking->end(), [ ]( const Distance lhs, const Distance rhs )
            {
                return lhs.value < rhs.value;
            });
        }
    }

    std::vector<std::vector<Distance>> Database::getDistanceMatrix() {
        std::vector<std::vector<Distance>> distanceMatrix;
        int itemsSize = this->getDatabaseSize();
        for(int x = 0; x < itemsSize; x++) {
            Item* item = items.at(x);
            distanceMatrix.push_back(item->getRanking());
        }
        return distanceMatrix;
    }

    int Database::getDatabaseSize() {
        return this->items.size();
    }

    void Database::clusterize(int t, int k, int ck) {
        this->k = k;
        this->l = 4 * k;
        this->runMethod(t, k, false, true, true);
        std::vector<std::vector<Distance>> distanceMatrix = this->distanceMatrix;
        this->runMethod(1, ck, true, false, false, true);
        this->distanceMatrix = distanceMatrix;
    }

    void Database::runMethod(int t, int k, bool saveResult, bool saveDistance, bool applyIncrement, bool normalize) {
        int databaseSize = this->getDatabaseSize();
        std::vector<Node> nodes;
        std::vector<int> ccMap;
        std::vector<std::unordered_set<int>> ccs;
        std::vector<std::vector<Distance>> similarityMatrix;
        for(int tn = 1; tn <= t; tn++) {
            if(normalize) {
                this->normalizeRanking();
            }
            initializePositionMap(databaseSize);
            similarityMatrix = getZeroMatrix(this->distanceMatrix, true, false);
            nodes = initializeNodes(databaseSize);
            ccs = initializeCCs(databaseSize, &ccMap);
            for(int tk = 1; tk <= k; tk ++) {
                createEdgesForTk(&nodes, &similarityMatrix, tk);
                updateCCs(&nodes, &ccs, &ccMap);

                // With the filled graph, increment the similarityMatrix by Edges and ConnectedComponents
                if(applyIncrement) {
                    int increment = k - tk + 1;
                    incrementByEdges(&similarityMatrix, &nodes, increment);
                    incrementByConnectedComponents(&similarityMatrix, ccs, increment);
                }
            }
            if(saveDistance) {
                this->distanceMatrix = convertToDistanceMatrix(&similarityMatrix);
                sortMatrixDesc(&similarityMatrix, databaseSize);
                this->similarityMatrix = similarityMatrix;
            }
        }
        if(saveResult) {
            this->ccMap = ccMap;
            this->ccs = ccs;
            this->makeOutput();
            this->clusterized = 1;
        }
    }

    float Database::calculateClustersDistance(int a, int b) {
        float totalDistances = 0.0;
        std::unordered_set<int> cca = this->outputCCItens[a], ccb = this->outputCCItens[b];
        float sizea = cca.size();

        for(int i : cca) {
            // std::cout << "Calculating for cc " << a << " item " << i << "\n";
            float itemDistance = 0.0;
            for(int j : ccb) {
                float value = this->distanceMatrix[i][positionMap[i][j]].getValue();
                // std::cout << "CC " << b << " item " << j << " with distance " << value << "\n";
                itemDistance += value;
            }
            totalDistances += itemDistance / (float)ccb.size();
        }

        return totalDistances / sizea;
    }

    void Database::fusionClusters(int a, int b, bool updateSilhouette) {
        std::unordered_set<int>* cca = &this->outputCCItens.at(a);
        std::unordered_set<int>* ccb = &this->outputCCItens.at(b);
        std::unordered_set<int>::iterator it;
        for(it = cca->begin(); it != cca->end(); ++it) {
            int item =  *it;
            ccb->insert(item);
            this->outputCCMap.at(item) = b;
        }
        cca->clear();
        this->outputCC[b] += this->outputCC[a];
        this->outputCC[a] = 0;

        // Aqui virá a atualização das distâncias e, consequentemente, do Silhouette Score dos items envolvidos 
    }

    void Database::uniteClustersByElemNumber(int elemNumber, bool updateSilhouette) {
        int numberOfClusters = this->outputCC.size();
        std::list<Distance> targetClusters;
        // Iterate over clusters, spliting into target and remaining clusters
        // based on i elements
        for(int j = 0; j < numberOfClusters; j++) {
            if(this->outputCC[j] == elemNumber) {
                Distance d (j, 0);
                targetClusters.push_back(d);
            }
        }
        for(const Distance d : targetClusters) {
            std::vector<Distance> distances;
            for(int y = 0; y < numberOfClusters; y++) {
                if(y != d.item && this->outputCC[y] != 0) {
                    float value = this->calculateClustersDistance(d.item, y);
                    Distance dist (y, value);
                    distances.push_back(dist);
                }
            }
            std::stable_sort(distances.begin(), distances.end(), [ ]( const Distance lhs, const Distance rhs )
            {
                return lhs.value < rhs.value;
            });

            int choosenCluster = distances[0].getItem();
            this->fusionClusters(d.item, choosenCluster);
        }
    }

    void Database::postProcessing(int minElem) {
        // Initialize variables
        if(positionMap == nullptr) {
            initializePositionMap(this->getDatabaseSize());
        }
        std::vector<std::vector<Distance>> similarityMatrix = getZeroMatrix(this->distanceMatrix, true, false);
        
        // Arbitrary fusion of unitary cluster to initiate the heuristic evaluation of the silhouette score
        this->uniteClustersByElemNumber(1);

        // Initialize i = 2, to iterate over the minors clusters, fusioning them until the Silhouette score gets worse
        int i = 2;
        while(i < minElem) {
            this->uniteClustersByElemNumber(i, false);
            i++;
        }

        this->ccs = this->outputCCItens;
        this->ccMap = this->outputCCMap;
        this->makeOutput();
    }

    void Database::makeOutput() {
        std::map<int,int> ccConvertMap;
        std::vector<int> outputCC, outputCCMap;
        std::vector<std::unordered_set<int>> outputCCItens;
        int counter = 0;
        int size = this->ccs.size();
        for(int i = 0; i < size; i++) {
            std::unordered_set<int> cc = this->ccs.at(i);
            if(!cc.empty()) {
                ccConvertMap.insert({i, counter});
                outputCC.push_back(cc.size());
                outputCCItens.push_back(cc);
                counter++;
            }
        }
        for(int j = 0; j < this->getDatabaseSize(); j++) {
            int ccNumber = this->ccMap.at(j);
            outputCCMap.push_back(ccConvertMap.at(ccNumber));
        }

        this->outputCCMap = outputCCMap;
        this->outputCC = outputCC;
        this->outputCCItens = outputCCItens;
    }

    std::vector<int> Database::getCCMap() {
        return this->outputCCMap;
    }

    std::vector<int> Database::getCC() {
        return this->outputCC;
    }
}

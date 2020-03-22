# distutils: sources = Rectangle.cpp

cdef extern from "Data.cpp":
    pass

# Use vector, necessary?
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as bool_t

# Declare Image with cdef
cdef extern from "Data.hpp" namespace "crecknn":
    cdef cppclass Distance:
        Distance() except +
        Distance(int, double) except +
        int getItem()
        setItem(int)
        double getValue()
        setValue(double value) except +
    cdef cppclass Item:
        Item() except +
        Item(int, vector[double]) except +
        Item(int, vector[double], vector[int]) except +
        setRanking(vector[Distance]) except +
        vector[Distance] getRanking()
        int getId()
        void orderRanking()
        float getAuthorityScore()
    cdef cppclass Database:
        Database() except +
        void addItem(Item*) except +
        vector[Item*] getItems()
        setItem(int, Item*) except +
        Item* getItem(int)
        void generateDistanceRanking()
        void passDistanceRanking() except +
        void normalizeRanking()
        void clusterize(int, int, int)
        bool_t isClusterized()
        vector[int] getCCMap()
        vector[int] getCC()
        void runMethod(int, int, bool, bool, bool, bool)
        void saveDistanceRanking();
        void loadDistanceRanking();
        void saveNormalizedRanking();
        void loadNormalizedRanking();
        void saveOriginalDistanceRanking()
        void loadOriginalDistanceRanking()
        void postProcessing(int);
        void uniteClustersByElemNumber(int, bool);
        void makeOutput();

#pragma once
#include <vector>

class ShapePredictorTrainer{
public:
    ShapePredictorTrainer(){
        cascadeDepth = 10;
        treeDepth = 4;
        numTreesPerCascadeLevel = 500;
        num = 0.1;
        oversamplingAmount = 20;
        featurePoolSize = 400;
        lambda = 0.1;
    }
    size_t cascadeDepth;
    size_t treeDepth;
    size_t numTreesPerCascadeLevel;
    double nu;
    size_t oversamplingAmount;
    size_t featurePoolSize;
    double lambda;
};

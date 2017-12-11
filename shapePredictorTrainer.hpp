#pragma once
#include <vector>

class ShapePredictorTrainer {
public:
    ShapePredictorTrainer() {
        cascadeDepth = 10;
        treeDepth = 4;
        numTreesPerCascadeLevel = 500;
        num = 0.1;
        oversamplingAmount = 20;
        featurePoolSize = 400;
        lambda = 0.1;
    }
private:
    size_t cascadeDepth;
    size_t treeDepth;
    size_t numTreesPerCascadeLevel;
    double nu;
    size_t oversamplingAmount;
    size_t featurePoolSize;
    double lambda;
};

class trainingSample {
public:
    size_t imageIdx;
    std::vector<size_t> rect;
    std::vector<float> targetShape;
    std::vector<float> present;
    std::vector<float> currentShape;
    std::vector<float> diffShape;
    std::vector<unsigned char> featurePixelValues;
};

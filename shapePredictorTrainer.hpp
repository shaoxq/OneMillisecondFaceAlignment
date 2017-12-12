#pragma once
#include <vector>
#include "rectangle.hpp"
#include "fullObjectDetection.hpp"

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
    void train(std::vector<Image>& images, std::vector<fullObjectDetection>& objects);
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
    std::vector<Rectangle> rect;
    std::vector<float> targetShape;
    std::vector<float> present;
    std::vector<float> currentShape;
    std::vector<float> diffShape;
    std::vector<unsigned char> featurePixelValues;
};

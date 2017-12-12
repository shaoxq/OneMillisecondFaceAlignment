#pragma once
#include <vector>
#include <cassert>
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
    void train(std::vector<Image>& images, std::vector<fullObjectDetection>& objects) {
        assert(images.size() == objects.size());

        // check whether all the objects have same number of parts
        size_t numParts = 0;
        for (size_t i = 0; i < objects.size(); i++) {
            if (numParts == 0){
                numParts = objects.parts.size();
                assert(numParts != 0);
            } else {
                assert(objects.parts.size() == numParts);
            }
        }

        

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
    std::vector<Rectangle> rect;
    std::vector<float> targetShape;
    std::vector<float> present;
    std::vector<float> currentShape;
    std::vector<float> diffShape;
    std::vector<unsigned char> featurePixelValues;
};

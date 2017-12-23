#pragma once
#include <vector>
#include <cassert>
#include "rectangle.hpp"
#include "fullObjectDetection.hpp"
#include "shapePredictor.hpp"

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
    void train(std::vector<Image>& images, std::vector<fullObjectDetection>& objects) const {
        assert(images.size() == objects.size());

        // check whether all the objects have same number of parts
        size_t numParts = objects[0].parts.size();
        assert(numParts != 0);
        std::vector<bool> partPresent(objects[0].parts.size(),false);
        for (size_t i = 0; i < objects.size(); i++) {
            assert(objects[i].parts.size() == numParts);
            for (size_t j = 0; j < objects[i].parts.size(); j++) {
                if (objects[i].isPresent[j] == true) {
                    partPresent[j] = true;
                }
            }
        }

        assert(numParts != 0);
        for (size_t i = 0; i < partPresent.size(); i++) {
            assert(partPresent[i] == true);
        }

        std::vector<TrainingSample> samples;
        std::vector<Point2d<float>> initialShape = populateTrainingSampleShapes(objects, samples);
        std::vector<std::vector<Point2d<float>>> pixelCoordinates = randomlySamplePixelCoordinates(initialShape);

        std::vector<std::vector<regressionTree>> forests(cascadeDepth);
        for (size_t cascade = 0; cascade < cascadeDepth; cascade++) {
            std::vector<size_t> anchorIdx;
            std::vector<Point2d<float>> deltas;
            createShapeRelativeEncoding(initialShape,pixelCoordinates[cascadeDepth],anchorIdx,deltas);

            for (size_t i = 0; i < samples.size(); i++) {
                extractFeaturePixelValues(images[samples[i].imageIdx],samples[i].rect,samples[i].currentShape,initialShape,anchorIdx,deltas,samples[i].featurePixelValues);
            }

            for (size_t i = 0; i < numTreesPerCascadeLevel; i++) {
                forests[cascade].push_back(makeRegressionTree(samples, pixelCoordinates[cascade]);
            }
        }
    }

    std::vector<Point2d<float>> populateTrainingSampleShapes(const std::vector<fullObjectDetection>& objects, std::vector<TrainingSample>& samples) const {
        std::vector<Point2d<float>> meanShape(objects[0].parts.size());
        std::vector<size_t> count(objects[0].parts.size(), 0);

        // first fill out the target shapes
        for (size_t i = 0; i < objects.size(); i++) {
            TrainingSample sample;
            sample.imageIdx = i;
            sample.rect = objects[i].rect;

            std::vector<Point2d<double>> sampleRect(4);
            sampleRect[0].x = sample.rect.leftTop.x;
            sampleRect[0].y = sample.rect.leftTop.y;
            sampleRect[1].x = sample.rect.rightBottom.x;
            sampleRect[1].y = sample.rect.leftTop.y;
            sampleRect[2].x = sample.rect.leftTop.x;
            sampleRect[2].y = sample.rect.rightBottom.y;
            sampleRect[3].x = sample.rect.rightBottom.x;
            sampleRect[3].y = sample.rect.rightBottom.y;

            std::vector<Point2d<double>> unitRect(4);
            unitRect[0].x = 0.0;
            unitRect[0].y = 0.0;
            unitRect[1].x = 1.0;
            unitRect[1].y = 0.0;
            unitRect[2].x = 0.0;
            unitRect[2].y = 1.0;
            unitRect[3].x = 1.0;
            unitRect[3].y = 1.0;

            // align object's rectangle to unit one
            PointTransformAffine tform = findAffinTransform(sampleRect, unitRect);
            for (size_t j = 0; j < objects[i].parts.size(); j++) {
                if (objects[i].isPresent[j] == true) {
                    sample.targetShape[j].x = tform.s * (tform.m[0][0] * objects[i].parts[j].x + tform[0][1] * objects[i].parts[j].y) + tform.b[0];
                    sample.targetShape[j].y = tform.s * (tform.m[1][0] * objects[i].parts[j].x + tform[1][1] * objects[i].parts[j].y) + tformb[1];
                    sample.isPresent[j] = true;
                } else {
                    sample.targetShape[j].x = 0;
                    sample.targetShape[j].y = 0;
                    sample.isPresent[j] = false;
                }
            }

            for (size_t j = 0; j < oversamplingAmount; j++) {
                samples.push_back(sample);
            }

            for (size_t j = 0; j < objects[i].parts.size(); j++) {
                if (sample.isPresent == true) {
                    meanShape[j].x += sample.targetShape[j].x;
                    meanShape[j].y += sample.targetShape[j].y;
                    count[j]++;
                }
            }
        }

        for (size_t j = 0; j < objects[i].parts.size(); j++) {
            meanShape[j].x /= count[j];
            meanShape[j].y /= count[j];
        }

        // new go pick random initial shapes
        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution<double> norm( 0.0, 1.0 );
        for (size_t i = 0; i < samples.size(); i++) {
            if (overSamplingAmount == 0) {
                samples.currentShape[i] = meanShape;
            } else {
                // pick a few samples at random and randomly average them together to make the initial shape
                std::vector<double> hits(samples[i].targetShape.size(),0);
                bool isHitAll = false;
                int itr = 0;
                while (isHitAll == false || itr < 2) {
                    isHitAll = true;
                    size_t rndIdx = mt() % samples.size();
                    double alpha = norm(mt) + 0.1;
                    for (size_t k = 0; k < samples[i].currentShape.size(); k++) {
                        samples[i].currentShape[k].x += alpha * samples[rndIdx].targetShape[k].x;
                        samples[i].currentShape[k].y += alpha * samples[rndIdx].targetShape[k].y;
                        if (samples[rndIdx].isPresent == true) {
                            hits[k] += alpha;
                        }
                        if (hits[k] < 0.1) {
                            isHitAll = false;
                        }
                    }
                    itr++;
                }
	        for (size_t k = 0; k < samples[i].currentShape.size(); k++) {
                    samples[i].currentShape[k].x /= hits[k];
                    samples[i].currentShape[k].y /= hits[k];
                }
            }
        }

        return meanShape;
    }

    std::vector<std::vector<Point2d<float>>> randomlySamplePixelCoordinates(const std::vector<Point2d<float>>& initialShape) {
        Point2d<float> minXY;
        Point2d<float> maxXY;
        minXY.x = 0.0;
        minXY.y = 0.0;
        maxXY.x = 1.0;
        maxXY.y = 1.0;

        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution<double> norm( 0.0, 1.0 );
        std::vector<std::vector<Point2d<float>>> pixelCoordinates(cascadeDepth,std::vector<Point2d<float>(featurePoolSize));
        for (size_t i = 0; i < cascadeDepth; i++) {
            for (size_t j = 0; j < featurePoolSize; j++) {
                pixelCoordinates[i][j].x = norm(mt) * (maxXY.x - minXY.x) + minXY.x;
                pixelCoordinates[i][j].y = norm(mt) * (maxXY.y - minXY.y) + minXY.y;
            }
        }

        return pixelCoordinates;
    }

    RegressionTree makeRegressionTree(std::vector<trainingSample>& samples, const std::vector<std::vector<Point2d<float>>& pixelCoordinates) {
        RegressionTree tree;

        const size_t numSplitsNodes = std::pow(2,treeDepth-1);
        std::vector<Point2d<float>> sums(numSplits*2+1);

        RegressionTree tree;

        for (size_t i = 0; i < samples.size(); i++) {
            samples[i].diffShape.x = samples[i].targetShape.x - samples[i].currentShape.x;
            samples[i].diffShape.y = samples[i].targetShape.y - samples[i].currentShape.y;
            sums[0].x += samples[i].diffShape.x;
            sums[0].y += samples[i].diffShape.y;
        }

        for (size_t i = 0; i < numSplitNodes; i++) {
            
        }
    }

    splitFeature generateSplit(const std::vector<trainingSample>& samples, size_t begin, size_t end, const std::vector<std::vector<Point2d<float>> pixelCoordinates, const std::vector<double>& sum, std::vector<double>& leftSum, std::vector<double>& rightSum ){
        std::vector<splitFeature> feats;
    }
};



private:
    size_t cascadeDepth;
    size_t treeDepth;
    size_t numTreesPerCascadeLevel;
    double nu;
    size_t oversamplingAmount;
    size_t featurePoolSize;
    double lambda;
    size_t numTestSplits;
};

class TrainingSample {
public:
    TrainingSample() {
        oversamplingAmount = 0;
    }
    size_t oversamplingAmount;
    size_t imageIdx;
    std::vector<Rectangle> rect;
    std::vector<float> targetShape;
    // isVisible or isNotOccluded
    std::vector<bool> isPresent;

    std::vector<float> currentShape;
    std::vector<float> diffShape;
    std::vector<unsigned char> featurePixelValues;
};



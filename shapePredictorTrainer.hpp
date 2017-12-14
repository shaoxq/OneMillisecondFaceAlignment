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
        std::vector<std::vector<Point2d>> pixelCoordinates;
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

private:
    size_t cascadeDepth;
    size_t treeDepth;
    size_t numTreesPerCascadeLevel;
    double nu;
    size_t oversamplingAmount;
    size_t featurePoolSize;
    double lambda;
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
    std::vector<bool> isPresent;

    std::vector<float> currentShape;
    std::vector<float> diffShape;
    std::vector<unsigned char> featurePixelValues;
};

class PointTransformAffine {
public:
    PointTransformAffine() {
        m = std::vector<std::vector<double>>(2,std::vector<double>(2,0));
        b = std::vector<double>(2,0);
        s = 1.0;
    }
    // rotation matrix
    std::vector<std::vector<double>> m;
    // scaling factor
    double s;
    // translation vector;
    std::vector<double> b;
};

// This is the implementation of "Appendix D Aligning Two Shapes" of "An Introduction to Active Shape Models" paper.
template <typename T>
PointTransformAffine findAffineTransform (const std::vector<Point2d<T>>& fromPoints, const std::vector<Point2d<T>>& toPoints) {
        assert(fromPoints.size() == toPoints.size());

        size_t pointsNum = fromPoints.size();

        Point2d<double> fromPointsMean;
        Point2d<double> toPointsMean;
        for (size_t i = 0; i < pointsNum; i++) {
            fromPointsMean.x += fromPoints[i].x / pointsNum;
            fromPointsMean.y += fromPoints[i].y / pointsNum;
            toPointsMean.x += toPoints[i].x / pointsNum;
            toPointsMean.y += toPoints[i].y /pointsNum;
        }

        std::vector<Point2d<double>> fromPointsZeroMean(fromPoints.size());
        std::vector<Point2d<double>> toPointsZeroMean(fromPoints.size());
        for (size_t i = 0; i < pointsNum; i++) {
            fromPointsZeroMean[i].x = fromPoints[i].x - fromPointsMean.x;
            fromPointsZeroMean[i].y = fromPoints[i].y - fromPointsMean.y;
            toPointsZeroMean[i].x = toPoints[i].x - toPointsMean.x;
            toPointsZeroMean[i].y = toPoints[i].y - toPointsMean.y;
        }

        double normToPointsZeroMean = 0;
        for (size_t i = 0; i < pointsNum; i++) {
            normToPointsZeroMean += pow(toPointsZeroMean[i].x,2);
            normToPointsZeroMean += pow(toPointsZeroMean[i].y,2);
        }
        normToPointsZeroMean = sqrt(normToPointsZeroMean);

        double dot = 0;
        for (size_t i = 0; i < pointsNum; i++) {
            dot += fromPointsZeroMean[i].x * toPointsZeroMean[i].x;
            dot += fromPointsZeroMean[i].y * toPointsZeroMean[i].y;
        }

        double cross = 0;
        for (size_t i = 0; i < pointsNum; i++) {
            cross +=toPointsZeroMean[i].x*fromPointsZeroMean[i].y;
            cross -=toPointsZeroMean[i].y*fromPointsZeroMean[i].x;
        }

        double a = dot/pow(normToPointsZeroMean,2);
        double b = cross/pow(normToPointsZeroMean,2);

        double scale = 1.0/sqrt(a*a + b*b);
        double angle = -atan2(b, a);
        Point2d<double> translate;
        translate.x = toPointsMean.x - scale * (cos(angle) * fromPointsMean.x - sin(angle) * fromPointsMean.y);
        translate.y = toPointsMean.y - scale * (sin(angle) * fromPointsMean.x + cos(angle) * fromPointsMean.y);

        PointTransformAffine tform;
        tform.m[0][0] = cos(angle);
        tform.m[0][1] = -sin(angle);
        tform.m[1][0] = sin(angle);
        tform.m[1]1[] = cos(angle);
        tform.s = scale;
        tform.b = translate;

        return tform;
     }

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
        size_t numParts = 0;
        for (size_t i = 0; i < objects.size(); i++) {
            if (numParts == 0){
                numParts = objects.parts.size();
                assert(numParts != 0);
            } else {
                assert(objects.parts.size() == numParts);
            }
        }

        assert(numParts != 0);

        std::vector<TrainingSample> samples;
        std::vector<Point2d<float>> initialShape;
        std::vector<std::vector<Point2d>> pixelCoordinates;
    }

    std::vector<Point2d<float>> populateTrainingSampleShapes(const std::vector<fullObjectDetection>& objects, std::vector<TrainingSample>& samples) const {
        std::vector<Point2d<float>> meanShape(objects[0].parts.size());

        for (size_t i = 0; i < objects.size(); i++) {
            TrainingSample sample;
            sample.imageIdx = i;
            sample.rect = objects[i].rect;
            for (size_t j = 0; j < oversamplingAmount; j++) {
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

class TrainingSample {
public:
    TrainingSample() {}
    size_t imageIdx;
    std::vector<Rectangle> rect;
    std::vector<float> targetShape;
    std::vector<float> present;
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

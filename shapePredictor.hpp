#pragma once

#include <vector>
#include "point.hpp"

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

class SplitFeature {
public:
    size_t idx1;
    size_t idx2;
    float threshold;
};

class RegressionTree {
public:
    std::vector<SplitFeature> splits;
    std::vector<Point2d<float>> leafValues;
};


inline void createShapeRelativeEncoding(const std::vector<Point2d<float>>& shape, const std::vector<Point2d<float>>& pixelCoordinates, std::vector<size_t>& anchorIdx, std::vector<Point2d<float>>& deltas) {
    for (size_t i = 1; i < pixelCoordinates.size(); i++) {
        anchorIdx[i] = nearestShapePoint(shape,pixelCoordinates[i];
        deltas[i].x = pixelCoordinates[i].x -shape[anchorIdx[i]].x;
        deltas[i].y = pixelCoordinates[i].y -shape[anchorIdx[i]].y;
    }
}

inline PointTransformAffine unnormalizingTform(const Rectangle& rect) {
    std::vector<Point2d<double>> fromPoints(4);
    fromPoints[0].x = 0.0;
    fromPoints[0].y = 0.0;
    fromPoints[1].x = 1.0;
    fromPoints[1].y = 0.0;
    fromPoints[2].x = 0.0;
    fromPoints[2].y = 1.0;
    fromPoints[3].x = 1.0;
    fromPoints[3].y = 1.0;

    std::vector<Point2d<double>> toPoints(4);
    toPoints[0].x = rect.leftTop.x;
    toPoints[0].y = rect.leftTop.y;
    toPoints[1].x = rect.rightBottom.x;
    toPoints[1].y = rect.leftTop.y;
    toPoints[2].x = rect.leftTop.x;
    toPoints[2].y = rect.rightBottom.y;
    toPoints[3].x = rect.rightBottom.x;
    toPoints[3].y = rect.rightBottom.y;

    return findAffine(fromPoints, toPoints);
}


// This is the implementation of "Appendix D Aligning Two Shapes" of "An Introduction to Active Shape Models" paper.
template <typename T>
PointTransformAffine findAffineTransform(const std::vector<Point2d<T>>& fromPoints, const std::vector<Point2d<T>>& toPoints) {
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

void extractFeaturePixelValues(const Image& img, const Rectangle& rect, const std::vector<Point2d<float>>& currentShape, const std::vector<Point2d<float>>& referehceShape, const std::vector<size_t>& referencePixelAnchorIdx, const std::vector<Point2d<float>>& referencePixelDeltas, std::vector<unsigned char>& featurePixelValues) {
    const PointTransformAffine tform = findAffineTransform(referenceShape, currentShape);
    const Rectangle area;
    area.leftTop.x = 0.0;
    area.leftTop.y = 0.0;
    area.rightBottom.x = img.col-1;
    area.rightBottom.y = img.row-1;

    for (size_t i = 0; i < referencePixelDeltas.size(); i++) {
        Point2d<float> d;
        d.x = referencePixelDeltax[i].x;
        d.y = referencePixelDeltax[i].y;

        size_t idx = referencePixelAnchorIdx[i];

        Point2d<float> p;
        p.x = tform.s*(tform.m[0][0]*d.x+tform.m[0][1]*d.y) + currentShape[idx].x;
        p.y = tform.s*(tform.m[1][0]*d.x+tform.m[1][1]*d.y) + currentShape[idx].y;

        if (0 < p.x && p.x < img.col-1 && 0 < p.y && img.row-1) {
            featurePixelValues[i] = img.data[p.y][p.x][0];
        } else {
            featurePixelValues[i] = 0;
        }
    }
}

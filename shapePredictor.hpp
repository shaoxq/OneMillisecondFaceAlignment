#pragma once

#include <vector>
#include "point.hpp"

inline void createShapeRelativeEncoding(const std::vector<Point2d<float>>& shape, const std::vector<Point2d<float>>& pixelCoordinates, std::vector<size_t>& anchorIdx, std::vector<Point2d<float>>& deltas) {
    for (size_t i = 1; i < pixelCoordinates.size(); i++) {
        anchorIdx[i] = nearestShapePoint(shape,pixelCoordinates[i];
        deltas[i].x = pixelCoordinates[i].x -shape[anchorIdx[i]].x;
        deltas[i].y = pixelCoordinates[i].y -shape[anchorIdx[i]].y;
    }
}


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

void extractFeaturePixelValues(const Image& img, const Rectangle& rect, const std::vector<Point2d<float>>& currentShape, const std::vector<Point2d<float>>& referehceShape, const std::vector<size_t>& referencePixelAnchorIdx, const std::vector<Point2d<float>>& referencePixelDeltas, std::vector<unsigned char>& featurePixelValues) {
    
}

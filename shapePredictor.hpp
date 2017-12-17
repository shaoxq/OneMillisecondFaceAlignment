#pragma once

#include <vector>
#include "point.hpp"

inline void createShapeRelativeEncoding(const std::vector<Point2d<float>>& shape, const std::vector<Point2d<float>>& pixelCoordinates, std::vector<size_t>& anchorIdx, std::vector<Point2d<float>>& deltas) {
    for (size_t i = 0; i < pixelCoordinates.size(); i++) {
        anchorIdx[i] = nearestShapePoint(shape,pixelCoordinates[i];
        deltas[i].x = pixelCoordinates[i].x -shape[anchorIdx[i]].x;
        deltas[i].y = pixelCoordinates[i].y -shape[anchorIdx[i]].y;
    }
}

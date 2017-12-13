#pragma once

#include <vector>
#include "point.hpp"
#include "rectangle.hpp"

class FullObjectDetection {
public:
    FullObjectDetection() {}
    std::vector<Rectangle> rect;
    std::vector<Point2d<double>> parts;
};

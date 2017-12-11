#pragma once

#include <vector>
#include "point.hpp"
#include "rectangle.hpp"

class FullObjectDetection {
public:
    FullObjectDetection() {}
private:
    std::vector<Rectangle> rect;
    std::vector<Point2d<double>> parts;
};

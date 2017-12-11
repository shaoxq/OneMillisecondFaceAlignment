#pragma once

#include "point.hpp"

class Rectangle {
public:
    Rectangle() {
    }
    Point<size_t> leftTop;
    Point<size_t> rightBottom;
};

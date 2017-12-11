#pragma once

template <typename T>
class Point2d {
public:
    Point2d() {
        x = 0;
        y = 0;
    }
    T x;
    T y;
};

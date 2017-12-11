#pragma once

#include <vector>

class FullObjectDetection {
public:
    FullObjectDetection() {}
private:
    std::vector<size_t> rect;
    std::vector<std::vector<double>> parts;
};

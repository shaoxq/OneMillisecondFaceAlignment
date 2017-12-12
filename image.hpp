#pragma once

#include <vector>

class Image {
public:
    Image() {}
    Image(size_t r, size_t c, size_t ch) {
        row = r;
        col = c;
        channel = ch;
    }
    std::vector<std::vector<std::vector<unsigned char>>> data;
    size_t row;
    size_t col;
    size_t channel
};

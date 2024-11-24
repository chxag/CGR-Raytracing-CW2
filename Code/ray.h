#pragma once

#include <vector>

class Ray {
  public:
    std::vector<float> origin;
    std::vector<float> direction;
    Ray(std::vector<float>& origin,std::vector<float>& direction)
        : origin(origin), direction(direction) {}
};
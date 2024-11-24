#ifndef LIGHT_H
#define LIGHT_H

#include <vector>
#include <string>

struct Light
{
    std::string light_type;
    std::vector<float> light_position;
    std::vector<float> intensity;

    Light(const std::string light_type, const std::vector<float> light_position, const std::vector<float> intensity)
        : light_type(light_type), light_position(light_position), intensity(intensity) {}
};

#endif
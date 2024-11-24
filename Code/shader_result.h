#ifndef SHADER_RESULT_H
#define SHADER_RESULT_H

#include <vector>
#include "material.h"

struct ShaderResult
{
    std::vector<float> color;
    bool intersected;
    std::vector<float> intersection_point;
    Material intersected_material;
    std::vector<float> normal;
    std::vector<float> uv_coordinates;
};

#endif

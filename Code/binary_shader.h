#ifndef BINARY_SHADER_H
#define BINARY_SHADER_H

#include <vector>
#include "material.h"
#include "ray.h"
#include "sphere.h"
#include "cylinder.h"
#include "triangle.h"
#include "shader_result.h"

class BinaryShader
{
public:
    static ShaderResult calculateColor(const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Cylinder> &cylinders, const std::vector<Triangle> &triangles, const std::vector<float> &backgroundcolor);
};

#endif
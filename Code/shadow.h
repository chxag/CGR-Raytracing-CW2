#ifndef SHADOW_H
#define SHADOW_H

#include <vector>
#include "light.h"
#include "sphere.h"
#include "cylinder.h"
#include "triangle.h"

class Shadow
{
public:
    static bool isInShadow(const std::vector<float>& point, const Light& light, const std::vector<Sphere>& spheres, const std::vector<Cylinder>& cylinders, const std::vector<Triangle>& triangles);
};

#endif
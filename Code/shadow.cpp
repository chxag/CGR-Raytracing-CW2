#include "shadow.h"
#include "vector_utils.h"

bool Shadow::isInShadow(const std::vector<float>& point, const Light& light, const std::vector<Sphere>& spheres, const std::vector<Cylinder>& cylinders, const std::vector<Triangle>& triangles)
{
    std::vector<float> lightDir = {
        light.light_position[0] - point[0],
        light.light_position[1] - point[1],
        light.light_position[2] - point[2]
    };
    normalize(lightDir);

    float shadowBias = 0.01f;
    std::vector<float> shadowRayOrigin = {
        point[0] + shadowBias * lightDir[0],
        point[1] + shadowBias * lightDir[1],
        point[2] + shadowBias * lightDir[2]
    };

    Ray shadowRay(shadowRayOrigin, lightDir);

    float t;
    for (const auto &sphere : spheres)
    {
        if (sphere.intersectSphere(shadowRay, t))
        {
            return true;
        }
    }

    for (const auto &cylinder : cylinders)
    {
        if (cylinder.intersectCylinder(shadowRay, t))
        {
            return true;
        }
    }

    for (const auto &triangle : triangles)
    {
        if (triangle.intersectTriangle(shadowRay, t))
        {
            return true;
        }
    }

    return false;
}
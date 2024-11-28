#include "shadow.h"
#include "vector_utils.h"

// Check if the point is in shadow
bool Shadow::isInShadow(const std::vector<float>& point, const Light& light, const std::vector<Sphere>& spheres, const std::vector<Cylinder>& cylinders, const std::vector<Triangle>& triangles)
{
    // Calculate the vector from the point to the light
    std::vector<float> lightDir = {
        light.light_position[0] - point[0],
        light.light_position[1] - point[1],
        light.light_position[2] - point[2]
    };
    normalize(lightDir);

    // Calculate the bias for the shadow ray
    float shadowBias = 0.01f;
    // Calculate the origin of the shadow ray
    std::vector<float> shadowRayOrigin = {
        point[0] + shadowBias * lightDir[0],
        point[1] + shadowBias * lightDir[1],
        point[2] + shadowBias * lightDir[2]
    };

    // Initialize the shadow ray
    Ray shadowRay(shadowRayOrigin, lightDir);

    float t;
    // Check for intersection with spheres
    for (const auto &sphere : spheres)
    {
        if (sphere.intersectSphere(shadowRay, t))
        {
            return true;
        }
    }

    // Check for intersection with cylinders
    for (const auto &cylinder : cylinders)
    {
        if (cylinder.intersectCylinder(shadowRay, t))
        {
            return true;
        }
    }

    // Check for intersection with triangles
    for (const auto &triangle : triangles)
    {
        if (triangle.intersectTriangle(shadowRay, t))
        {
            return true;
        }
    }

    // If there is no intersection, the point is not in shadow
    return false;
}
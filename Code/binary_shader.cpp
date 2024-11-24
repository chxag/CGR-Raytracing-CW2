#include "binary_shader.h"
#include <limits>

ShaderResult BinaryShader::calculateColor(const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Cylinder> &cylinders, const std::vector<Triangle> &triangles, const std::vector<float> &backgroundcolor)
{
    float closestT = std::numeric_limits<float>::max();
    bool intersected = false;
    std::vector<float> intersected_color = backgroundcolor;

    // Check intersection with spheres
    for (const auto &sphere : spheres)
    {
        float t;
        if (sphere.intersectSphere(ray, t) && t < closestT)
        {
            closestT = t;
            intersected = true;
            intersected_color = {1.0f, 0.0f, 0.0f}; // Hardcoded color for binary mode
        }
    }

    // Check intersection with cylinders
    for (const auto &cylinder : cylinders)
    {
        float t;
        if (cylinder.intersectCylinder(ray, t) && t < closestT)
        {
            closestT = t;
            intersected = true;
            intersected_color = {1.0f, 0.0f, 0.0f};
        }
    }

    // Check intersection with triangles
    for (const auto &triangle : triangles)
    {
        float t;
        if (triangle.intersectTriangle(ray, t) && t < closestT)
        {
            closestT = t;
            intersected = true;
            intersected_color = {1.0f, 0.0f, 0.0f};
        }
    }

    return {intersected_color, intersected, {0.0f, 0.0f, 0.0f}, Material(), {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
}
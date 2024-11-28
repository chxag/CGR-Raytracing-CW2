#include "binary_shader.h"
#include <limits>

// Function to calculate the color of a ray intersecting with the objects in the scene
ShaderResult BinaryShader::calculateColor(const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Cylinder> &cylinders, const std::vector<Triangle> &triangles, const std::vector<float> &backgroundcolor)
{
    float closestT = std::numeric_limits<float>::max(); // Initialize closest intersection distance to a large value
    bool intersected = false; // Initialize intersected flag to false
    std::vector<float> intersected_color = backgroundcolor; // Initialize intersected color to background color 

    
    // Check intersection with spheres
    for (const auto &sphere : spheres)
    {
        float t;
        //If the ray intersects the sphere and is closer
        if (sphere.intersectSphere(ray, t) && t < closestT)
        {
            closestT = t; // Update closest intersection distance
            intersected = true; // Set intersected flag to true
            intersected_color = {1.0f, 0.0f, 0.0f}; // Set intersected color to red
        }
    }

    
    // Check intersection with cylinders
    for (const auto &cylinder : cylinders)
    {
        float t;
        //If the ray intersects the cylinder and is closer
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
        //If the ray intersects the triangle and is closer
        if (triangle.intersectTriangle(ray, t) && t < closestT)
        {
            closestT = t; // Update closest intersection distance
            intersected = true;
            intersected_color = {1.0f, 0.0f, 0.0f};
        }
    }

    // Return intersection color, flag, and default material properties
    return {intersected_color, intersected, {0.0f, 0.0f, 0.0f}, Material(), {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
}
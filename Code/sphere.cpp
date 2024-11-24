#include "sphere.h"
#include <cmath>
#include <limits>
#include "vector_utils.h"

Sphere::Sphere(const std::vector<float> center, float radius, Material material)
    : center(center), radius(radius), material(material) {}

float Sphere::find_root(const Ray &ray) const
{
    std::vector<float> oc = {ray.origin[0] - center[0], ray.origin[1] - center[1], ray.origin[2] - center[2]};
    float a = ray.direction[0] * ray.direction[0] + ray.direction[1] * ray.direction[1] + ray.direction[2] * ray.direction[2];
    float b = 2.0f * (ray.direction[0] * oc[0] + ray.direction[1] * oc[1] + ray.direction[2] * oc[2]);
    float c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0)
    {
        return -1.0f;
    }
    
    float root1 = (-b - sqrt(discriminant)) / (2.0f * a);
    float root2 = (-b + sqrt(discriminant)) / (2.0f * a);

    if (root1 > 0 && root2 > 0)
    {
        return std::min(root1, root2);
    }
    else if (root1 > 0)
    {
        return root1;
    }
    else if (root2 > 0)
    {
        return root2;
    }
    return -1.0f;
}

std::vector<float> Sphere::getUV(const std::vector<float> &intersectionPoint) const
{
    std::vector<float> localPoint = {
        intersectionPoint[0] - center[0],
        intersectionPoint[1] - center[1],
        intersectionPoint[2] - center[2]};
        
    normalize(localPoint);

    float u = 0.5f + atan2(localPoint[2], localPoint[0]) / (2.0f * M_PI);
    float v = 0.5f - asin(localPoint[1]) / M_PI;

    return {u, v};
}

bool Sphere::intersectSphere(const Ray &ray, float &t) const
{
    float root = find_root(ray);
    if (root < 0)
    {
        return false;
    }
    t = root;
    return true;
}
#ifndef SPHERE_H
#define SPHERE_H

#include "ray.h"
#include "material.h"
#include <vector>

class Sphere{
    public:
        Sphere(const std::vector<float> center, float radius, Material material);
        bool intersectSphere(const Ray& ray, float& t) const;
        float find_root(const Ray& ray) const;
        std::vector<float> getUV(const std::vector<float>& intersectionPoint) const;
        std::vector<float> center;
        float radius;
        Material material;
    private:

};


#endif
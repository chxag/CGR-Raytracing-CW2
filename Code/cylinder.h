#ifndef CYLINDER_H
#define CYLINDER_H

#include "ray.h"
#include <vector>
#include "material.h"

class Cylinder {
    public:

        Cylinder(const std::vector<float>& center, float radius, const std::vector<float>& axis, float height, Material material);
        bool intersectCylinder(const Ray& ray, float& t) const;
        // float find_root(const Ray& ray) const;
        // bool check_cap_intersection(const Ray& ray, float& t, const std::vector<float>& cap_center, const std::vector<float>& cap_normal) const;
        std::vector<float> getUV(const std::vector<float>& intersectionPoint) const;
        std::vector<float> center;
        float radius;
        std::vector<float> axis;
        float height;
        Material material;
    
    private:
};

#endif
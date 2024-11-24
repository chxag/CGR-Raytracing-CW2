#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "ray.h"
#include <vector>
#include "material.h"

class Triangle {
    public:

        Triangle(const std::vector<float>& v0, 
                 const std::vector<float>& v1, 
                 const std::vector<float>& v2, 
                 Material material);
        bool intersectTriangle(const Ray& ray, float& t) const;
        std::vector<float> getUV(const std::vector<float>& intersectionPoint) const;
        std::vector<float> v0;
        std::vector<float> v1;
        std::vector<float> v2;
        Material material;
    
    private:
};
#endif
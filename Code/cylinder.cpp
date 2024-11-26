#include "cylinder.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>
#include "vector_utils.h"

Cylinder::Cylinder(const std::vector<float> &center, float radius, const std::vector<float> &axis, float height, Material material)
    : center(center), radius(radius), axis(axis), height(height*2), material(material){
    float axis_length = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    this->axis = {axis[0] / axis_length, axis[1] / axis_length, axis[2] / axis_length};

}

bool Cylinder::intersectCylinder(const Ray& ray, float& t) const {
    std::vector<float> oc = {ray.origin[0] - center[0],
                             ray.origin[1] - center[1],
                             ray.origin[2] - center[2]};
     float d_dot_a = ray.direction[0] * axis[0] +
                    ray.direction[1] * axis[1] +
                    ray.direction[2] * axis[2];
    float oc_dot_a = oc[0] * axis[0] +
                     oc[1] * axis[1] +
                     oc[2] * axis[2];
    float a = ray.direction[0] * ray.direction[0] +
              ray.direction[1] * ray.direction[1] +
              ray.direction[2] * ray.direction[2] -
              d_dot_a * d_dot_a;
    float half_b = oc[0] * ray.direction[0] +
                   oc[1] * ray.direction[1] +
                   oc[2] * ray.direction[2] -
                   oc_dot_a * d_dot_a;
    float c = oc[0] * oc[0] +
              oc[1] * oc[1] +
              oc[2] * oc[2] -
              oc_dot_a * oc_dot_a -
              radius * radius;
    
    float t_side = std::numeric_limits<float>::max();
    float discriminant = static_cast<double>(half_b) * half_b - a * c;

    if (discriminant >= 0){
        float sqrt_discriminant = sqrt(discriminant);
        float root1 = (-half_b - sqrt_discriminant) / a;
        float root2 = (-half_b + sqrt_discriminant) / a;
        
        for (float root : {root1, root2}){
            if (root > 0.0f){
                std::vector<float> point = {
                    ray.origin[0] + root * ray.direction[0],
                    ray.origin[1] + root * ray.direction[1],
                    ray.origin[2] + root * ray.direction[2]
                };
                std::vector<float> toPoint = {
                    point[0] - center[0],
                    point[1] - center[1],
                    point[2] - center[2]
                };
                float projection = toPoint[0] * axis[0] +
                                   toPoint[1] * axis[1] +
                                   toPoint[2] * axis[2];
                if (projection >= -height / 2 && projection <= height / 2){
                    t_side = std::min(t_side, root);
                }
            }
        }
    }

    float t_caps = std::numeric_limits<float>::max();

    for (float h : {-height / 2, height / 2}){
        std::vector<float> capCenter = {center[0] + axis[0] * h,
                                        center[1] + axis[1] * h,
                                        center[2] + axis[2] * h};
        
        float denom = ray.direction[0] * axis[0] +
                      ray.direction[1] * axis[1] +
                      ray.direction[2] * axis[2];
        if (fabs(denom) > 1e-6){
            std::vector<float> oc = {capCenter[0] - ray.origin[0],
                                     capCenter[1] - ray.origin[1],
                                     capCenter[2] - ray.origin[2]};
            float t_cap = (oc[0] * axis[0] +
                           oc[1] * axis[1] +
                           oc[2] * axis[2]) / denom;
            
            if (t_cap > 0.0f){
                std::vector<float> p = {
                    ray.origin[0] + t_cap * ray.direction[0],
                    ray.origin[1] + t_cap * ray.direction[1],
                    ray.origin[2] + t_cap * ray.direction[2]
                };

                float dist_sq = (p[0] - capCenter[0]) * (p[0] - capCenter[0]) +
                                (p[1] - capCenter[1]) * (p[1] - capCenter[1]) +
                                (p[2] - capCenter[2]) * (p[2] - capCenter[2]);
                if (dist_sq <= radius * radius){
                    t_caps = std::min(t_caps, t_cap);
                }
            }
        }
    }

    t = std::min(t_side, t_caps);
    return t < std::numeric_limits<float>::max();
}

std::vector<float> Cylinder::getUV(const std::vector<float>& intersectionPoint) const{
    std::vector<float> toPoint = {
        intersectionPoint[0] - center[0],
        intersectionPoint[1] - center[1],
        intersectionPoint[2] - center[2]
    };

    float projection = toPoint[0] * axis[0] +
                       toPoint[1] * axis[1] +
                       toPoint[2] * axis[2];

    float v = (projection + height / 2) / height;

    std::vector<float> axisComponent = {
        axis[0] * projection,
        axis[1] * projection,
        axis[2] * projection
    };

    std::vector<float> projectedToPoint = {
        toPoint[0] - axisComponent[0],
        toPoint[1] - axisComponent[1],
        toPoint[2] - axisComponent[2]
    };

    float angle = atan2(projectedToPoint[2], projectedToPoint[0]);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    float u = angle / (2 * M_PI);

    return {u,v};
}

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

// Intersect the ray with the cylinder
bool Cylinder::intersectCylinder(const Ray& ray, float& t) const {

    // Calculate the vector from the ray origin to the cylinder center
    std::vector<float> oc = {ray.origin[0] - center[0],
                             ray.origin[1] - center[1],
                             ray.origin[2] - center[2]};

    // Calculate the dot product of the ray direction and the cylinder axis
     float d_dot_a = ray.direction[0] * axis[0] +
                    ray.direction[1] * axis[1] +
                    ray.direction[2] * axis[2];
    
    // Calculate the dot product of the vector from the ray origin to the cylinder center and the cylinder axis
    float oc_dot_a = oc[0] * axis[0] +
                     oc[1] * axis[1] +
                     oc[2] * axis[2];
    
    // Calculate the squared length of the ray direction
    float a = ray.direction[0] * ray.direction[0] +
              ray.direction[1] * ray.direction[1] +
              ray.direction[2] * ray.direction[2] -
              d_dot_a * d_dot_a;
    
    // Calculate the dot product of the vector from the ray origin to the cylinder center and the ray direction
    float half_b = oc[0] * ray.direction[0] +
                   oc[1] * ray.direction[1] +
                   oc[2] * ray.direction[2] -
                   oc_dot_a * d_dot_a;
    
    // Calculate the squared length of the vector from the ray origin to the cylinder center
    float c = oc[0] * oc[0] +
              oc[1] * oc[1] +
              oc[2] * oc[2] -
              oc_dot_a * oc_dot_a -
              radius * radius;
    
    // Initialize the closest intersection distance
    float t_side = std::numeric_limits<float>::max();

    // Calculate the discriminant
    float discriminant = static_cast<double>(half_b) * half_b - a * c;

    // If the discriminant is non-negative
    if (discriminant >= 0){
        // Calculate the square root of the discriminant
        float sqrt_discriminant = sqrt(discriminant);
        float root1 = (-half_b - sqrt_discriminant) / a;
        float root2 = (-half_b + sqrt_discriminant) / a;
        
        // Check both roots
        for (float root : {root1, root2}){
            // If the root is positive
            if (root > 0.0f){
                // Calculate the intersection point
                std::vector<float> point = {
                    ray.origin[0] + root * ray.direction[0],
                    ray.origin[1] + root * ray.direction[1],
                    ray.origin[2] + root * ray.direction[2]
                };

                // Calculate the vector from the intersection point to the cylinder center
                std::vector<float> toPoint = {
                    point[0] - center[0],
                    point[1] - center[1],
                    point[2] - center[2]
                };
                //Check if intersection point is within the height of the cylinder
                float projection = toPoint[0] * axis[0] +
                                   toPoint[1] * axis[1] +
                                   toPoint[2] * axis[2];
                if (projection >= -height / 2 && projection <= height / 2){
                    t_side = std::min(t_side, root);
                }
            }
        }
    }

    // Initialize the closest intersection distance for the caps
    float t_caps = std::numeric_limits<float>::max();

    // Check both caps
    for (float h : {-height / 2, height / 2}){
        // Calculate the center of the cap
        std::vector<float> capCenter = {center[0] + axis[0] * h,
                                        center[1] + axis[1] * h,
                                        center[2] + axis[2] * h};
        
        // Calculate the denominator
        float denom = ray.direction[0] * axis[0] +
                      ray.direction[1] * axis[1] +
                      ray.direction[2] * axis[2];
        // If the denominator is not zero
            if (fabs(denom) > 1e-6){
            // Calculate the vector from the ray origin to the cap center
            std::vector<float> oc = {capCenter[0] - ray.origin[0],
                                     capCenter[1] - ray.origin[1],
                                     capCenter[2] - ray.origin[2]};
            // Calculate the intersection distance
            float t_cap = (oc[0] * axis[0] +
                           oc[1] * axis[1] +
                           oc[2] * axis[2]) / denom;
            
            // If the intersection distance is positive
            if (t_cap > 0.0f){
                // Calculate the intersection point
                std::vector<float> p = {
                    ray.origin[0] + t_cap * ray.direction[0],
                    ray.origin[1] + t_cap * ray.direction[1],
                    ray.origin[2] + t_cap * ray.direction[2]
                };

                // Calculate the squared distance between the intersection point and the cap center
                float dist_sq = (p[0] - capCenter[0]) * (p[0] - capCenter[0]) +
                                (p[1] - capCenter[1]) * (p[1] - capCenter[1]) +
                                (p[2] - capCenter[2]) * (p[2] - capCenter[2]);
                // If the distance is less than or equal to the radius squared
                if (dist_sq <= radius * radius){
                    // Update the closest intersection distance for the caps
                    t_caps = std::min(t_caps, t_cap);
                }
            }
        }
    }

    // Update the closest intersection distance
    t = std::min(t_side, t_caps);
    // Return true if there is an intersection
    return t < std::numeric_limits<float>::max();
}

// Calculate the UV coordinates of the intersection point
std::vector<float> Cylinder::getUV(const std::vector<float>& intersectionPoint) const{

    // Calculate the vector from the intersection point to the cylinder center
    std::vector<float> toPoint = {
        intersectionPoint[0] - center[0],
        intersectionPoint[1] - center[1],
        intersectionPoint[2] - center[2]
    };

    // Calculate the projection of the vector from the intersection point to the cylinder center onto the cylinder axis
    float projection = toPoint[0] * axis[0] +
                       toPoint[1] * axis[1] +
                       toPoint[2] * axis[2];

    // Calculate the V coordinate
    float v = (projection + height / 2) / height;

    // Calculate the axis component of the vector from the intersection point to the cylinder center
    std::vector<float> axisComponent = {
        axis[0] * projection,
        axis[1] * projection,
        axis[2] * projection
    };

    // Calculate the vector from the intersection point to the cylinder center minus the axis component
    std::vector<float> projectedToPoint = {
        toPoint[0] - axisComponent[0],
        toPoint[1] - axisComponent[1],
        toPoint[2] - axisComponent[2]
    };

    // Calculate the angle for U coordinate
    float angle = atan2(projectedToPoint[2], projectedToPoint[0]);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    // Normalize the U coordinate to the range [0, 1]
    float u = angle / (2 * M_PI);

    return {u,v};
}

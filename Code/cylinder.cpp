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
     double d_dot_a = ray.direction[0] * axis[0] +
                    ray.direction[1] * axis[1] +
                    ray.direction[2] * axis[2];
    double oc_dot_a = oc[0] * axis[0] +
                     oc[1] * axis[1] +
                     oc[2] * axis[2];
    double a = ray.direction[0] * ray.direction[0] +
              ray.direction[1] * ray.direction[1] +
              ray.direction[2] * ray.direction[2] -
              d_dot_a * d_dot_a;
    double half_b = oc[0] * ray.direction[0] +
                   oc[1] * ray.direction[1] +
                   oc[2] * ray.direction[2] -
                   oc_dot_a * d_dot_a;
    double c = oc[0] * oc[0] +
              oc[1] * oc[1] +
              oc[2] * oc[2] -
              oc_dot_a * oc_dot_a -
              radius * radius;
    double discriminant = static_cast<double>(half_b) * half_b - a * c;

    double sqrt_discriminant = sqrt(discriminant);
    double root1 = (-half_b - sqrt_discriminant) / a;
    double root2 = (-half_b + sqrt_discriminant) / a;

    float t_min = std::numeric_limits<float>::max();
    for (float t_candidate : {root1, root2}) {
        if (t_candidate > 0) {
            std::vector<float> point = {
                ray.origin[0] + t_candidate * ray.direction[0],
                ray.origin[1] + t_candidate * ray.direction[1],
                ray.origin[2] + t_candidate * ray.direction[2]};
            std::vector<float> to_point = {
                point[0] - center[0],
                point[1] - center[1],
                point[2] - center[2]};
            float projection = to_point[0] * axis[0] +
                               to_point[1] * axis[1] +
                               to_point[2] * axis[2];
            if (projection >= -height / 2 && projection <= height / 2 && t_candidate < t_min) {
                t_min = t_candidate;
            }
        }
    }

    if (t_min < std::numeric_limits<float>::max()) {
        t = t_min;
        return true;
    }

    return false;
    
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

// bool Cylinder::intersectCylinder(const Ray &ray, float &t) const
// {

//     float side_t = std::numeric_limits<float>::max();
//     float cap_t = std::numeric_limits<float>::max();

//     float side_root = find_root(ray);

//     if (side_root >= 0)
//     {

//         std::vector<float> point = {
//             ray.origin[0] + side_root * ray.direction[0],
//             ray.origin[1] + side_root * ray.direction[1],
//             ray.origin[2] + side_root * ray.direction[2]};
//         std::vector<float> to_point = {
//             point[0] - center[0],
//             point[1] - center[1],
//             point[2] - center[2]};
//         float projection = to_point[0] * axis[0] +
//                            to_point[1] * axis[1] +
//                            to_point[2] * axis[2];
//         if (projection >= -height && projection <= height)
//         {
//             side_t = side_root;
//         }
//     }

//     std::vector<float> bottom_center = center;
//     float bottom_cap_t = std::numeric_limits<float>::max();
//     if (check_cap_intersection(ray, bottom_cap_t, bottom_center, {-axis[0], -axis[1], -axis[2]}))
//     {
//         cap_t = std::min(cap_t, bottom_cap_t);
//     }

//     std::vector<float> top_center = {center[0] + axis[0] * height, center[1] + axis[1] * height, center[2] + axis[2] * height};
//     float top_cap_t = std::numeric_limits<float>::max();
//     if (check_cap_intersection(ray, top_cap_t, top_center, axis))
//     {
//         cap_t = std::min(cap_t, top_cap_t);
//     }

//     t = std::numeric_limits<float>::max();
//     if(side_t > 0.0f){
//         t = side_t;
//     } 
//     if(cap_t > 0.0f && cap_t < t){
//         t = cap_t;
//     }
//     return t < std::numeric_limits<float>::max();
    

// }

// bool Cylinder::check_cap_intersection(const Ray &ray, float &t, const std::vector<float> &cap_center, const std::vector<float> &cap_normal) const
// {
//     float denom = ray.direction[0] * cap_normal[0] +
//                   ray.direction[1] * cap_normal[1] +
//                   ray.direction[2] * cap_normal[2];

//     if (fabs(denom) < 1e-5)
//     {
//         return false;
//     }
//     float t_cap = ((cap_center[0] - ray.origin[0]) * cap_normal[0] +
//                    (cap_center[1] - ray.origin[1]) * cap_normal[1] +
//                    (cap_center[2] - ray.origin[2]) * cap_normal[2]) /
//                   denom;
//     if (t_cap < 0)
//     {
//         return false;
//     }

//     std::vector<float> point = {
//         ray.origin[0] + t_cap * ray.direction[0],
//         ray.origin[1] + t_cap * ray.direction[1],
//         ray.origin[2] + t_cap * ray.direction[2]};

//     float dist_squared = (point[0] - cap_center[0]) * (point[0] - cap_center[0]) +
//                          (point[1] - cap_center[1]) * (point[1] - cap_center[1]) +
//                          (point[2] - cap_center[2]) * (point[2] - cap_center[2]);
//     if (dist_squared > radius * radius + 1e-5)
//     {
//         return false;
//     }

//     t = t_cap;
//     return true;
// }

// float Cylinder::find_root(const Ray &ray) const
// {
//     std::vector<float> oc = {ray.origin[0] - center[0],
//                              ray.origin[1] - center[1],
//                              ray.origin[2] - center[2]};

//     double d_dot_a = ray.direction[0] * axis[0] +
//                     ray.direction[1] * axis[1] +
//                     ray.direction[2] * axis[2];
//     double oc_dot_a = oc[0] * axis[0] +
//                      oc[1] * axis[1] +
//                      oc[2] * axis[2];
//     double a = ray.direction[0] * ray.direction[0] +
//               ray.direction[1] * ray.direction[1] +
//               ray.direction[2] * ray.direction[2] -
//               d_dot_a * d_dot_a;
//     double half_b = oc[0] * ray.direction[0] +
//                    oc[1] * ray.direction[1] +
//                    oc[2] * ray.direction[2] -
//                    oc_dot_a * d_dot_a;
//     double c = oc[0] * oc[0] +
//               oc[1] * oc[1] +
//               oc[2] * oc[2] -
//               oc_dot_a * oc_dot_a -
//               radius * radius;
//     double discriminant = static_cast<double>(half_b) * half_b - a * c;

//     if (discriminant < 1e-6 || fabs(a) < 1e-6)
//     {
//         return -1.0f;
//     }

//     double sqrt_discriminant = sqrt(discriminant);
//     double root1 = (-half_b - sqrt_discriminant) / a;
//     double root2 = (-half_b + sqrt_discriminant) / a;
//     if (root1 > 0 && root2 > 0)
//     {
//         return std::min(root1, root2);
//     } else if (root1 > 0)
//     {
//         return root1;
//     } else if (root2 > 0)
//     {
//         return root2;
//     }
//     return -1.0f;
// }
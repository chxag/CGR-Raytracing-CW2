#include "triangle.h"
#include <cmath>

Triangle::Triangle(const std::vector<float>& v0, 
                   const std::vector<float>& v1, 
                   const std::vector<float>& v2,
                   Material material)
    : v0(v0), v1(v1), v2(v2), material(material) {}

bool Triangle::intersectTriangle(const Ray& ray, float& t) const{

    std::vector<float> e1 = {v1[0] - v0[0], v1[1] - v0[1],  v1[2] - v0[2]};
    std::vector<float> e2 = {v2[0] - v0[0], v2[1] - v0[1],  v2[2] - v0[2]};

    std::vector<float> p = {ray.direction[1] * e2[2] - ray.direction[2] * e2[1],
                            ray.direction[2] * e2[0] - ray.direction[0] * e2[2],
                            ray.direction[0] * e2[1] - ray.direction[1] * e2[0]};

    float a = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];

    if (std::fabs(a) < 1e-8){
        return false;
    }

    float f = 1.0 / a;

    std::vector<float> s = {ray.origin[0] - v0[0], ray.origin[1] - v0[1], ray.origin[2] - v0[2]};

    float u = f * (s[0] * p[0] +  s[1] * p[1] + s[2] * p[2]);

    if (u < 0.0 || u > 1.0){
        return false;
    }

    std::vector<float> q = {s[1] * e1[2] - s[2] * e1[1],
                            s[2] * e1[0] - s[0] * e1[2],
                            s[0] * e1[1] - s[1] * e1[0]};
    double v = f * (ray.direction[0] * q[0] + ray.direction[1] * q[1] + ray.direction[2] * q[2]);

    if (v < 0.0 || u + v > 1.0){
        return false;
    }

    t = f * e2[0] * q[0] + f * e2[1] * q[1] + f * e2[2] * q[2];

    return t > 1e-8;
}

std::vector<float> Triangle::getUV(const std::vector<float>& intersectionPoint) const{
    std::vector<float> e1 = {v1[0] - v0[0], v1[1] - v0[1],  v1[2] - v0[2]};
    std::vector<float> e2 = {v2[0] - v0[0], v2[1] - v0[1],  v2[2] - v0[2]};

    std::vector<float> v0p = {intersectionPoint[0] - v0[0], intersectionPoint[1] - v0[1], intersectionPoint[2] - v0[2]};

    float d_e1_e1 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
    float d_e1_e2 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
    float d_e2_e2 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
    float d_e1_v0p = e1[0] * v0p[0] + e1[1] * v0p[1] + e1[2] * v0p[2];
    float d_e2_v0p = e2[0] * v0p[0] + e2[1] * v0p[1] + e2[2] * v0p[2];

    float denom = d_e1_e1 * d_e2_e2 - d_e1_e2 * d_e1_e2;

    float u = (d_e1_e2 * d_e2_v0p - d_e2_e2 * d_e1_v0p) / denom;
    float v = (d_e1_e1 * d_e2_v0p - d_e1_e2 * d_e1_v0p) / denom;

    return {u, v};
}
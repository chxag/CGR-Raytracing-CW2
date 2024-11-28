#include "triangle.h"
#include <cmath>

Triangle::Triangle(const std::vector<float>& v0, 
                   const std::vector<float>& v1, 
                   const std::vector<float>& v2,
                   Material material)
    : v0(v0), v1(v1), v2(v2), material(material) {}

// Intersect the triangle with the ray
bool Triangle::intersectTriangle(const Ray& ray, float& t) const{

    // Calculate the edge vectors
    std::vector<float> e1 = {v1[0] - v0[0], v1[1] - v0[1],  v1[2] - v0[2]};
    std::vector<float> e2 = {v2[0] - v0[0], v2[1] - v0[1],  v2[2] - v0[2]};

    // Calculate the plane normal
    std::vector<float> p = {ray.direction[1] * e2[2] - ray.direction[2] * e2[1],
                            ray.direction[2] * e2[0] - ray.direction[0] * e2[2],
                            ray.direction[0] * e2[1] - ray.direction[1] * e2[0]};

    // Calculate the determinant
    float a = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];

    // If the determinant is zero, the ray is parallel to the triangle
    if (std::fabs(a) < 1e-8){
        return false;
    }

    // Calculate the inverse of the determinant
    float f = 1.0 / a;

    // Calculate the vector from the ray origin to the triangle vertex
    std::vector<float> s = {ray.origin[0] - v0[0], ray.origin[1] - v0[1], ray.origin[2] - v0[2]};

    // Calculate the u coordinate
    float u = f * (s[0] * p[0] +  s[1] * p[1] + s[2] * p[2]);

    // If the u coordinate is outside the triangle, there is no intersection
    if (u < 0.0 || u > 1.0){
        return false;
    }

    // Calculate the vector q
    std::vector<float> q = {s[1] * e1[2] - s[2] * e1[1],
                            s[2] * e1[0] - s[0] * e1[2],
                            s[0] * e1[1] - s[1] * e1[0]};

    // Calculate the v coordinate
    double v = f * (ray.direction[0] * q[0] + ray.direction[1] * q[1] + ray.direction[2] * q[2]);

    // If the v coordinate is outside the triangle, there is no intersection
    if (v < 0.0 || u + v > 1.0){
        return false;
    }

    // Calculate the t coordinate   
    t = f * e2[0] * q[0] + f * e2[1] * q[1] + f * e2[2] * q[2];

    // If the t coordinate is outside the triangle, there is no intersection
    return t > 1e-8;
}

// Get the uv coordinates of the intersection point
std::vector<float> Triangle::getUV(const std::vector<float>& intersectionPoint) const {
    // Calculate vectors from v0 to intersection point and other vertices
    std::vector<float> v0p = {intersectionPoint[0] - v0[0], intersectionPoint[1] - v0[1], intersectionPoint[2] - v0[2]};
    std::vector<float> e1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    std::vector<float> e2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};

    // Calculate dot products
    float d00 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
    float d01 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
    float d11 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
    float d20 = v0p[0] * e1[0] + v0p[1] * e1[1] + v0p[2] * e1[2];
    float d21 = v0p[0] * e2[0] + v0p[1] * e2[1] + v0p[2] * e2[2];

    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    return {u, v};
}
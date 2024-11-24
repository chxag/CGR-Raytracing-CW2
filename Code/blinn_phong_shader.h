#ifndef BLINN_PHONG_SHADER_H
#define BLINN_PHONG_SHADER_H

#include <vector>
#include "material.h"
#include "ray.h"
#include "sphere.h"
#include "cylinder.h"
#include "triangle.h"
#include "light.h"
#include "shader_result.h"

class BlinnPhongShader
{
public:
    static std::vector<float> calculateColor(const std::vector<float> &intersectionPoint, 
                                             const std::vector<float> &normal, 
                                             const std::vector<float> &viewDir, 
                                             const Material &material, 
                                             const std::vector<Light> &lights, 
                                             const std::vector<Sphere> &spheres, 
                                             const std::vector<Cylinder> &cylinders, 
                                             const std::vector<Triangle> &triangles, 
                                             const std::vector<float> uv_coordinates);
    static ShaderResult intersectionTests(const Ray &ray, 
                                          const std::vector<Sphere> &spheres, 
                                          const std::vector<Cylinder> &cylinders, 
                                          const std::vector<Triangle> &triangles, 
                                          std::vector<float> &backgroundcolor);
};

#endif
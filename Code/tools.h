#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <vector>
#include "sphere.h"
#include "cylinder.h"
#include "triangle.h"
#include "ppmWriter.h"
#include "light.h"

class Tools

{
public:
    void readConfig(const std::string &filename);
    void render(PPMWriter& ppmwriter, std::string rendermode);
    std::vector<float> traceRay(const Ray& ray, int depth, const std::string& rendermode);
    std::vector<float> handleReflection(const Ray &ray, const std::vector<float> &intersectionPoint, const std::vector<float> &normal, int depth, const std::string &rendermode);
    std::vector<float> handleRefraction(const Ray &ray, const std::vector<float> &intersectionPoint, std::vector<float> &normal, const Material &material, float cos_theta, int depth, const std::string &rendermode);
    std::vector<float> combineColors(const std::vector<float>& phongColor, const std::vector<float>& reflectionColor, const std::vector<float>& refractionColor, const Material& material, const float effectiveReflectivity, float transparency);

private:

    int nbounces;
    std::string rendermode;

    std::string camera_type;
    int width;
    int height;
    std::vector<float>  position;
    std::vector<float>  lookAt;
    std::vector<float>  upVector;
    float fov;
    float exposure;

    std::vector<float>  backgroundcolor;
    std::vector<Sphere> spheres;
    std::vector<Cylinder> cylinders;
    std::vector<Triangle> triangles;
    std::vector<Light> lightsources;

    float max_value = 0.0f;


};

#endif
#include <iostream>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <vector>
#include "tools.h"
#include <stdexcept>
#include <limits>
#include <iostream>
#include <algorithm>
#include "material.h"
#include "binary_shader.h"
#include "blinn_phong_shader.h"
#include "vector_utils.h"
#include "tone_mapping.h"
#include "shadow.h"
#include <chrono>
#include <random>

using json = nlohmann::json;

float pi = 3.14159265358979323846;

void Tools::readConfig(const std::string &filename)
{
    std::ifstream file(filename);
    json j;
    file >> j;

    if (j.contains("nbounces") && j["nbounces"].is_number())
    {
        nbounces = j["nbounces"];
    }
    if (j.contains("rendermode") && j["rendermode"].is_string())
    {
        rendermode = j["rendermode"];
    }


    camera_type = j["camera"]["type"];
    width = j["camera"]["width"].get<int>();
    height = j["camera"]["height"].get<int>();
    position = j["camera"]["position"].get<std::vector<float>>();
    lookAt = j["camera"]["lookAt"].get<std::vector<float>>();
    upVector = j["camera"]["upVector"].get<std::vector<float>>();
    fov = j["camera"]["fov"].get<float>();
    exposure = j["camera"]["exposure"].get<float>();

    if (j.contains("lens_sampling") && j["lens_sampling"].is_boolean()){
        lens_sampling = j["lens_sampling"].get<bool>();
    }

    backgroundcolor = j["scene"]["backgroundcolor"].get<std::vector<float>>();

    if (j["scene"].contains("lightsources"))
    {
        for (const auto &light : j["scene"]["lightsources"])
        {
            std::string light_type = light["type"].get<std::string>();
            std::vector<float> light_position = light["position"].get<std::vector<float>>();
            std::vector<float> intensity = light["intensity"].get<std::vector<float>>();

            lightsources.emplace_back(light_type, light_position, intensity);
        }
    }

    for (const auto &shape : j["scene"]["shapes"])
    {
        if (shape["type"].get<std::string>() == "sphere")
        {
            std::vector<float> center = {shape["center"][0].get<float>(), shape["center"][1].get<float>(), shape["center"][2].get<float>()};
            float radius = shape["radius"].get<float>();

            if (rendermode == "phong" && shape.contains("material"))
            {
                material = readMaterial(shape["material"]);
            } 
            spheres.emplace_back(center, radius, material);

        }

        if (shape["type"].get<std::string>() == "cylinder")
        {
            std::vector<float> center = {shape["center"][0].get<float>(), shape["center"][1].get<float>(), shape["center"][2].get<float>()};
            float radius = shape["radius"].get<float>();
            std::vector<float> axis = {shape["axis"][0].get<float>(), shape["axis"][1].get<float>(), shape["axis"][2].get<float>()};
            float height = shape["height"].get<float>();

            if (rendermode == "phong" && shape.contains("material"))
            {
                
                material = readMaterial(shape["material"]);
            } 
            
            cylinders.emplace_back(center, radius, axis, height, material);
        }

        if (shape["type"].get<std::string>() == "triangle")
        {
            std::vector<float> v0 = {shape["v0"][0].get<float>(), shape["v0"][1].get<float>(), shape["v0"][2].get<float>()};
            std::vector<float> v1 = {shape["v1"][0].get<float>(), shape["v1"][1].get<float>(), shape["v1"][2].get<float>()};
            std::vector<float> v2 = {shape["v2"][0].get<float>(), shape["v2"][1].get<float>(), shape["v2"][2].get<float>()};
        
            if (rendermode == "phong" && shape.contains("material"))
            {
                material = readMaterial(shape["material"]); 
            } 

            triangles.emplace_back(v0, v1, v2, material);
        }
    }

    if (j.contains("bvh") && j["bvh"].is_boolean())
    {
        useBVH = j["bvh"].get<bool>();
        if (useBVH && !spheres.empty() && !cylinders.empty() && !triangles.empty()){
            try {
                std::cout << "Building BVH..." << std::endl;
                bvh = std::make_unique<BVH>(spheres, cylinders, triangles);
                std::cout << "BVH built successfully!" << std::endl;
            } catch (const std::exception &e) {
                std::cerr << "Error building BVH: " << e.what() << std::endl;
                useBVH = false;
            }
        }
    }
};

Material Tools::readMaterial(const json &material_json){
    try{
        float ks_coeffcient = material_json.value("ks", 0.0f);
        float kd_coeffcient = material_json.value("kd", 0.0f);
        float specular_exponent = material_json.value("specularexponent", 0.0f);
        std::vector<float> specular_color = material_json.value("specularcolor", std::vector<float>{0.0f, 0.0f, 0.0f});
        std::vector<float> diffuse_color = material_json.value("diffusecolor", std::vector<float>{0.0f, 0.0f, 0.0f});
        bool is_reflective = material_json.value("isreflective", false);
        float reflectivity = material_json.value("reflectivity", 0.0f);
        bool is_refractive = material_json.value("isrefractive", false);
        float refractive_index = material_json.value("refractiveindex", 0.0f);

        Texture *texture = nullptr;
        if (material_json.contains("texture") && !material_json["texture"].is_null())
        {
            std::string texture_path = material_json["texture"].get<std::string>();
            texture = new Texture(texture_path);
        }

         return Material(ks_coeffcient, kd_coeffcient, specular_exponent, diffuse_color, specular_color, is_reflective, reflectivity, is_refractive, refractive_index, texture);
    } catch (const std::exception &e){
        std::cerr << "Error reading material: " << e.what() << std::endl;
        return Material();
    }
}

std::vector<float> Tools::handleReflection(const Ray &ray, const std::vector<float> &intersectionPoint, const std::vector<float> &normal, int depth, const std::string &rendermode)
{
    std::vector<float> reflectionDir = reflect(ray.direction, normal);
    normalize(reflectionDir);
    std::vector<float> temp = {intersectionPoint[0] + 0.001f * reflectionDir[0],
                               intersectionPoint[1] + 0.001f * reflectionDir[1],
                               intersectionPoint[2] + 0.001f * reflectionDir[2]};
    Ray reflectionRay(temp, reflectionDir);
    return traceRay(reflectionRay, depth + 1, rendermode);
};

std::vector<float> Tools::handleRefraction(const Ray &ray, const std::vector<float> &intersectionPoint, std::vector<float> &normal, const Material &material, float cos_theta, int depth, const std::string &rendermode)
{
    float eta_ratio = material.refractive_index;
    if (cos_theta < 0.0f)
    {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
        eta_ratio = 1.0f / eta_ratio;
    }

    std::vector<float> refractedDir = refract(ray.direction, normal, eta_ratio);
    if (!refractedDir.empty())
    {
        normalize(refractedDir);
        std::vector<float> refractionPoint{
            intersectionPoint[0] + 0.001f * refractedDir[0],
            intersectionPoint[1] + 0.001f * refractedDir[1],
            intersectionPoint[2] + 0.001f * refractedDir[2]};
        Ray refractionRay(refractionPoint, refractedDir);
        return traceRay(refractionRay, depth + 1, rendermode);
    }
    return {0.0f, 0.0f, 0.0f};
};

std::vector<float> Tools::combineColors(const std::vector<float> &phongColor, const std::vector<float> &reflectionColor, const std::vector<float> &refractionColor, const Material &material, float effectiveReflectivity, float transparency)
{
    std::vector<float> finalColor(3);

    float reflectivity = material.is_reflective ? material.reflectivity : 0.0f;

    float refractionFactor = material.is_refractive ? (1.0f - reflectivity) : 0.0f;

    float phongFactor = 1.0f - effectiveReflectivity - transparency;
    if (phongFactor < 0.0f)
    {
        phongFactor = 0.0f;
    }

    for (int i = 0; i < 3; ++i)
    {
        finalColor[i] = phongColor[i] * (1.0f - reflectivity - refractionFactor) + reflectionColor[i] * reflectivity + refractionColor[i] * refractionFactor;

        finalColor[i] = std::min(std::max(finalColor[i], 0.0f), 1.0f);
    };

    return finalColor;
}

std::vector<float> Tools::traceRay(const Ray &ray, int depth, const std::string &rendermode)
{

    if (depth > nbounces)
    {
        return backgroundcolor;
    }

    std::vector intersection_color = backgroundcolor;

    if (rendermode == "phong")
    {
        ShaderResult result = BlinnPhongShader::intersectionTests(ray, spheres, cylinders, triangles, backgroundcolor, useBVH ? bvh.get() : nullptr);
        std::vector<float> uv_coordinates = result.uv_coordinates;
        intersection_color = result.color;
        bool intersected = result.intersected;
        std::vector<float> intersectionPoint = result.intersection_point;
        Material intersectedMaterial = result.intersected_material;
        std::vector<float> normal = result.normal;

        if (intersected)
        {
            std::vector<float> viewDir = {
                position[0] - intersectionPoint[0],
                position[1] - intersectionPoint[1],
                position[2] - intersectionPoint[2]};
            normalize(viewDir);
            float cos_theta = -(ray.direction[0] * normal[0] +
                                ray.direction[1] * normal[1] +
                                ray.direction[2] * normal[2]);

            std::vector<float> phong_color = BlinnPhongShader::calculateColor(intersectionPoint, normal, viewDir, intersectedMaterial, lightsources, spheres, cylinders, triangles, uv_coordinates, useBVH ? bvh.get() : nullptr);

            std::vector<float> reflectionColor = {0.0f, 0.0f, 0.0f};
            std::vector<float> refractionColor = {0.0f, 0.0f, 0.0f};

            if (intersectedMaterial.is_reflective)
            {
                reflectionColor = handleReflection(ray, intersectionPoint, normal, depth, rendermode);
            }
            if (intersectedMaterial.is_refractive)
            {
                refractionColor = handleRefraction(ray, intersectionPoint, normal, intersectedMaterial, cos_theta, depth, rendermode);
            }
            float reflectivity = intersectedMaterial.is_reflective ? intersectedMaterial.reflectivity : 0.0f;
            float transparency = intersectedMaterial.is_refractive ? (1.0f - reflectivity) : 0.0f;
            float fresnel = pow(1.0f - fabs(cos_theta), 5.0f);
            float effectiveReflectivity = reflectivity * fresnel;

            intersection_color = combineColors(phong_color, reflectionColor, refractionColor, intersectedMaterial, effectiveReflectivity, transparency);
        }
    }

    if (rendermode == "binary")
    {
        ShaderResult result = BinaryShader::calculateColor(ray, spheres, cylinders, triangles, backgroundcolor);
        intersection_color = result.color;
    }

    return intersection_color;
}

void Tools::render(PPMWriter &ppmwriter, std::string rendermode)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<float> forward = {lookAt[0] - position[0], lookAt[1] - position[1], lookAt[2] - position[2]};
    normalize(forward);
    std::vector<float> right = {upVector[1] * forward[2] - upVector[2] * forward[1],
                                upVector[2] * forward[0] - upVector[0] * forward[2],
                                upVector[0] * forward[1] - upVector[1] * forward[0]};
    normalize(right);
    std::vector<float> up = {forward[1] * right[2] - forward[2] * right[1],
                             forward[2] * right[0] - forward[0] * right[2],
                             forward[0] * right[1] - forward[1] * right[0]};

    normalize(up);
    float aspectRatio = static_cast<float>(width) / height;
    float scale = tan(fov * 0.5 * pi / 180.0f);

    const int samples = 16;
    const int gridSize = static_cast<int>(sqrt(samples));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);


    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            std::vector<float> pixel_color = {0.0f, 0.0f, 0.0f};

            if (lens_sampling){
                aperture = 0.01f;
                focalDistance = 5.0f;
                
                for(int py = 0; py < gridSize; py++){
                    for(int px = 0; px < gridSize; px++){

                        float cellX = px + 0.5f;
                        float cellY = py + 0.5f;
                        
                        float offsetX = (dis(gen) - 0.5f) * 0.9f;
                        float offsetY = (dis(gen) - 0.5f) * 0.9f;

                        float sampleX = (cellX + offsetX) / gridSize;
                        float sampleY = (cellY + offsetY) / gridSize;

                        float u = (2 * (x + sampleX) / width - 1) * aspectRatio * scale;
                        float v = (1 - 2 * (y + sampleY) / height) * scale;

                        std::vector<float> focusPoint = {
                            position[0] + (right[0] * u + up[0] * v + forward[0]) * focalDistance,
                            position[1] + (right[1] * u + up[1] * v + forward[1]) * focalDistance,
                            position[2] + (right[2] * u + up[2] * v + forward[2]) * focalDistance
                        };
                        
                        float r = aperture * sqrt(dis(gen));
                        float theta = 2 * M_PI * dis(gen);

                        std::vector<float> rayOrigin = {
                            position[0] + r * cos(theta) * right[0] + r * sin(theta) * up[0],
                            position[1] + r * cos(theta) * right[1] + r * sin(theta) * up[1],
                            position[2] + r * cos(theta) * right[2] + r * sin(theta) * up[2]
                        };

                        std::vector<float> direction = {
                            focusPoint[0] - rayOrigin[0],
                            focusPoint[1] - rayOrigin[1],
                            focusPoint[2] - rayOrigin[2]
                        };
                        normalize(direction);

                        Ray ray(rayOrigin, direction);
                        std::vector<float> intersection_color_sample = traceRay(ray, 0, rendermode);

                        pixel_color[0] += intersection_color_sample[0];
                        pixel_color[1] += intersection_color_sample[1];
                        pixel_color[2] += intersection_color_sample[2];
                    }
                }
            } else {

                for (int py = 0; py < gridSize; py++){
                    for (int px = 0; px < gridSize; px++){
                        float cellX = px + 0.5f;
                        float cellY = py + 0.5f;

                        float offsetX = (dis(gen) - 0.5f) * 0.9f;
                        float offsetY = (dis(gen) - 0.5f) * 0.9f;

                        float sampleX = (cellX + offsetX) / gridSize;
                        float sampleY = (cellY + offsetY) / gridSize;

                        float u = (2 * (x + sampleX)/ width - 1) * aspectRatio * scale;
                        float v = (1 - 2 * (y + sampleY) / height) * scale;

                        std::vector<float> direction = {right[0] * u + up[0] * v + forward[0],
                                                        right[1] * u + up[1] * v + forward[1],
                                                        right[2] * u + up[2] * v + forward[2]};
                        normalize(direction);
                        Ray ray(position, direction);

                        std::vector<float> intersection_color_sample = traceRay(ray, 0, rendermode);

                        pixel_color[0] += intersection_color_sample[0];
                        pixel_color[1] += intersection_color_sample[1];
                        pixel_color[2] += intersection_color_sample[2];
                    }
                }
            }

            float inv_samples = 1.0f / (gridSize * gridSize);
            pixel_color[0] *= inv_samples;
            pixel_color[1] *= inv_samples;
            pixel_color[2] *= inv_samples;

            pixel_color = linearToneMapping(pixel_color);

            ppmwriter.getPixelData(x, y, {static_cast<unsigned char>(pixel_color[0] * 255), static_cast<unsigned char>(pixel_color[1] * 255), static_cast<unsigned char>(pixel_color[2] * 255)});
        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Render Time:" << duration.count()/ 1000.0 << "seconds" << std::endl;
}
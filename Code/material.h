#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include "texture.h"
#include <memory>

struct Material
{
    float ks_coeffcient;
    float kd_coeffcient;
    float specular_exponent;
    std::vector<float> diffuse_color;
    std::vector<float> specular_color;
    bool is_reflective;
    float reflectivity;
    bool is_refractive;
    float refractive_index;
    Texture *texture = nullptr;

    ~Material()
    {
        if (texture){
            delete texture;
        }
    }

    Material() : ks_coeffcient(0.0f), kd_coeffcient(0.0f), specular_exponent(0.0f), diffuse_color({0.0f, 0.0f, 0.0f}), specular_color({0.0f, 0.0f, 0.0f}), is_reflective(false), reflectivity(0.0f), is_refractive(false), refractive_index(0.0f), texture{nullptr} {}

    Material(float ks_coeffcient, float kd_coeffcient, float specular_exponent, std::vector<float> diffuse_color, std::vector<float> specular_color, bool is_reflective, float reflectivity, bool is_refractive, float refractive_index, Texture* texture)
        : ks_coeffcient(ks_coeffcient), kd_coeffcient(kd_coeffcient), specular_exponent(specular_exponent), diffuse_color(diffuse_color), specular_color(specular_color), is_reflective(is_reflective), reflectivity(reflectivity), is_refractive(is_refractive), refractive_index(refractive_index), texture(texture) {}

    Material(const Material &other_material) : ks_coeffcient(other_material.ks_coeffcient), kd_coeffcient(other_material.kd_coeffcient), specular_exponent(other_material.specular_exponent), diffuse_color(other_material.diffuse_color), specular_color(other_material.specular_color), is_reflective(other_material.is_reflective), reflectivity(other_material.reflectivity), is_refractive(other_material.is_refractive), refractive_index(other_material.refractive_index)
    {
        if (other_material.texture)
        {
            texture = new Texture (*other_material.texture);
        } else {
            texture = nullptr;
        }
    }

    Material& operator=(const Material& other_material) {
        if (this == &other_material) {
            return *this;
        }

        ks_coeffcient = other_material.ks_coeffcient;
        kd_coeffcient = other_material.kd_coeffcient;
        specular_exponent = other_material.specular_exponent;
        diffuse_color = other_material.diffuse_color;
        specular_color = other_material.specular_color;
        is_reflective = other_material.is_reflective;
        reflectivity = other_material.reflectivity;
        is_refractive = other_material.is_refractive;
        refractive_index = other_material.refractive_index;

        delete texture;
        if (other_material.texture)
        {
            texture = new Texture (*other_material.texture);
        } else {
            texture = nullptr;
        }
        return *this;
    }
};


#endif
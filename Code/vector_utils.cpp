#include "vector_utils.h"

void normalize(std::vector<float> &vec)
{
    float length = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (length > 1e-6)
    {
        for (float &v : vec)
        {
            v = v / length;
        }
    }
}

std::vector<float> reflect(const std::vector<float> &incident, const std::vector<float> &normal)
{
    float dot = incident[0] * normal[0] + incident[1] * normal[1] + incident[2] * normal[2];
    std::vector<float> reflected = {
        incident[0] - 2 * dot * normal[0],
        incident[1] - 2 * dot * normal[1],
        incident[2] - 2 * dot * normal[2]};
    return reflected;
}

std::vector<float> refract(const std::vector<float> &incident, const std::vector<float> &normal, float eta_ratio)
{
    float cos_theta = -(incident[0] * normal[0] + incident[1] * normal[1] + incident[2] * normal[2]);
    float k = eta_ratio * eta_ratio * (1.0f - cos_theta * cos_theta);
    if (k < 0.0f) return {};
    float sqrt_k = sqrt(k);
    return {
        eta_ratio * incident[0] + (eta_ratio * cos_theta - sqrt_k) * normal[0],
        eta_ratio * incident[1] + (eta_ratio * cos_theta - sqrt_k) * normal[1],
        eta_ratio * incident[2] + (eta_ratio * cos_theta - sqrt_k) * normal[2]};
}
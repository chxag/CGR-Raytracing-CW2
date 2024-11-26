#include "tone_mapping.h"
#include <algorithm>

std::vector<float> linearToneMapping(const std::vector<float> &color)
{
    std::vector<float> tone_mapped_color(3);
    float exposure = 1.0f;
    
    // Scale each color channel by the max_value and clamp it between [0, 1]
    for (size_t i = 0; i < 3; ++i)
    {

        float exposed = color[i] * exposure;
        tone_mapped_color[i] = exposed / (1.0f + exposed);

        tone_mapped_color[i] = std::min(std::max(tone_mapped_color[i], 0.0f), 1.0f);
    }

    return tone_mapped_color;
}

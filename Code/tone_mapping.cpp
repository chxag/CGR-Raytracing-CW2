#include "tone_mapping.h"
#include <algorithm>

std::vector<float> linearToneMapping(const std::vector<float> &color)
{
    std::vector<float> tone_mapped_color(3);
    // float luminance = 0.2126f * color[0] + 0.7152f * color[1] + 0.0722f * color[2];

    const float boost = 1.2f;
    for (size_t i = 0; i < 3; ++i){
        tone_mapped_color[i] = std::min(std::max(color[i] * boost, 0.0f), 1.0f);
    }

    return tone_mapped_color;
}

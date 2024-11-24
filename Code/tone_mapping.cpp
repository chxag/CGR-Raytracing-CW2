#include "tone_mapping.h"
#include <algorithm>

std::vector<float> linearToneMapping(const std::vector<float> &color, float max_value)
{
    std::vector<float> tone_mapped_color(3);
    
    // Scale each color channel by the max_value and clamp it between [0, 1]
    for (size_t i = 0; i < 3; ++i)
    {
        tone_mapped_color[i] = std::min(color[i] / max_value, 1.0f);
    }

    return tone_mapped_color;
}

#ifndef TEXTURE_H
#define TEXTURE_H

#include <string>
#include <vector>

class Texture
{
public:
    Texture(const std::string &filename);
    std::vector<unsigned char> getTexel(float u, float v) const;

private:
    int width, height;
    std::vector<unsigned char> pixelData;
    void loadTexture(const std::string &filename);
};

#endif
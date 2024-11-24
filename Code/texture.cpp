#include "texture.h"
#include <fstream>
#include <stdexcept>
#include <sstream>

Texture::Texture(const std::string &filename)
{
    loadTexture(filename);
}

void Texture::loadTexture(const std::string &filename)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        throw std::runtime_error("Could not open texture file: " + filename);
    }

    std::string line;
    std::getline(file, line);
    if (line != "P6")
    {
        throw std::runtime_error("Unsupported file format (only P6 PPM is supported)");
    }

    do
    {
        std::getline(file, line);
    } while (line[0] == '#');

    std::istringstream ss(line);
    ss >> width >> height;

    int maxColorValue;
    file >> maxColorValue;
    file.get();

    pixelData.resize(height * width * 3);

    file.read(reinterpret_cast<char *>(pixelData.data()), pixelData.size());
    if (file.gcount() != static_cast<std::streamsize>(pixelData.size()))
    {
        throw std::runtime_error("Error reading texture file: " + filename);
    }
    file.close();
}

std::vector<unsigned char> Texture::getTexel(float u, float v) const
{
    u = std::max(0.0f, std::min(1.0f, u));
    v = std::max(0.0f, std::min(1.0f, v));

    int x = static_cast<int>(u * (width - 1));
    int y = static_cast<int>(v * (height - 1));
    int index = (y * width + x) * 3;

    if (index < 2 || index + 2 >= static_cast<int>(pixelData.size())) 
    {
        return {0, 0, 0};
    }

    return {pixelData[index], pixelData[index + 1], pixelData[index + 2]};
}

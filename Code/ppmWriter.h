#ifndef PPM_WRITER_H
#define PPM_WRITER_H

#include <string>
#include <vector>

class PPMWriter
{
    public:
        PPMWriter(int width, int height, const std::vector<unsigned char>& backgrounddata);
        void getPixelData(int x, int y, const std::vector<unsigned char>& colordata);
        void writePPM(const std::string& filename) const;
    
    private:
        int width;
        int height;
        std::vector<unsigned char> backgrounddata;
        std::vector<unsigned char> pixeldata;
};

#endif
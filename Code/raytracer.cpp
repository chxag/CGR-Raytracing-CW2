#include "tools.h"
#include "ppmWriter.h"
#include <iostream>

int main()
{
    // Initialize the tools
    Tools tools;
    // Read the configuration file
    tools.readConfig("Jsons/scene.json");
    // Set the width and height of the image
    int width = 1200;
    int height = 800;
    // Set the background color
    std::vector<unsigned char> backgrounddata = {64, 64, 64};
    // Initialize the PPM writer
    PPMWriter ppmwriter(width, height, backgrounddata);
    // Render the scene
    tools.render(ppmwriter, "phong");
    // Write the PPM file
    ppmwriter.writePPM("output.ppm");
    return 0;
}
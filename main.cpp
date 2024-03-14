/*
Lab 2 - Ray Tracer
Victor Fenton Aguilar
CS 34300 Advanced Computer Graphics
Purdue University
*/

// GLM
#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_transform.hpp"

// Image output
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Default includes
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cctype>

using namespace std;

typedef struct
{
    glm::vec3 origin;
    glm::vec3 direction;
} Ray;

class Intersection
{
public:
    float t;
    glm::vec3 hitLocation;
    glm::vec3 normal;
};

class Light
{
public:
    glm::vec3 pos;
    glm::vec3 diffuse;
    glm::vec3 specular;
};

// Base Class
class Surface
{
public:
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;

    virtual bool intersect(Intersection res, Ray r) = 0;
};

class Sphere : public Surface
{
public:
    glm::vec3 pos;
    float radius;

    bool intersect(Intersection res, Ray r)
    {
        return true;
    }
};

pair<Sphere, float> intersect(Ray r, float tMin, float tMax)
{
    Sphere s = Sphere();
    pair<Sphere, float> pair1;
    pair1 = make_pair(s, 0.0);
    return pair1;
}

glm::vec3 TraceRay(Ray r, int depth)
{
    return glm::vec3(0, 0, 0);
}

Ray CalculateRay()
{
    Ray r;
    return r;
}

inline void ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
                                    { return !std::isspace(ch); }));
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "Error: specify an input file for ray tracing." << endl;
        return 1;
    }

    // Set up scene, define and save all objects in scene
    int resX;
    int resY;
    int MAXDEPTH;
    int ANTIALIAS;
    glm::vec3 bg_color;
    const int channels = 3; // Red, Green, Blue
    vector<Surface> objects;
    vector<Light> lights;

    // Parse Input File from command line arg filename
    string line;
    ifstream file(argv[1]);
    if (file.is_open())
    {
        while (getline(file, line))
        {
            ltrim(line);

            if (line[0] == '/' && line[1] == '/')
                continue;

            if (line.rfind("ANTIALIAS ", 0) == 0)
            {
                sscanf(line.substr(10).c_str(), "%d", &ANTIALIAS);
                continue;
            }

            if (line.rfind("BACKGROUND ", 0) == 0)
            {
                int r, g, b;
                sscanf(line.substr(11).c_str(), "%d %d %d", &r, &g, &b);
                bg_color = glm::vec3(r, g, b);
                continue;
            }

            if (line.rfind("MAXDEPTH ", 0) == 0)
            {
                sscanf(line.substr(9).c_str(), "%d", &MAXDEPTH);
                continue;
            }

            if (line.rfind("RESOLUTION ", 0) == 0)
            {
                sscanf(line.substr(11).c_str(), "%d %d", &resX, &resY);
                continue;
            }

            if (line.rfind("LIGHT", 0) == 0)
            {
                string format[3] = {"POS ", "DIFF ", "SPEC "};
                for (int i = 0; i < 3; i++)
                {
                    if (!getline(file, line))
                    {
                        cout << "Invalid input file format" << endl;
                        return 2;
                    }

                    ltrim(line);

                    if (line.rfind(format[i], 0) != 0)
                    {
                        cout << "Invalid Light Format in input file" << endl;
                        return 2;
                    }

                    float x, y, z;
                    sscanf(line.substr(format[i].length()).c_str(), "%f %f %f", &x, &y, &z);

                    Light l;
                    if (i == 0)
                        l.pos = glm::vec3(x, y, z);
                    else if (i == 1)
                        l.diffuse = glm::vec3(x, y, z);
                    else
                        l.specular = glm::vec3(x, y, z);

                    lights.push_back(l);
                }
            }
        }
    }

    cout << lights << endl;

    std::vector<unsigned char> image_data(resX * resY * channels, 0);

    for (int i = 0; i < resX; ++i)
        for (int j = 0; j < resY; ++j)
        {
            // Ray ray = CalculateRay(); // Get the primary ray
            // glm::vec3 color = TraceRay(ray, MAXDEPTH);
            glm::vec3 color = glm::vec3(255, 0, 0);
            int index = (j * resX + i) * channels;
            image_data[index] = color[0];
            image_data[index + 1] = color[1];
            image_data[index + 2] = color[2];
        }

    // Output Image
    stbi_write_png("ray_output.png", resX, resY, channels, image_data.data(), resX * channels);

    return 0;
}

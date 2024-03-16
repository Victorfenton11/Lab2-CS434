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

    virtual float intersect(Intersection *res, Ray r) = 0;
};

class Sphere : public Surface
{
public:
    glm::vec3 pos;
    float radius;

    float intersect(Intersection *res, Ray r)
    {

        return true;
    }
};

class Quad : public Surface
{
public:
    glm::vec3 v1;
    glm::vec3 v2;
    glm::vec3 v3;

    float intersect(Intersection *res, Ray r)
    {
        return true;
    }
};

Ray Reflected(Ray r, glm::vec3 p)
{
    Ray r;
    return r;
}

Ray Refracted(Ray r, glm::vec3 p)
{
    Ray r;
    return r;
}

glm::vec3 Phong(glm::vec3 p, Light l)
{
    return glm::vec3(0.0f, 0.0f, 0.0f);
}

int CastShadowRays(glm::vec3 p)
{
    return 0;
}

bool FirstIntersection(Ray r, glm::vec3 *p, vector<Sphere> spheres, vector<Quad> quads)
{
    Intersection first;
    first.t = 1000000.0f;
    bool hit = false;
    // Check Spheres
    for (int i = 0; i < spheres.size(); i++)
    {
        Intersection current;
        if (spheres[i].intersect(&current, r))
        {
            hit = true;
            if (current.t < first.t)
                first = current;
        }
    }

    // Check quads
    for (int j = 0; j < quads.size(); j++)
    {
        Intersection current;
        if (quads[j].intersect(&current, r))
        {
            hit = true;
            if (current.t < first.t)
                first = current;
        }
    }

    if (!hit)
        return false;

    *p = first.hitLocation;
    return true;
}

void TraceRay(Ray r, int depth, glm::vec3 *color, vector<Light> lights, vector<Sphere> spheres, vector<Quad> quads)
{
    glm::vec3 intersection;
    if (!FirstIntersection(r, &intersection, spheres, quads))
        return;
    int contributedLights = CastShadowRays(intersection);
    for (int i = 0; i < contributedLights; i++)
        *color += Phong(intersection, lights[i]);
    if (depth <= 0)
        return;
    Ray reflected = Reflected(r, intersection);
    Ray refracted = Refracted(r, intersection);
    TraceRay(reflected, depth - 1, color, lights, spheres, quads);
    TraceRay(refracted, depth - 1, color, lights, spheres, quads);
    return;
}

Ray CalculateRay(glm::vec3 cameraPos, glm::vec3 up, int i, int j, int resX, int resY, float fov)
{
    glm::vec3 lookAt = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 l = glm::normalize(lookAt - cameraPos);
    glm::vec3 v = glm::normalize(glm::cross(l, up));
    glm::vec3 u = glm::cross(v, l);

    float fov_rad = glm::radians(fov);
    float aspect_ratio = (float)resX / (float)resY;
    glm::vec3 ll = cameraPos + (1 / glm::tan(fov_rad / 2)) * l - aspect_ratio * v - u;
    glm::vec3 p = ll + ((2.0f * aspect_ratio * v * (float)i) / (float)resX) + (2.0f * u * (float)j);
    glm::vec3 d = glm::normalize(p - cameraPos);

    Ray r;
    r.origin = cameraPos;
    r.direction = d;

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
    float fov = 60.0;
    vector<Light> lights;
    vector<Sphere> spheres;
    vector<Quad> quads;

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

            if (line.rfind("FOV ", 0) == 0)
            {
                sscanf(line.substr(4).c_str(), "%f", &fov);
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
                Light l;
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

                    if (i == 0)
                        l.pos = glm::vec3(x, y, z);
                    else if (i == 1)
                        l.diffuse = glm::vec3(x, y, z);
                    else
                        l.specular = glm::vec3(x, y, z);
                }

                lights.push_back(l);
            }

            if (line.rfind("SPHERE", 0) == 0)
            {
                string format[4] = {"POS ", "DIFF ", "SPEC ", "SHININESS "};
                Sphere s;
                for (int i = 0; i < 4; i++)
                {
                    if (!getline(file, line))
                    {
                        cout << "Invalid input file format" << endl;
                        return 2;
                    }

                    ltrim(line);

                    if (line.rfind(format[i], 0) != 0)
                    {
                        cout << "Invalid Sphere Format in input file" << endl;
                        return 2;
                    }

                    if (i == 3)
                    {
                        float shininess;
                        sscanf(line.substr(format[i].length()).c_str(), "%f", &shininess);
                        s.shininess = shininess;
                        break;
                    }

                    float x, y, z;
                    sscanf(line.substr(format[i].length()).c_str(), "%f %f %f", &x, &y, &z);

                    if (i == 0)
                    {
                        s.pos = glm::vec3(x, y, z);
                        int rad_index = line.find("RADIUS");
                        if (rad_index == string::npos)
                        {
                            cout << "Invalid Sphere format in input file" << endl;
                            return 2;
                        }
                        sscanf(line.substr(rad_index + 6).c_str(), "%f", &s.radius);
                    }
                    else if (i == 1)
                        s.diffuse = glm::vec3(x, y, z);
                    else
                        s.specular = glm::vec3(x, y, z);
                }

                spheres.push_back(s);
            }

            if (line.rfind("QUAD", 0) == 0)
            {
                string format[6] = {"POS ", "POS ", "POS ", "DIFF ", "SPEC ", "SHININESS "};
                Quad q;
                for (int i = 0; i < 6; i++)
                {
                    if (!getline(file, line))
                    {
                        cout << "Invalid input file format" << endl;
                        return 2;
                    }

                    ltrim(line);

                    if (line.rfind(format[i], 0) != 0)
                    {
                        cout << "Invalid Quad Format in input file" << endl;
                        return 2;
                    }

                    if (i == 5)
                    {
                        float shininess;
                        sscanf(line.substr(format[i].length()).c_str(), "%f", &shininess);
                        q.shininess = shininess;
                        break;
                    }

                    float x, y, z;
                    sscanf(line.substr(format[i].length()).c_str(), "%f %f %f", &x, &y, &z);

                    if (i == 0)
                        q.v1 = glm::vec3(x, y, z);
                    else if (i == 1)
                        q.v2 = glm::vec3(x, y, z);
                    else if (i == 2)
                        q.v3 = glm::vec3(x, y, z);
                    else if (i == 3)
                        q.diffuse = glm::vec3(x, y, z);
                    else
                        q.specular = glm::vec3(x, y, z);
                }

                quads.push_back(q);
            }
        }
    }

    std::vector<unsigned char> image_data(resX * resY * channels, 0);

    glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 0.0f);
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    for (int i = 0; i < resX; ++i)
        for (int j = 0; j < resY; ++j)
        {
            Ray ray = CalculateRay(cameraPos, up, i, j, resX, resY, fov); // Get the primary ray
            glm::vec3 color;
            TraceRay(ray, MAXDEPTH, &color);
            int index = (j * resX + i) * channels;
            image_data[index] = color[0];
            image_data[index + 1] = color[1];
            image_data[index + 2] = color[2];
        }

    // Output Image
    stbi_write_png("ray_output.png", resX, resY, channels, image_data.data(), resX * channels);

    return 0;
}

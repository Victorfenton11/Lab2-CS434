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
#include <memory>
#include <omp.h>

using namespace std;

int MAXDEPTH = 1;

class Surface;
class Sphere;
class Triangle;

typedef struct
{
    glm::vec3 origin;
    glm::vec3 direction;
} Ray;

class Intersection
{
public:
    float t;
    Surface *obj;
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
    Surface() {}
    virtual ~Surface() {}
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;

    virtual bool intersect(Intersection *res, Ray r) = 0;
};

class Sphere : public Surface
{
public:
    glm::vec3 pos;
    float radius;

    bool intersect(Intersection *res, Ray r)
    {
        glm::vec3 oc = r.origin - pos;

        // Assuming ray direction will always be normalized (a = 1)
        float b = 2.0f * glm::dot(r.direction, oc);
        float c = glm::dot(oc, oc) - radius * radius;

        float discriminant = b * b - 4 * c;
        if (discriminant < 0)
            return false;

        float t1 = (-b - sqrt(discriminant)) / 2.0f;
        if (t1 > 0)
        {
            res->t = t1;
            res->hitLocation = r.origin + t1 * r.direction;
            res->normal = glm::normalize(res->hitLocation - pos);
            res->obj = this;

            // Add bias to the intersection point to avoid self-intersection
            res->hitLocation += res->normal * 0.001f;

            return true;
        }

        float t2 = (-b + sqrt(discriminant)) / 2.0f;
        if (t2 > 0)
        {
            res->t = t2;
            res->hitLocation = r.origin + t2 * r.direction;
            res->normal = glm::normalize(res->hitLocation - pos);
            res->obj = this;

            // Add bias to the intersection point
            res->hitLocation += res->normal * 0.001f;

            return true;
        }

        return false;
    }
};

class Triangle : public Surface
{
public:
    Triangle() {}
    virtual ~Triangle() {}
    glm::vec3 p1;
    glm::vec3 p2;
    glm::vec3 p3;

    bool intersect(Intersection *res, Ray r)
    {
        glm::vec3 p = p2 - p1;
        glm::vec3 q = p3 - p1;
        glm::vec3 tmp1 = glm::cross(p, q);

        float dot1 = glm::dot(tmp1, r.direction);
        float eps = 0.00000001;
        if (fabs(dot1) < eps)
            return false;

        float d = glm::dot(-tmp1, p1);
        float t = (-1.0f * (glm::dot(tmp1, r.origin) + d)) / dot1;

        if (t < 0.f)
            return false; // triangle is behind the ray

        glm::vec3 P = r.origin + (t * r.direction);

        // edge 1
        glm::vec3 vp1 = P - p1;
        glm::vec3 C = glm::cross(p, vp1);
        if (glm::dot(tmp1, C) < 0.f)
            return false;

        // edge 2
        glm::vec3 edge2 = p3 - p2;
        glm::vec3 vp2 = P - p2;
        C = glm::cross(edge2, vp2);
        if (glm::dot(tmp1, C) < 0.f)
            return false;

        // edge 3
        glm::vec3 edge3 = p1 - p3;
        glm::vec3 vp3 = P - p3;
        C = glm::cross(edge3, vp3);
        if (glm::dot(tmp1, C) < 0)
            return false;

        res->t = t;
        res->hitLocation = P;
        res->normal = glm::normalize(tmp1);       // Assuming counter-clockwise winding order
        res->hitLocation += 0.001f * res->normal; // Add bias to avoid self-intersection
        res->obj = this;
        return true;
    }
};

glm::vec3 Phong(Intersection p, Light l, glm::vec3 cameraPos)
{
    glm::vec3 diffuse = glm::vec3(0.0f);
    glm::vec3 specular = glm::vec3(0.0f);

    glm::vec3 N = p.normal;
    glm::vec3 L = glm::normalize(l.pos - p.hitLocation);
    glm::vec3 V = glm::normalize(cameraPos - p.hitLocation);
    glm::vec3 R = glm::reflect(-L, N);

    float diffuse_factor = std::max(0.0f, glm::dot(N, L));
    diffuse = p.obj->diffuse * l.diffuse * diffuse_factor;

    float specular_factor = std::pow(std::max(0.0f, glm::dot(R, V)), p.obj->shininess);
    specular = p.obj->specular * l.specular * specular_factor;

    return diffuse + specular;
}

bool FirstIntersection(Ray r, Intersection *first, vector<unique_ptr<Surface>> &objects)
{
    first->t = 100000.0f;
    bool hit = false;

    for (unsigned int i = 0; i < objects.size(); i++)
    {
        Intersection current;
        if (objects[i]->intersect(&current, r))
        {
            if (current.t < first->t)
            {
                hit = true;
                *first = current;
            }
        }
    }

    return hit;
}

vector<int> CastShadowRays(Intersection intersection, vector<Light> lights, vector<unique_ptr<Surface>> &objects)
{
    vector<int> visible_lights;
    Intersection p;
    for (unsigned int i = 0; i < lights.size(); i++)
    {
        Ray shadow_ray;
        shadow_ray.origin = intersection.hitLocation;
        shadow_ray.direction = glm::normalize(lights[i].pos - intersection.hitLocation);
        bool hit = FirstIntersection(shadow_ray, &p, objects);
        if (!hit)
        {
            visible_lights.push_back(i);
        }
        continue;
    }
    return visible_lights;
}

void TraceRay(Ray r, int depth, glm::vec3 cameraPos, glm::vec3 bg, glm::vec3 *color, vector<Light> lights, vector<unique_ptr<Surface>> &objects, float ambient, float k_spec)
{
    if (depth <= 0)
        return; // Stop if too deep recursion

    Intersection intersection;
    if (!FirstIntersection(r, &intersection, objects))
    {
        if (depth == MAXDEPTH) // If no intersection on main ray, set to background color
            *color = bg;
        return;
    }

    glm::vec3 contribution(0.f);
    vector<int> contributedLights = CastShadowRays(intersection, lights, objects);
    for (auto i : contributedLights)
        contribution += Phong(intersection, lights[i], cameraPos);

    float contribution_k = 1;
    if (depth != MAXDEPTH)
        contribution_k = k_spec;

    // Always include ambient contribution (not included in phong)
    glm::vec3 ambient_contrib = glm::vec3(ambient) * intersection.obj->diffuse;
    *color += contribution_k * 255.f * glm::clamp(ambient_contrib + contribution, 0.f, 1.f);

    float eps_cutoff = 0.00001f; // cutoff for too small contribution
    if (glm::dot(contribution, contribution) < eps_cutoff)
        return; // Stop if contribution is too small

    // Stop if we hit a diffuse surface
    if (intersection.obj->specular == glm::vec3(0.f, 0.f, 0.f))
        return;

    Ray reflected;
    reflected.origin = intersection.hitLocation;
    reflected.direction = glm::reflect(r.direction, intersection.normal);

    TraceRay(reflected, depth - 1, cameraPos, bg, color, lights, objects, ambient, contribution_k * k_spec);

    return;
}

Ray CalculateRay(glm::vec3 cameraPos, glm::vec3 up, glm::vec3 lookAt, int i, int j, int resX, int resY, float fov)
{
    float scale = glm::tan(glm::radians(fov * 0.5f));
    float aspect_ratio = (float)resX / (float)resY;

    glm::vec3 w = glm::normalize(cameraPos - lookAt);
    glm::vec3 u = glm::normalize(glm::cross(up, w));
    glm::vec3 v = glm::cross(w, u);

    // Calculate ray direction
    float x = (2 * (i + 0.5f) / static_cast<float>(resX) - 1) * aspect_ratio * scale;
    float y = (1 - 2 * (j + 0.5f) / static_cast<float>(resY)) * scale;
    glm::vec3 direction = glm::normalize(x * u + y * v - w);

    Ray ray;
    ray.origin = cameraPos;
    ray.direction = direction;

    return ray;
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
    int resX = 800;
    int resY = 800;
    int ANTIALIAS = 1;
    glm::vec3 bg_color = glm::vec3(0.f);
    const int channels = 3; // Red, Green, Blue
    float fov = 90.0f;
    float ambient = 0.0f;
    float k_spec = 0.5; // reflection coefficient (how much the contribution of a ray decreases after each reflection)
    glm::vec3 cameraPos = glm::vec3(0.f);
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 lookAt = glm::vec3(0.0f, 0.0f, -1.0f);
    vector<Light> lights;
    vector<unique_ptr<Surface>> objects;

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
                ANTIALIAS = max(1, ANTIALIAS); // Minimum factor is 1 (no antialias)
                continue;
            }

            if (line.rfind("KSPEC ", 0) == 0)
            {
                sscanf(line.substr(6).c_str(), "%f", &k_spec);
                if (k_spec < 0.f || k_spec > 1.f)
                {
                    cout << "Invalid Reflection Coefficient Entered. Tracing will proceed using KSPEC = 0.5" << endl;
                    k_spec = 0.5;
                }
                continue;
            }

            if (line.rfind("FOV ", 0) == 0)
            {
                sscanf(line.substr(4).c_str(), "%f", &fov);
                continue;
            }

            if (line.rfind("AMBIENT ", 0) == 0)
            {
                sscanf(line.substr(8).c_str(), "%f", &ambient);
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
                MAXDEPTH = max(1, MAXDEPTH); // Minimum depth is 1 (ray casting)
                continue;
            }

            if (line.rfind("RESOLUTION ", 0) == 0)
            {
                sscanf(line.substr(11).c_str(), "%d %d", &resX, &resY);
                continue;
            }

            if (line.rfind("CAMERA ", 0) == 0)
            {
                float x, y, z;
                sscanf(line.substr(7).c_str(), "%f %f %f", &x, &y, &z);
                cameraPos = glm::vec3(x, y, z);
                continue;
            }

            if (line.rfind("UP ", 0) == 0)
            {
                float x, y, z;
                sscanf(line.substr(3).c_str(), "%f %f %f", &x, &y, &z);
                up = glm::vec3(x, y, z);
                continue;
            }

            if (line.rfind("LOOKAT ", 0) == 0)
            {
                float x, y, z;
                sscanf(line.substr(7).c_str(), "%f %f %f", &x, &y, &z);
                lookAt = glm::vec3(x, y, z);
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
                Sphere *s = new Sphere();
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
                        s->shininess = shininess;
                        break;
                    }

                    float x, y, z;
                    sscanf(line.substr(format[i].length()).c_str(), "%f %f %f", &x, &y, &z);

                    if (i == 0)
                    {
                        s->pos = glm::vec3(x, y, z);
                        int rad_index = line.find("RADIUS");
                        if (rad_index == static_cast<int>(string::npos))
                        {
                            cout << "Invalid Sphere format in input file" << endl;
                            return 2;
                        }
                        sscanf(line.substr(rad_index + 6).c_str(), "%f", &s->radius);
                    }
                    else if (i == 1)
                        s->diffuse = glm::vec3(x, y, z);
                    else
                        s->specular = glm::vec3(x, y, z);
                }

                objects.push_back(unique_ptr<Sphere>(s));
            }

            if (line.rfind("TRIANGLE", 0) == 0)
            {
                string format[6] = {"POS ", "POS ", "POS ", "DIFF ", "SPEC ", "SHININESS "};
                Triangle *s = new Triangle();
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
                        cout << "Invalid Triangle Format in input file" << endl;
                        return 2;
                    }

                    if (i == 5)
                    {
                        float shininess;
                        sscanf(line.substr(format[i].length()).c_str(), "%f", &shininess);
                        s->shininess = shininess;
                        break;
                    }

                    float x, y, z;
                    sscanf(line.substr(format[i].length()).c_str(), "%f %f %f", &x, &y, &z);

                    if (i == 0)
                        s->p1 = glm::vec3(x, y, z);
                    else if (i == 1)
                        s->p2 = glm::vec3(x, y, z);
                    else if (i == 2)
                        s->p3 = glm::vec3(x, y, z);
                    else if (i == 3)
                        s->diffuse = glm::vec3(x, y, z);
                    else
                        s->specular = glm::vec3(x, y, z);
                }

                objects.push_back(unique_ptr<Triangle>(s));
            }

            if (line.rfind("QUAD", 0) == 0)
            {
                string format[6] = {"POS ", "POS ", "POS ", "DIFF ", "SPEC ", "SHININESS "};
                Triangle *s1 = new Triangle();
                Triangle *s2 = new Triangle();
                glm::vec3 p1, p2, p3;
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
                        s1->shininess = shininess;
                        s2->shininess = shininess;
                        break;
                    }

                    float x, y, z;
                    sscanf(line.substr(format[i].length()).c_str(), "%f %f %f", &x, &y, &z);

                    if (i == 0)
                        p1 = glm::vec3(x, y, z);
                    else if (i == 1)
                        p2 = glm::vec3(x, y, z);
                    else if (i == 2)
                        p3 = glm::vec3(x, y, z);
                    else if (i == 3)
                    {
                        s1->diffuse = glm::vec3(x, y, z);
                        s2->diffuse = glm::vec3(x, y, z);
                    }
                    else
                    {
                        s1->specular = glm::vec3(x, y, z);
                        s2->specular = glm::vec3(x, y, z);
                    }
                }

                // Break quad into two triangles
                // Triangle 1
                s1->p1 = p1;
                s1->p2 = p2;
                s1->p3 = p3;
                objects.push_back(unique_ptr<Triangle>(s1));

                // Triangle 2
                s2->p1 = p1;
                s2->p2 = p3;
                s2->p3 = p3 + (p1 - p2);
                objects.push_back(unique_ptr<Triangle>(s2));
            }
        }
    }

    std::vector<unsigned char> image_data(resX * resY * channels, 0);

    int supersample_resX = resX * ANTIALIAS;
    int supersample_resY = resY * ANTIALIAS;

#pragma omp parallel for schedule(dynamic) // paralellize for performance spped boost
    for (int i = 0; i < resX; ++i)
        for (int j = 0; j < resY; ++j)
        {
            // Anti-aliasing
            glm::vec3 color_sum(0.f);
            for (int m = 0; m < ANTIALIAS; m++)
                for (int n = 0; n < ANTIALIAS; n++)
                {
                    Ray ray = CalculateRay(cameraPos, up, lookAt, i * ANTIALIAS + m, j * ANTIALIAS + n, supersample_resX, supersample_resY, fov); // Get the primary ray
                    glm::vec3 color(0.f);
                    TraceRay(ray, MAXDEPTH, cameraPos, bg_color, &color, lights, objects, ambient, k_spec);
                    color_sum += color;
                }

            glm::vec3 final_color = color_sum / (float)(ANTIALIAS * ANTIALIAS);

            int index = (j * resX + i) * channels;
            image_data[index] = final_color[0];
            image_data[index + 1] = final_color[1];
            image_data[index + 2] = final_color[2];
        }

    // Output Image
    stbi_write_png("ray_output.png", resX, resY, channels, image_data.data(), resX * channels);

    return 0;
}

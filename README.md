# Lab2-CS434-RayTracing

Lab 2 - Ray Tracer
Victor Fenton Aguilar
CS 34300 Advanced Computer Graphics
Purdue University

# Running Instructions

cd into the project's directory and compile the code using g++ as follows (Note: valid installation of g++ supporting openMP required):

g++ -o ray_tracer main.cpp -fopenmp

once compiled, run the executable as follows
.\ray_tracer.exe [relative path to input file]

The output of running the ray tracer with the specified input file will be seen by the generated image 'ray_output.png'. The image will be saved to the same directory of the source code file. If there is an error in processing the input file, a corresponding message may print to the console where the program was ran, altough input file error checking was not thoroughly implemented.

# Testing

run the compiled executable specifying one of the provided files in the input_files folder, or create your own input file.

To run provided input files:
.\ray_tracer.exe input_files/shapes.txt
.\ray_tracer.exe input_files/cornell_diffuse.txt
.\ray_tracer.exe input_files/mirror_sphere.txt
.\ray_tracer.exe input_files/shadow_test.txt
.\ray_tracer.exe input_files/recursive_test.txt
.\ray_tracer.exe input_files/full_test.txt

To create your own input file, note that the parsing was tailored for a specific format. You can learn the conventions of the format by reading the provided input files and the parser code within the main function of main.cpp. A few details to keep in mind when writing an input file:

- The winding order used for triangles is counter-clockwise.
- Spheres must specify the position and radius in the same line as: POS x y z RADIUS r
- Quadrilaterals are also specified by providing 3 points. These points must also follow counter-clockwise winding order, since a quadrilateral is created by creating 2 triangles. The missing unspecified point of the quadrilateral is infered by adding the vector BA to point C.
- The recursive decay coefficient can be specified with a line: KSPEC 0.8
- Aspect ratio may not be provided in the input file, instead, resolution is expected as RESOLUTION x y. Default is 800x800 if not specified.
- Refer to parser code starting on line 335 for more details.

#include "raytracer.h"
#include "BVH.h"
#include <iostream>
#include <ctime>
#include <chrono>
int main(int argc, char* argv[]) {
    // Check if a file was provided on the command line
    // printf("Size of Triange: %lu\n", sizeof(Triangle));
    // printf("Size of Triange*: %lu\n", sizeof(Triangle*));
    if(argc != 2) {
        std::cout << "Usage: ./rayTracer -scenefile" << std::endl;
        exit(1);
    }
    // A file was provided, so open it up to extract scene info
    RayTracer tracer;
    tracer.InitFromFile(argv[1]);
    clock_t s = std::clock();
    tracer.RayTrace();
    clock_t e = std::clock();
    printf("Render Time: %.4fs\n", (float)(e - s)/ CLOCKS_PER_SEC);
    return 0;
}
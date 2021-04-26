#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>
#include <algorithm>
#include <utility>
#include "image_lib.h"
#include "PGA_3D.h"
using namespace std;

struct MaterialGL {
    float ka[3];
    float padding1;
    float kd[3];
    float padding2;
    float ks[3];
    float padding3;
    float kt[3];
    float padding4;
    // 3 vec4
    float ns;
    float ior;
    float padding[2];
    // 2 vec
    MaterialGL() {
        ka[0] = 1.0f;
        ka[1] = 1.0f;
        ka[2] = 1.0f;
        kd[0] = 1.0f;
        kd[1] = 1.0f;
        kd[2] = 1.0f;
        ks[0] = 0.0f;
        ks[1] = 0.0f;
        ks[2] = 0.0f;
        kt[0] = 0.0f;
        kt[1] = 0.0f;
        kt[2] = 0.0f;
        ns = 5.0f;
        ior = 1.0f;
    };
};

struct TriangleGL {
    float p1[3]; // 1 vec3
    float padding1;
    float p2[3]; // 1 vec3
    float padding2;
    float p3[3]; // 1 vec3
    float padding3;
    float n1[3]; // 1 vec3
    float padding4;
    float n2[3]; // 1 vec3
    float padding5;
    float n3[3]; // 1 vec3
    float padding6;
    // right now every triangle has a material, I would like to change this
    // to an offset in a material array, but I'm experiencing alignment issues
    // when I replace this with an int
    MaterialGL mat;
};

struct DimensionGL {
    float min[4];
    float max[4];
    //float padding[4]; // This needs to pad the struct to the nearest offset
};

struct nodeGL {
    DimensionGL dim;
    //float padding1; // This needs to pad the struct
    int l_child_offset;
    int r_child_offset;
    int triangle_offset;
    float padding2; // this also needs to pad the struct
};

/**
 * LightGL - this struct represents the light struct which is defined in the GLSL
 * compute shader. The float arrays ensure the struct is tightly packed.
 */ 
struct LightGL {
    float pos[4];
    float dir[4];
    float clr[4];
    int type[4];
    // int fluff[3];
    // assume point light,
    // try adding a type specifier later
};

struct Dimension {
    float min_x;
    float min_y;
    float min_z;
    float max_x;
    float max_y;
    float max_z;

    Dimension() {
        min_x = INFINITY;
        min_y = INFINITY;
        min_z = INFINITY;
        max_x = -INFINITY;
        max_y = -INFINITY;
        max_z = -INFINITY;
    };
    
    Dimension(const Dimension& d) {
        min_x = d.min_x;
        min_y = d.min_y;
        min_z = d.min_z;
        max_x = d.max_x;
        max_y = d.max_y;
        max_z = d.max_z;
    };
    Dimension Union(const Dimension& rhs) {
        Dimension join;
        join.min_x = min(rhs.min_x, min_x);
        join.min_y = min(rhs.min_y, min_y);
        join.min_z = min(rhs.min_z, min_z);
        join.max_x = max(rhs.max_x, max_x);
        join.max_y = max(rhs.max_y, max_y);
        join.max_z = max(rhs.max_z, max_z);
        return join;
    };
    pair<float, float> operator[](int idx) {
        if (idx == 0) {
            // return 0;
            return pair<float, float>(min_x, max_x);
        }
        else if (idx == 1) {
            //return 1;
            return pair<float, float>(min_y, max_y);
        }
        else {
            //return 2;
            return pair<float, float>(min_z, max_z);
        }
    };
};

#endif  // STRUCTS_H
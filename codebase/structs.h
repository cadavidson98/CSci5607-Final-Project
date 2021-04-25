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

struct Material {
    Color a_;   // ambient
    Color d_;   // diffuse
    Color s_;   // specular
    Color t_;   // transmissive
    float ns_;  // cosine power thingy
    float ior_; // index of refraction
    // Default is a matte white
    Material(Color a=Color(0,0,0), Color d=Color(1,1,1), Color s=Color(0,0,0), 
            float ns=5, Color t=Color(0,0,0), float ior=1) : a_(a),
            d_(d), s_(s), ns_(ns), t_(t), ior_(ior) {};
};

struct Sphere {
    Point3D center_;
    float radius_;
    Material mat_;
    Sphere(Point3D c=Point3D(0,0,0), float r=0.0) : center_(c), radius_(r) {};
};

struct Triangle {
    Point3D p1_;
    Point3D p2_;
    Point3D p3_;

    Dir3D n1_;
    Dir3D n2_;
    Dir3D n3_;
    // pointer to hopefully offset massive memory
    // requirements
    Material *mat_;
};

struct AABB {
    Dimension dim;
    AABB * l_child;
    AABB * r_child;
    std::vector<Sphere> s_;
    // dynamically allocated to reduce memory requirements
    std::vector<Triangle*> t_;
    
    AABB() {
        l_child = nullptr;
        r_child = nullptr;
    }
};

struct Hit {
    Point3D pos;
    Dir3D norm;
    Material m;
    float t = INFINITY;

    friend bool operator<(const Hit& lhs, const Hit& rhs) {
        return lhs.t < rhs.t;
    };
};

enum type {point, dir, spot, emissive};

struct Light {
    type type_;
    Color clr_;
    // optional values
    Point3D pos_;
    Dir3D dir1_;
    float angle1_;
    float angle2_;
};

#endif  // STRUCTS_H
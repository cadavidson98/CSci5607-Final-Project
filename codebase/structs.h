#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>
#include "image_lib.h"
#include "PGA_3D.h"

inline Color operator+(Color lhs, Color rhs) {
    return Color(lhs.r+rhs.r, lhs.g+rhs.g, lhs.b+rhs.b);
}

inline Color operator*(Color lhs, Color rhs) {
    return Color(lhs.r * rhs.r, lhs.g*rhs.g, lhs.b*rhs.b);
}

inline Color operator*(Color lhs, float rhs) {
    return Color(lhs.r * rhs, lhs.g * rhs, lhs.b * rhs);
}

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
    Dimension(Dimension& d) {
        min_x = d.min_x;
        min_y = d.min_y;
        min_z = d.min_z;
        max_x = d.max_x;
        max_y = d.max_y;
        max_z = d.max_z;
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
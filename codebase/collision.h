#ifndef COLLISION_H
#define COLLISION_H

#include "PGA_3D.h"
#include "structs.h"
#include <algorithm>

inline bool TriangleIntersect(Point3D p, Dir3D d, Triangle t, Hit &hit) {
  // Get the plane normal
  Dir3D norm = cross((t.p3_ - t.p1_), (t.p2_ - t.p1_)).normalized();

  float dt = fabs(dot(d, norm));

  if(dt < .001) {
    // the ray is parallel, so return false
    return false;
  }

  float denom = dot(d, norm);
  float hit_t = dot(t.p1_ - p, norm) / denom;
  if(hit_t < 0) {
    // only go forward in the direction
    return false;
  }
  
  Point3D hit_pos = p + d * hit_t;
  // now check if it is inside the triangle
  float a, b, c;
  float tri_mag = vee(t.p1_, t.p2_, t.p3_).magnitude();
  a = vee(t.p2_, t.p3_, hit_pos).magnitude()/tri_mag;
  b = vee(t.p3_, t.p1_, hit_pos).magnitude()/tri_mag;
  c = vee(t.p1_, t.p2_, hit_pos).magnitude()/tri_mag;
  //float tri_mag = .5* cross(t.p3_ - t.p1_, t.p2_ - t.p1_).magnitude();
  //a = .5*cross(t.p3_ - hit_pos, t.p2_ - hit_pos).magnitude()/tri_mag;
  //b = .5*cross(t.p3_ - hit_pos, t.p1_ - hit_pos).magnitude()/tri_mag;
  //c = .5*cross(t.p1_ - hit_pos, t.p2_ - hit_pos).magnitude()/tri_mag;

  if(a <= 1.0001 && b <= 1.0001 && c <= 1.0001 && a+b+c <= 1.0001) {
    hit.pos = hit_pos;
    hit.t = hit_t;
    hit.norm = (a*t.n1_ + b*t.n2_ + c*t.n3_).normalized();
    //hit.norm = (cross(t.p2_ - t.p1_, t.p3_ - t.p1_)).normalized();
    return true;
  }
  return false;
}

inline bool sphereIntersect(Point3D p_not, Dir3D v, Sphere s, Hit &hit) {
  Dir3D to_c = (p_not - s.center_);
  float a = dot(v, v);
  float b = 2.0 * dot(v, to_c);
  float c = dot(to_c, to_c) - s.radius_ * s.radius_;
  float root = b * b - 4.0 * a * c;
  if (root < 0) {
    return false;
  }
  float t_1 = (-b + sqrt(root)) / (2.0*a);
  float t_2 = (-b - sqrt(root)) / (2.0*a);
  if(t_1 > 0.0 && t_2 <= 0.0) {
    hit.t = t_1;
    hit.pos = p_not + t_1 * v;
    return true;
  }
  else if(t_2 > 0.0 && t_1 <= 0.0) {
    hit.t = t_2;  
    hit.pos = p_not + t_2 * v;
    return true;
  }
  if(t_1 > 0.0 && t_2 > 0.0) {
    hit.t = std::min(t_1, t_2);
    hit.pos = p_not + hit.t * v;
    return true;
  }
  return false;
}

#endif  // COLLISION_H
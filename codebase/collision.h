#ifndef COLLISION_H
#define COLLISION_H

#include "PGA_3D.h"
#include "structs.h"
#include <algorithm>

inline bool TriangleIntersect(Point3D p, Dir3D d, TriangleGL t, Hit &hit) {
  // Get the plane normal
  Point3D p1(t.p1[0], t.p1[1], t.p1[2]), p2(t.p2[0], t.p2[1], t.p2[2]), p3(t.p3[0], t.p3[1], t.p3[2]);
  Dir3D norm = cross((p3 - p1), (p2 - p1)).normalized();

  float dt = fabs(dot(d, norm));

  if(dt < .001) {
    // the ray is parallel, so return false
    return false;
  }

  float denom = dot(d, norm);
  float hit_t = dot(p1 - p, norm) / denom;
  if(hit_t < 0) {
    // only go forward in the direction
    return false;
  }
  
  Point3D hit_pos = p + d * hit_t;
  // now check if it is inside the triangle
  float a, b, c;
  float tri_mag = vee(p1, p2, p3).magnitude();
  a = vee(p2, p3, hit_pos).magnitude()/tri_mag;
  b = vee(p3, p1, hit_pos).magnitude()/tri_mag;
  c = vee(p1, p2, hit_pos).magnitude()/tri_mag;
  
  if(a <= 1.0001 && b <= 1.0001 && c <= 1.0001 && a+b+c <= 1.0001) {
      Dir3D n1(t.n1[0], t.n1[1], t.n1[2]), n2(t.n2[0], t.n2[1], t.n2[2]), n3(t.n3[0], t.n3[1], t.n3[2]);
      hit.pos = hit_pos;
      hit.t = hit_t;
      hit.norm = (a*n1 + b*n2 + c*n3).normalized();
      hit.norm = (dot(hit.norm, d) > 0.0) ? -1*hit.norm: hit.norm;
      //hit.m = t.mat;
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

inline bool AABBIntersect(Point3D p, Dir3D d, Dimension AABB, float &i_t) {
Dimension dim = AABB;
    float tmin(-INFINITY), tmax(INFINITY);
        float tx1 = (dim.min_x - p.x)/d.x;
        float tx2 = (dim.max_x - p.x)/d.x;
        tmin = std::max(tmin, std::min(tx1, tx2));
        tmax = std::min(tmax, std::max(tx1, tx2));
        float ty1 = (dim.min_y - p.y)/d.y;
        float ty2 = (dim.max_y - p.y)/d.y;
        tmin = std::max(tmin, std::min(ty1, ty2));
        tmax = std::min(tmax, std::max(ty1, ty2));
        float tz1 = (dim.min_z - p.z)/d.z;
        float tz2 = (dim.max_z - p.z)/d.z;
        tmin = std::max(tmin, std::min(tz1, tz2));
        tmax = std::min(tmax, std::max(tz1, tz2));
    if(tmax > 0 && tmax >= tmin) {
        i_t = (tmin >= 0) ? tmin : tmax;
        return true;
    }
    return false;
}
#endif  // COLLISION_H
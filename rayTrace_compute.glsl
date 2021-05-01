// An opengl compute shader for RayTracing
#version 430

layout(local_size_x = 10, local_size_y = 10) in;

struct Material {
    vec3 ka;
    vec3 kd;
    vec3 ks;
    vec3 kt;
    float ns;
    float ior;
};

struct HitInfo {
    vec3 pos;
    vec3 norm;
    float time;
    Material mat;
    bool hit;
};

struct Triangle {
    vec3 p1;
    vec3 p2;
    vec3 p3;
    vec3 n1;
    vec3 n2;
    vec3 n3;
    Material mat;
};

struct Light {
    vec3 pos;
    vec3 dir;
    vec3 clr;
    ivec3 type;
};

struct Dimension {
  vec3 min_pt;
  vec3 max_pt;
};

struct Node {
  Dimension dim;
  int l_child;
  int r_child;
  int tri_offset;
};

// these are structs defined for easy organization
layout(binding = 0, rgba32f) uniform writeonly image2D result;
layout(binding = 1, std430) buffer triangles
{
    Triangle tris[];
};
layout(binding = 2, std430) buffer light_buf {
    Light lights[];
};
layout(binding = 3, std430) buffer bvh_scene {
    Node nodes[];
};

// these are for the camera
uniform vec3 eye;
uniform vec3 forward;
uniform vec3 right;
uniform vec3 up;
uniform vec3 background_clr;

uniform int num_tris;
uniform int num_lights;
// these are the precomputed values
uniform float d;
uniform float width;
uniform float height;
uniform float half_width;
uniform float half_height;

void rayRecurse(in vec3 pos, in vec3 dir, in int depth, out vec4 color);
void rayRecurse2(in vec3 pos, in vec3 dir, in int depth, out vec4 color);
void rayRecurse3(in vec3 pos, in vec3 dir, in int depth, out vec4 color);
void rayRecurse4(in vec3 pos, in vec3 dir, in int depth, out vec4 color);
void rayRecurse5(in vec3 pos, in vec3 dir, in int depth, out vec4 color);
void sceneIntersect(in vec3 pos, in vec3 dir, inout HitInfo hit);
void triangleIntersect(in vec3 pos, in vec3 dir, in Triangle tri, inout HitInfo hit);
void AABBIntersect(in vec3 pos, in vec3 dir, in Dimension dim, inout HitInfo hit);
void lightPoint(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec4 color);
void light(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec4 color);

void main () {
    ivec2 id = ivec2(gl_GlobalInvocationID.xy);
    // compute pixel offset
    float u = half_width - (id.x + 0.5);
    float v = half_height - (id.y + 0.5);
    vec3 ray_pt = eye - d * (forward) + u * (right) + v * (up);
    
    vec4 clr = vec4(0,0,0,1);
    
    vec3 ray_dir = (ray_pt - eye);
    HitInfo hit;
    hit.hit = false;
    hit.time = 1.0 / 0.0;
    rayRecurse(eye, ray_dir, 1, clr);
    
    imageStore(result, id, clr);
}

void rayRecurse(in vec3 pos, in vec3 dir, in int depth, out vec4 color) {
  HitInfo hit;
  hit.time = 1.0 / 0.0;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  vec4 clr = vec4(background_clr, 1.0);
  if (hit.hit) {
    light(hit.pos, dir, hit.norm, hit.mat, clr);
    if(hit.mat.ks.x + hit.mat.ks.y + hit.mat.ks.z > 0.0) {
      HitInfo reflect_hit;
      reflect_hit.hit = false;
      reflect_hit.time = 1.0 / 0.0;
      vec4 reflect_clr = vec4(0,0,0,1);
      vec3 r = normalize(reflect(dir, normalize(hit.norm)));
      vec3 wiggle = hit.pos + .001 * (r);
      rayRecurse2(wiggle, r, 2, reflect_clr);
      clr = vec4(clr.rgb + hit.mat.ks * reflect_clr.rgb, 1);
      //sceneIntersect(wiggle, r, reflect_hit);
      //if(reflect_hit.hit) {
      //  light(reflect_hit.pos, r, reflect_hit.norm, reflect_hit.mat, reflect_clr);
      //  reflect_clr = vec4(hit.mat.ks * reflect_clr.rgb, 0);
      //  clr = vec4(clr.rgb + reflect_clr.rgb, 1);
      //}
    }
  }
  color = clr;
}

void rayRecurse2(in vec3 pos, in vec3 dir, in int depth, out vec4 color) {
  HitInfo hit;
  hit.time = 1.0 / 0.0;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  vec4 clr = vec4(background_clr, 1.0);
  if (hit.hit) {
    light(hit.pos, dir, hit.norm, hit.mat, clr);
    if(hit.mat.ks.x + hit.mat.ks.y + hit.mat.ks.z > 0.0) {
      HitInfo reflect_hit;
      reflect_hit.hit = false;
      reflect_hit.time = 1.0 / 0.0;
      vec4 reflect_clr = vec4(0,0,0,1);
      vec3 r = normalize(reflect(dir, normalize(hit.norm)));
      vec3 wiggle = hit.pos + .001 * (r);
      rayRecurse3(wiggle, r, 2, reflect_clr);
      clr = vec4(clr.rgb + hit.mat.ks * reflect_clr.rgb, 1);
      //sceneIntersect(wiggle, r, reflect_hit);
      //if(reflect_hit.hit) {
      //  light(reflect_hit.pos, r, reflect_hit.norm, reflect_hit.mat, reflect_clr);
      //  reflect_clr = vec4(hit.mat.ks * reflect_clr.rgb, 0);
      //  clr = vec4(clr.rgb + reflect_clr.rgb, 1);
      //}
    }
  }
  color = clr;
}

void rayRecurse3(in vec3 pos, in vec3 dir, in int depth, out vec4 color) {
  HitInfo hit;
  hit.time = 1.0 / 0.0;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  vec4 clr = vec4(background_clr, 1.0);
  if (hit.hit) {
    light(hit.pos, dir, hit.norm, hit.mat, clr);
    if(hit.mat.ks.x + hit.mat.ks.y + hit.mat.ks.z > 0.0) {
      HitInfo reflect_hit;
      reflect_hit.hit = false;
      reflect_hit.time = 1.0 / 0.0;
      vec4 reflect_clr = vec4(0,0,0,1);
      vec3 r = normalize(reflect(dir, normalize(hit.norm)));
      vec3 wiggle = hit.pos + .001 * (r);
      rayRecurse4(wiggle, r, 2, reflect_clr);
      clr = vec4(clr.rgb + hit.mat.ks * reflect_clr.rgb, 1);
      //sceneIntersect(wiggle, r, reflect_hit);
      //if(reflect_hit.hit) {
      //  light(reflect_hit.pos, r, reflect_hit.norm, reflect_hit.mat, reflect_clr);
      //  reflect_clr = vec4(hit.mat.ks * reflect_clr.rgb, 0);
      //  clr = vec4(clr.rgb + reflect_clr.rgb, 1);
      //}
    }
  }
  color = clr;
}

void rayRecurse4(in vec3 pos, in vec3 dir, in int depth, out vec4 color) {
  HitInfo hit;
  hit.time = 1.0 / 0.0;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  vec4 clr = vec4(background_clr, 1.0);
  if (hit.hit) {
    light(hit.pos, dir, hit.norm, hit.mat, clr);
    if(hit.mat.ks.x + hit.mat.ks.y + hit.mat.ks.z > 0.0) {
      HitInfo reflect_hit;
      reflect_hit.hit = false;
      reflect_hit.time = 1.0 / 0.0;
      vec4 reflect_clr = vec4(0,0,0,1);
      vec3 r = normalize(reflect(dir, normalize(hit.norm)));
      vec3 wiggle = hit.pos + .001 * (r);
      rayRecurse5(wiggle, r, 2, reflect_clr);
      clr = vec4(clr.rgb + hit.mat.ks * reflect_clr.rgb, 1);
      //sceneIntersect(wiggle, r, reflect_hit);
      //if(reflect_hit.hit) {
      //  light(reflect_hit.pos, r, reflect_hit.norm, reflect_hit.mat, reflect_clr);
      //  reflect_clr = vec4(hit.mat.ks * reflect_clr.rgb, 0);
      //  clr = vec4(clr.rgb + reflect_clr.rgb, 1);
      //}
    }
  }
  color = clr;
}


void rayRecurse5(in vec3 pos, in vec3 dir, in int depth, out vec4 color) {
  color = vec4(0,0,0,1);
}
/**
 * Iteratively traverse the BVH using DFS to find a triangle collision.
 * Since we don't have any fancy data structures in GLSL, we
 * shall represent a stack using a finite array. This may seem
 * like a troublesome idea at first, but recall that DFS has a space complexity
 * of O(d). Therefore even if we limit the stack to a finite size of 20, we can
 * support up to 1,048,576 triangles!
 */
    
void sceneIntersect(in vec3 pos, in vec3 dir, inout HitInfo hit) {
    int index = 0;
    int stack[20];
    stack[0] = 0;
    HitInfo cur_hit;
    cur_hit.hit = false;
    cur_hit.time = 1.0 / 0.0;
    while(index >= 0) {
      // pop off the "top" node
      int cur_node_idx = stack[index];
      index = index - 1;
      if(cur_node_idx == -1) {
        // this is a "null" node
        continue;
      }
      Node cur_node = nodes[cur_node_idx];
      // Check if the node intersects the ray
      HitInfo box_hit;
      box_hit.hit = false;
      box_hit.time = 1.0 / 0.0;
      AABBIntersect(pos, dir, cur_node.dim, box_hit);
      if(!box_hit.hit) {
        continue;
      }
      if(cur_node.l_child == -1 && cur_node.r_child == -1) {
        // This is a leaf node, so check if the ray intersects the triangle
        HitInfo tri_hit;
        tri_hit.time = 1.0/0.0;
        tri_hit.hit = false;
        triangleIntersect(pos, dir, tris[cur_node.tri_offset], tri_hit);
        if(tri_hit.hit && tri_hit.time < hit.time) {
          // this triangle is closest, so keep track of it
          hit.hit = true;
          hit.norm = tri_hit.norm;
          hit.pos = tri_hit.pos;
          hit.time = tri_hit.time;
          hit.mat = tri_hit.mat;
        }
      }
      else {
        // add the child nodes to the stack (right child first then left)
        index = index + 1;
        stack[index] = cur_node.r_child;
        index = index + 1;
        stack[index] = cur_node.l_child;
      }
    }
}

/**
 * Check if the ray (pos, dir) intersects the AABB dim, and store the results in hit
 */
void AABBIntersect(in vec3 pos, in vec3 dir, in Dimension dim, inout HitInfo hit) {
  float tmin = -1.0 / 0.0;
  float tmax = 1.0 / 0.0;
  float tx1 = (dim.min_pt.x - pos.x) / dir.x;
  float tx2 = (dim.max_pt.x - pos.x) / dir.x;
  tmin = max(tmin, min(tx1, tx2));
  tmax = min(tmax, max(tx1, tx2));

  float ty1 = (dim.min_pt.y - pos.y) / dir.y;
  float ty2 = (dim.max_pt.y - pos.y) / dir.y;
  tmin = max(tmin, min(ty1, ty2));
  tmax = min(tmax, max(ty1, ty2));

  float tz1 = (dim.min_pt.z - pos.z) / dir.z;
  float tz2 = (dim.max_pt.z - pos.z) / dir.z;
  tmin = max(tmin, min(tz1, tz2));
  tmax = min(tmax, max(tz1, tz2));

  if(tmax > 0 && tmax >= tmin) {
    if(tmin >= 0) {
      hit.time = tmin;
    }
    else {
      hit.time = tmax;
    }
    hit.hit = true;
    return;
  }
  hit.hit = false;
  return;
}

void triangleIntersect(in vec3 pos, in vec3 dir, in Triangle tri, inout HitInfo hit) {
    // get the plane normal
    vec3 to_plane = vec3(tri.p1.x - pos.x, tri.p1.y - pos.y, tri.p1.z - pos.z);
    vec3 norm  = cross(tri.p3 - tri.p1, tri.p2 - tri.p1);
    float denom = dot(norm, dir);
    // the ray is parallel, so return nothing
    if (abs(denom) < .001) {
      hit.hit = false;
      hit.time = 1.0;
      return;
    }
    float t = dot(to_plane, norm) / denom;
    hit.time = t;
    hit.pos = pos + dir * hit.time;
    if (t < 0.0) {
      hit.hit = false;
      hit.time = 2.0;
      return;
    }
    vec3 plane_pos = vec3(pos.x + (dir.x * t), pos.y + (dir.y * t), pos.z + (dir.z * t));
    // no half cause it will cancel out 
    float tri_area = length(cross(tri.p2 - tri.p1, tri.p3 - tri.p1));
    float a = length(cross(tri.p3 - plane_pos, tri.p2 - plane_pos)) / tri_area;
    float b = length(cross(tri.p3 - plane_pos, tri.p1 - plane_pos)) / tri_area;
    float c = length(cross(tri.p1 - plane_pos, tri.p2 - plane_pos)) / tri_area;
    if(a <= 1.0001 && b <= 1.0001 && c <= 1.0001 && (a + b + c) <= 1.0001) {
      hit.hit = true;
      hit.time = (t);
      hit.pos = plane_pos;
      hit.norm = normalize(a * tri.n1 + b * tri.n2 + c * tri.n3);
      if(dot(hit.norm, dir) > 0.0) {
        hit.norm = -1.0 * hit.norm;
      }
      hit.mat = tri.mat;
      return;
    }
    else {
      hit.time = 3.0;
      hit.hit = false;
      return;
    }
}

void lightPoint(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec4 color) {
  vec3 to_eye = normalize(eye - pos);
  vec3 ambient = vec3(0.25*mat.ka.x, 0.25*mat.ka.y, 0.25*mat.ka.z);
  vec3 tot_clr = ambient;
  for(int i = 0; i < num_lights; i++) {
    vec3 to_light = (lights[i].pos - pos);
    float dist = to_light.length();
    to_light = normalize(to_light);
    HitInfo shadow_hit;
    shadow_hit.time = 1.0 / 0.0;
    shadow_hit.hit = false;
    vec3 wiggle = vec3(pos.x + 0.1 * to_light.x, pos.y + 0.1 * to_light.y, pos.z + 0.1 * to_light.z);
    sceneIntersect(wiggle, to_light, shadow_hit);
    if(shadow_hit.hit && abs(shadow_hit.time) < dist) {
      continue;
    }
    vec3 n = normalize(norm);
    float dir = dot(n, to_light);
    vec3 r = normalize(reflect(to_light, n));
    float kd = max(0.0, dir);
    float ks = max(0.0, pow(dot(r, to_eye), 5));
    
    vec3 diffuse = 1.0/(dist*dist) * vec3(lights[i].clr.x*mat.kd.x*kd, lights[i].clr.y*mat.kd.y*kd, lights[i].clr.z*mat.kd.z*kd);
    vec3 specular = vec3(ks*mat.ks.x, ks*mat.ks.y, ks*mat.ks.z);
    tot_clr = tot_clr+diffuse+specular;
  }
  color = vec4(tot_clr, 1);
}

void light(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec4 color) {
  vec3 to_eye = normalize(eye - pos);
  vec3 ambient = vec3(0.25*mat.ka.x, 0.25*mat.ka.y, 0.25*mat.ka.z);
  vec3 tot_clr = ambient;
  for(int i = 0; i < num_lights; i++) {
    float attenuation = 1.0;
    float dist = 99999;
    vec3 to_light = vec3(0.0,0.0,0.0);
    int type = lights[i].type.x;
    switch (type){
      case 0: //point light
        to_light = (lights[i].pos - pos);
        dist = to_light.length();
        attenuation = 1.0/(dist*dist);
        break;
      case 1: //directional light
        to_light = -lights[i].dir;
        break;
    }
    to_light = normalize(to_light);
    HitInfo shadow_hit;
    shadow_hit.time = 1.0 / 0.0;
    shadow_hit.hit = false;
    vec3 wiggle = vec3(pos.x + 0.001 * to_light.x, pos.y + 0.001 * to_light.y, pos.z + 0.001 * to_light.z);
    sceneIntersect(wiggle, to_light, shadow_hit);
    if(shadow_hit.hit && shadow_hit.time < dist) {
      // in depth shadow testing
      // get distance from point of origin to shadow hit
      continue;
    }
    vec3 n = normalize(norm);
    float dir = dot(n, to_light);
    vec3 r = normalize(reflect(to_light, n));
    float kd = max(0.0, dir);
    float ks = max(0.0, pow(dot(r, to_eye), 5));
    
    vec3 diffuse = attenuation * vec3(lights[i].clr.x*mat.kd.x*kd, lights[i].clr.y*mat.kd.y*kd, lights[i].clr.z*mat.kd.z*kd);
    vec3 specular = vec3(ks*mat.ks.x, ks*mat.ks.y, ks*mat.ks.z);
    tot_clr = tot_clr+diffuse+specular;
  }
  color = vec4(tot_clr, 1);
}
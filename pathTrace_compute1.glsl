// An opengl compute shader for RayTracing
#version 430
#define MAX_T 99999
#define MIN_T 0.0001
#define USELIGHT
int MAX_RAND_VALUE = 100;

layout(local_size_x = 10, local_size_y = 10) in;

struct Material {
    vec3 ka;
    vec3 kd;
    vec3 ks;
    vec3 kt;
    vec3 ke;
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

int rand_pos = 0;
int rand_img_pos;
int max_img_pos;
vec4 rands;
// these are structs defined for easy organization
layout(binding = 0, rgba32f) uniform image2D result;
layout(binding = 4, rgba32f) uniform image2D random;
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

const int numSamples = 100;

vec3 TracePath0(in vec3 pos, in vec3 dir);
vec3 TracePath1(in vec3 pos, in vec3 dir);
vec3 TracePath2(in vec3 pos, in vec3 dir);
vec3 TracePath3(in vec3 pos, in vec3 dir);
vec3 TracePath4(in vec3 pos, in vec3 dir);
void sceneIntersect(in vec3 pos, in vec3 dir, inout HitInfo hit);
void triangleIntersect(in vec3 pos, in vec3 dir, in Triangle tri, inout HitInfo hit);
void AABBIntersect(in vec3 pos, in vec3 dir, in Dimension dim, inout HitInfo hit);
void lightPoint(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec4 color);
void light(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec3 color);
float rand01();

void main () {
    ivec2 id = ivec2(gl_GlobalInvocationID.xy);
    rand_img_pos = int(gl_LocalInvocationIndex);
    max_img_pos = int(gl_WorkGroupSize.x * gl_NumWorkGroups.x + gl_WorkGroupSize.y * gl_NumWorkGroups.y) - 1;
    // compute pixel offset
    
    // vec3 clr_sum = vec3(0.0,0.0,0.0);
    // for (int i=0; i<numSamples; i++){
      // float u = half_width - (id.x + rand01() - 0.5);
      // float v = half_height - (id.y + rand01() - 0.5);
      float u = half_width - (id.x + 0.5);
    float v = half_height - (id.y + 0.5);
      vec3 ray_pt = eye - d * (forward) + u * (right) + v * (up);

      
      float seed = u + v * 3.43121412313;
      vec3 clr = vec3(0.0,0.0,0.0);
      vec3 ray_dir = (ray_pt - eye);
      HitInfo hit;
      hit.hit = false;
      hit.time = MAX_T;
      clr = TracePath0(eye, ray_dir);
      // clr_sum=clr_sum+clr;
    // }
    // clr_sum = clr_sum * (1.0/float(numSamples));
    imageStore(result, id, vec4(clr,1.0));
}

//Sample uniform disk, then extrude for cosine weighted hemisphere
vec3 sampleHemisphereCosine(in float u1, in float u2){
    const float r = sqrt(u1);
    const float theta = 2.0 * 3.14 * u2;
    const float x = r * cos(theta);
    const float z = r * sin(theta);
    return vec3(x, sqrt(max(0.0f, 1.0 - u1)), z);
};

float rand01(){
  // get the random value in the random image
  if(rand_pos > 3) {
    rand_pos = 0;
    rand_img_pos = rand_img_pos + 1;
  }
  if(rand_img_pos > max_img_pos) {
    rand_img_pos = 0;
  }
  return rands[rand_pos];
}

// float rand(vec2 co){
//     return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
// }

//Change of basis from BTN space to world space
vec3 sampleAroundNormalCosine(in vec3 normal, in vec3 bitangent, in vec3 tangent){
  float r1 = rand01();
  float r2 = rand01();;
  vec3 s = sampleHemisphereCosine(r1, r2);
  return vec3( 
    dot(s,vec3(bitangent.x,normal.x,tangent.x)), 
    dot(s,vec3(bitangent.y,normal.y,tangent.y)), 
    dot(s,vec3(bitangent.z,normal.z,tangent.z)));
};

float luminance(in vec3 v){
  return 0.299*v.x + 0.587*v.y + 0.114*v.z;
}

vec3 TracePath0(in vec3 pos, in vec3 dir) {
  HitInfo hit;
  hit.time = MAX_T;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  // return vec3(1.0,0.0,0.0);
  if (!hit.hit) return background_clr;
  // return vec3(1.0,0.0,0.0);
  vec3 clr_sum = hit.mat.ke;

  vec3 clr_diffuse = vec3(0.0,0.0,0.0);
  vec3 clr_reflective = vec3(0.0,0.0,0.0);
  vec3 clr_refractive = vec3(0.0,0.0,0.0);

  float diffusive_l = luminance(hit.mat.kd);
  float reflective_l = luminance(hit.mat.ks);
  float refractive_l = luminance(hit.mat.kt);
  float l_sum = diffusive_l + reflective_l + refractive_l;
  float diffusive_a = diffusive_l/l_sum;
  float reflective_a = refractive_l/l_sum;
  float refractive_a = refractive_l/l_sum;

  if (diffusive_l>0.0){
    for (int i=0; i<numSamples; i++){
      vec3 t1 = normalize(cross(hit.norm,dir));//(-hi.hitNorm.z, 0.0, hi.hitNorm.x);
      vec3 t2 = normalize(cross(hit.norm,t1)); 
      vec3 dir_indirect = sampleAroundNormalCosine(hit.norm,t1,t2);
      clr_diffuse = clr_diffuse + TracePath1(hit.pos, dir_indirect);    
    }
    clr_diffuse = clr_diffuse * hit.mat.kd * diffusive_a * (1.0/float(numSamples));
  }

  if (reflective_l>0.0){
    vec3 dir_reflect = normalize(reflect(dir, normalize(hit.norm)));
    clr_reflective = TracePath1(hit.pos, dir_reflect) * hit.mat.ks * reflective_a; 
  }

  if (refractive_l>0.0){
    float ctheta = dot(dir, hit.norm);
    float eta = 1.0/hit.mat.ior;
    vec3 n = hit.norm;
    if (ctheta>0.0){
      ctheta = -ctheta;
      eta = 1.0/eta;
      n = -n;
    }
    float eta_sq = eta * eta;
    float rf_term = 1.0 - eta_sq + ctheta * ctheta * eta_sq;
    if (rf_term > 0.0) { //there is refration
        vec3 dir_refract = (-sqrt(rf_term) - ctheta * eta) * n + eta * dir;
        clr_refractive = TracePath1(hit.pos, dir_refract) * hit.mat.kt * refractive_a; 
    }
  }
  clr_sum = clr_sum + clr_diffuse + clr_reflective + clr_refractive;
  #ifdef USELIGHT
  vec3 clr;
  light(hit.pos, dir, hit.norm, hit.mat, clr);
  clr_sum = clr_sum + clr;
  #endif
  return clr_sum;
}
vec3 TracePath1(in vec3 pos, in vec3 dir) {
  HitInfo hit;
  hit.time = MAX_T;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  if (!hit.hit) return background_clr;
  
  vec3 clr_sum = hit.mat.ke;

  vec3 clr_diffuse = vec3(0.0,0.0,0.0);
  vec3 clr_reflective = vec3(0.0,0.0,0.0);
  vec3 clr_refractive = vec3(0.0,0.0,0.0);

  float diffusive_l = luminance(hit.mat.kd);
  float reflective_l = luminance(hit.mat.ks);
  float refractive_l = luminance(hit.mat.kt);
  float l_sum = diffusive_l + reflective_l + refractive_l;
  float diffusive_a = diffusive_l/l_sum;
  float reflective_a = refractive_l/l_sum;
  float refractive_a = refractive_l/l_sum;

  if (diffusive_l>0.0){
      vec3 t1 = normalize(cross(hit.norm,dir));//(-hi.hitNorm.z, 0.0, hi.hitNorm.x);
      vec3 t2 = normalize(cross(hit.norm,t1)); 
      vec3 dir_indirect = sampleAroundNormalCosine(hit.norm,t1,t2);
      clr_diffuse = clr_diffuse + TracePath2(hit.pos, dir_indirect);    

    clr_diffuse = clr_diffuse * hit.mat.kd * diffusive_a;// * (1.0/float(numSamples));
  }

  if (reflective_l>0.0){
    vec3 dir_reflect = normalize(reflect(dir, normalize(hit.norm)));
    clr_reflective = TracePath2(hit.pos, dir_reflect) * hit.mat.ks * reflective_a; 
  }

  if (refractive_l>0.0){
    float ctheta = dot(dir, hit.norm);
    float eta = 1.0/hit.mat.ior;
    vec3 n = hit.norm;
    if (ctheta>0.0){
      ctheta = -ctheta;
      eta = 1.0/eta;
      n = -n;
    }
    float eta_sq = eta * eta;
    float rf_term = 1.0 - eta_sq + ctheta * ctheta * eta_sq;
    if (rf_term > 0.0) { //there is refration
        vec3 dir_refract = (-sqrt(rf_term) - ctheta * eta) * n + eta * dir;
        clr_refractive = TracePath2(hit.pos, dir_refract) * hit.mat.kt * refractive_a; 
    }
  }
  clr_sum = clr_sum + clr_diffuse + clr_reflective + clr_refractive;
  #ifdef USELIGHT
  vec3 clr;
  light(hit.pos, dir, hit.norm, hit.mat, clr);
  clr_sum = clr_sum + clr;
  #endif
  return clr_sum;
}
vec3 TracePath2(in vec3 pos, in vec3 dir) {
  HitInfo hit;
  hit.time = MAX_T;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  if (!hit.hit) return background_clr;
  
  vec3 clr_sum = hit.mat.ke;

  vec3 clr_diffuse = vec3(0.0,0.0,0.0);
  vec3 clr_reflective = vec3(0.0,0.0,0.0);
  vec3 clr_refractive = vec3(0.0,0.0,0.0);

  float diffusive_l = luminance(hit.mat.kd);
  float reflective_l = luminance(hit.mat.ks);
  float refractive_l = luminance(hit.mat.kt);
  float l_sum = diffusive_l + reflective_l + refractive_l;
  float diffusive_a = diffusive_l/l_sum;
  float reflective_a = refractive_l/l_sum;
  float refractive_a = refractive_l/l_sum;

  if (diffusive_l>0.0){
      vec3 t1 = normalize(cross(hit.norm,dir));//(-hi.hitNorm.z, 0.0, hi.hitNorm.x);
      vec3 t2 = normalize(cross(hit.norm,t1)); 
      vec3 dir_indirect = sampleAroundNormalCosine(hit.norm,t1,t2);
      clr_diffuse = clr_diffuse + TracePath3(hit.pos, dir_indirect);    

    clr_diffuse = clr_diffuse * hit.mat.kd * diffusive_a;// * (1.0/float(numSamples));
  }

  if (reflective_l>0.0){
    vec3 dir_reflect = normalize(reflect(dir, normalize(hit.norm)));
    clr_reflective = TracePath3(hit.pos, dir_reflect) * hit.mat.ks * reflective_a; 
  }

  if (refractive_l>0.0){
    float ctheta = dot(dir, hit.norm);
    float eta = 1.0/hit.mat.ior;
    vec3 n = hit.norm;
    if (ctheta>0.0){
      ctheta = -ctheta;
      eta = 1.0/eta;
      n = -n;
    }
    float eta_sq = eta * eta;
    float rf_term = 1.0 - eta_sq + ctheta * ctheta * eta_sq;
    if (rf_term > 0.0) { //there is refration
        vec3 dir_refract = (-sqrt(rf_term) - ctheta * eta) * n + eta * dir;
        clr_refractive = TracePath3(hit.pos, dir_refract) * hit.mat.kt * refractive_a; 
    }
  }
  clr_sum = clr_sum + clr_diffuse + clr_reflective + clr_refractive;
  #ifdef USELIGHT
  vec3 clr;
  light(hit.pos, dir, hit.norm, hit.mat, clr);
  clr_sum = clr_sum + clr;
  #endif
  return clr_sum;
}

vec3 TracePath3(in vec3 pos, in vec3 dir) {
  HitInfo hit;
  hit.time = MAX_T;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  if (!hit.hit) return background_clr;
  
  vec3 clr_sum = hit.mat.ke;

  vec3 clr_diffuse = vec3(0.0,0.0,0.0);
  vec3 clr_reflective = vec3(0.0,0.0,0.0);
  vec3 clr_refractive = vec3(0.0,0.0,0.0);

  float diffusive_l = luminance(hit.mat.kd);
  float reflective_l = luminance(hit.mat.ks);
  float refractive_l = luminance(hit.mat.kt);
  float l_sum = diffusive_l + reflective_l + refractive_l;
  float diffusive_a = diffusive_l/l_sum;
  float reflective_a = refractive_l/l_sum;
  float refractive_a = refractive_l/l_sum;

  if (diffusive_l>0.0){
      vec3 t1 = normalize(cross(hit.norm,dir));//(-hi.hitNorm.z, 0.0, hi.hitNorm.x);
      vec3 t2 = normalize(cross(hit.norm,t1)); 
      vec3 dir_indirect = sampleAroundNormalCosine(hit.norm,t1,t2);
      clr_diffuse = clr_diffuse + TracePath4(hit.pos, dir_indirect);    

    clr_diffuse = clr_diffuse * hit.mat.kd * diffusive_a;// * (1.0/float(numSamples));
  }

  if (reflective_l>0.0){
    vec3 dir_reflect = normalize(reflect(dir, normalize(hit.norm)));
    clr_reflective = TracePath4(hit.pos, dir_reflect) * hit.mat.ks * reflective_a; 
  }

  if (refractive_l>0.0){
    float ctheta = dot(dir, hit.norm);
    float eta = 1.0/hit.mat.ior;
    vec3 n = hit.norm;
    if (ctheta>0.0){
      ctheta = -ctheta;
      eta = 1.0/eta;
      n = -n;
    }
    float eta_sq = eta * eta;
    float rf_term = 1.0 - eta_sq + ctheta * ctheta * eta_sq;
    if (rf_term > 0.0) { //there is refration
        vec3 dir_refract = (-sqrt(rf_term) - ctheta * eta) * n + eta * dir;
        clr_refractive = TracePath4(hit.pos, dir_refract) * hit.mat.kt * refractive_a; 
    }
  }
  clr_sum = clr_sum + clr_diffuse + clr_reflective + clr_refractive;
  #ifdef USELIGHT
  vec3 clr;
  light(hit.pos, dir, hit.norm, hit.mat, clr);
  clr_sum = clr_sum + clr;
  #endif
  return clr_sum;
}

vec3 TracePath4(in vec3 pos, in vec3 dir) {
  return vec3(0.0,0.0,0.0);
}

// vec3 TracePath4(in vec3 pos, in vec3 dir) {
//   HitInfo hit;
//   hit.time = MAX_T;
//   hit.hit = false;
//   sceneIntersect(pos, dir, hit);
//   if (!hit.hit) return background_clr;
//   return hit.mat.ke;
// }


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
    cur_hit.time = MAX_T;
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
      box_hit.time = MAX_T;
      AABBIntersect(pos, dir, cur_node.dim, box_hit);
      if(!box_hit.hit) {
        continue;
      }
      if(cur_node.l_child == -1 && cur_node.r_child == -1) {
        // This is a leaf node, so check if the ray intersects the triangle
        HitInfo tri_hit;
        tri_hit.time = MAX_T;
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
  float tmin = -MAX_T;
  float tmax = MAX_T;
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

  if(tmax > MIN_T && tmax >= tmin) {
    if(tmin >= MIN_T) {
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
    if (abs(denom) < MIN_T) {
      hit.hit = false;
      hit.time = 1.0;
      return;
    }
    float t = dot(to_plane, norm) / denom;
    hit.time = t;
    hit.pos = pos + dir * hit.time;
    if (t < MIN_T) {
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
    shadow_hit.time = MAX_T;
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

void light(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec3 color) {
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
    shadow_hit.time = MAX_T;
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
  color = tot_clr;
}
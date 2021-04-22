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
void sceneIntersect(in vec3 pos, in vec3 dir, inout HitInfo hit);
void triangleIntersect(in vec3 pos, in vec3 dir, in Triangle tri, inout HitInfo hit);
void lightPoint(in vec3 pos, in vec3 dir, in vec3 norm, in Material mat, out vec4 color);

void main () {
    ivec2 id = ivec2(gl_GlobalInvocationID.xy);
    // compute pixel offset
    float u = half_width - (id.x + 0.5);
    float v = half_height - (id.y + 0.5);
    vec3 ray_pt = eye - d * (forward) + u * (right) + v * (up);
    
    vec4 clr1 = vec4(0,0,0,1);
    vec4 clr2 = vec4(0,0,0,1);
    vec4 clr3 = vec4(0,0,0,1);
    vec4 clr4 = vec4(0,0,0,1);
    vec4 clr5 = vec4(0,0,0,1);

    vec3 ray_dir = (ray_pt - eye);
    //vec3 ray_dir2 = eye - d * (forward) + (u - 0.5) * (right) + (v - 0.5) * (up);
    //vec3 ray_dir3 = eye - d * (forward) + (u + 0.5) * (right) + (v - 0.5) * (up);
    //vec3 ray_dir4 = eye - d * (forward) + (u - 0.5) * (right) + (v + 0.5) * (up);
    //vec3 ray_dir5 = eye - d * (forward) + (u + 0.5) * (right) + (v + 0.5) * (up);
    
    rayRecurse(eye, ray_dir, 1, clr1);
    //rayRecurse(eye, ray_dir, 1, clr2);
    //rayRecurse(eye, ray_dir, 1, clr3);
    //rayRecurse(eye, ray_dir, 1, clr4);
    //rayRecurse(eye, ray_dir, 1, clr5);
    imageStore(result, id, clr1);
}

void rayRecurse(in vec3 pos, in vec3 dir, in int depth, out vec4 color) {
  if(depth > 5) {
    return;
  }

  HitInfo hit;
  hit.time = 1.0 / 0.0;
  hit.hit = false;
  sceneIntersect(pos, dir, hit);
  vec4 clr = vec4(background_clr, 1.0);
  if (hit.hit) {
    lightPoint(hit.pos, dir, hit.norm, hit.mat, clr);
    if(hit.mat.ks.x + hit.mat.ks.y + hit.mat.ks.z > 0.0) {
      HitInfo reflect_hit;
      reflect_hit.hit = false;
      reflect_hit.time = 1.0 / 0.0;
      vec4 reflect_clr = vec4(0,0,0,1);
      vec3 r = normalize(reflect(dir, normalize(hit.norm)));
      vec3 wiggle = hit.pos + .1 * (r);
      sceneIntersect(wiggle, r, reflect_hit);
      if(reflect_hit.hit) {
        lightPoint(reflect_hit.pos, r, reflect_hit.norm, reflect_hit.mat, reflect_clr);
        reflect_clr = vec4(hit.mat.ks * reflect_clr.rgb, 0);
        clr = vec4(reflect_clr.rgb, 1);
      }
    }
  }
  color = clr;
}

void sceneIntersect(in vec3 pos, in vec3 dir, inout HitInfo hit) {
    HitInfo cur_hit;
    cur_hit.hit = false;
    cur_hit.time = 1.0 / 0.0;
    for(int i = 0; i < num_tris; i++) {
      triangleIntersect(pos, dir, tris[i], cur_hit);
      if(cur_hit.hit && cur_hit.time < hit.time) {
        hit.hit = cur_hit.hit;
        hit.norm = cur_hit.norm;
        hit.pos = cur_hit.pos;
        hit.time = cur_hit.time;
        hit.mat = cur_hit.mat;
      }
    }
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
#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "raytracer.h"
#include "collision.h"
#define USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <chrono>

#include <algorithm>

RayTracer::RayTracer(Point3D eye, Dir3D forward, Dir3D right, Dir3D up, float half_FOV) {
    eye_ = eye;
    forward_ = forward;
    right_ = right;
    up_ = up;
    HalfAngleFOV_ = half_FOV;
    max_depth_ = 5;
    verts_ = nullptr;
    norms_ = nullptr;
    scene = nullptr;
    img_ = nullptr;
}

RayTracer::~RayTracer() {
  if (scene != nullptr)
    delete scene;
  for(int i = 0; i < mats_.size(); i++) {
    delete mats_[i];
  }
  for(int i = 0; i < tris_.size(); i++) {
    delete tris_[i];
  }
  if (img_ != nullptr)
    delete img_;
}

void RayTracer::InitFromFile(std::string input_file_name) {
  // parse the file to get all the relevant information
  // TODO - save the spheres in the scene
  FILE* in_file;
  output_name_ = "rayTraced.bmp";
  width_ = 640; 
  height_ = 480;
  Material * cur_mat = new Material;
  mats_.push_back(cur_mat);
  in_file = fopen(input_file_name.c_str(), "r");
  if(in_file == NULL) {
    std::cerr << "Couldn't open file: " << input_file_name << std::endl;
    exit(1);
  }
  // The file was opened, so we are good to parse data now!
  clock_t s = std::clock();
  char arg[1024];
  Point3D min_pt(INFINITY, INFINITY, INFINITY), max_pt(-INFINITY, -INFINITY, -INFINITY);
  int max_vert(-1), max_norm(-1), nindex(0), vindex(0);
  while(fgets(arg, 1024, in_file)) {
    // Check if it is a comment, if so skip!
    // std::cout << arg << std::endl;
    if(arg[0] == '#') {
      continue;
    }
    // This contains valid data, so extract the command
    char command[100];
    int fieldsread = sscanf(arg, "%s ", command);
    std::string commandstr = command;
    if(fieldsread < 1) {
      continue;
    }

    if(commandstr == "sphere:") {
      // read in the sphere contents
      Point3D p;
      float r;
      sscanf(arg, "sphere: %f %f %f %f", &p.x, &p.y, &p.z, &r);
      Sphere s = Sphere(p, r);
      s.mat_ = *cur_mat;
      balls_.push_back(s);
      // get min and max corners of sphere and check scene dimensions
      Point3D mn = p + Dir3D(-r, -r, -r);
      Point3D mx = p + Dir3D(r, r, r);
      min_pt.x = std::min(mn.x, min_pt.x);
      min_pt.y = std::min(mn.y, min_pt.y);
      min_pt.z = std::min(mn.z, min_pt.z);

      max_pt.x = std::max(mx.x, max_pt.x);
      max_pt.y = std::max(mx.y, max_pt.y);
      max_pt.z = std::max(mx.z, max_pt.z);
    }
    else if(commandstr == "film_resolution:") {
      sscanf(arg, "film_resolution: %d %d", &width_, &height_);
    }
    else if(commandstr == "output_image:") {
      char name[100];
      sscanf(arg, "output_image: %s", name);
      output_name_ = name;
    }
    else if(commandstr == "camera_pos:") {
      sscanf(arg, "camera_pos: %f %f %f", &eye_.x, &eye_.y, &eye_.z);
    }
    else if(commandstr == "camera_fwd:") {
      sscanf(arg, "camera_fwd: %f %f %f", &forward_.x, &forward_.y, &forward_.z);
    }
    else if(commandstr == "camera_up:") {
      sscanf(arg, "camera_up: %f %f %f", &up_.x, &up_.y, &up_.z);
    }
    else if(commandstr == "camera_fov_ha:") {
      sscanf(arg, "camera_fov_ha: %f", &HalfAngleFOV_);
    }
    else if(commandstr == "max_depth:") {
      sscanf(arg, "max_depth: %d", &max_depth_);
    }
    else if(commandstr == "background:") {
      sscanf(arg, "background: %f %f %f", &background_clr_.r, &background_clr_.g, &background_clr_.b);
    }
    else if(commandstr == "material:") {
      float a_r, a_g, a_b, d_r, d_g, d_b, s_r, s_g, s_b, t_r, t_g, t_b, ns, ior;
      sscanf(arg, "material: %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &a_r, &a_g, &a_b, 
              &d_r, &d_g, &d_b, &s_r, &s_g, &s_b, &ns, &t_r, &t_g, &t_b, &ior);
      cur_mat = new Material;
      cur_mat->a_ = Color(a_r, a_g, a_b);
      cur_mat->d_ = Color(d_r, d_g, d_b);
      cur_mat->s_ = Color(s_r, s_g, s_b);
      cur_mat->t_ = Color(t_r, t_g, t_b);
      cur_mat->ns_ = ns;
      cur_mat->ior_ = ior;
      mats_.push_back(cur_mat);
    }
    else if(commandstr == "directional_light:" ) {
      Light l;
      float r, g, b, x, y, z;
      sscanf(arg, "directional_light: %f %f %f %f %f %f", &r, &g, &b, &x, &y, &z);
      l.clr_ = Color(r,g,b);
      l.dir1_ = Dir3D(x,y,z);
      l.type_ = dir;
      lights_.push_back(l);
    }
    else if(commandstr == "point_light:" ) {
      Light l;
      float r, g, b, x, y, z;
      sscanf(arg, "point_light: %f %f %f %f %f %f", &r, &g, &b, &x, &y, &z);
      l.clr_ = Color(r,g,b);
      l.pos_ = Point3D(x,y,z);
      l.type_ = point;
      lights_.push_back(l);
    }
    else if(commandstr == "spot_light:" ) {
      Light l;
      float r, g, b, px, py, pz, dx, dy, dz, a1, a2;
      sscanf(arg, "spot_light: %f %f %f %f %f %f %f %f %f %f %f", 
              &r, &g, &b, &px, &py, &pz, &dx, &dy, &dz, &a1, &a2);
      l.clr_ = Color(r,g,b);
      l.pos_ = Point3D(px, py, pz);
      l.dir1_ = Dir3D(dx,dy,dz);
      l.angle1_ = a1 * M_PI / 180.0;
      l.angle2_ = a2 * M_PI / 180.0;
      l.type_ = spot;
      lights_.push_back(l);
    }
    else if(commandstr == "ambient_light:" ) {
      float r, g, b;
      sscanf(arg, "ambient_light: %f %f %f", &r, &g, &b);
      amb_light_ = Color(r,g,b);
    }
    else if(commandstr == "max_vertices:") {
      sscanf(arg, "max_vertices: %d", &max_vert);
      verts_ = new Point3D[max_vert];
    }
    else if(commandstr == "max_normals:") {
      sscanf(arg, "max_normals: %d", &max_norm);
      norms_ = new Dir3D[max_norm];
    }
    else if(commandstr == "vertex:") {
      Point3D pt;
      sscanf(arg, "vertex: %f %f %f", &pt.x, &pt.y, &pt.z);
      // check scene dimensions
      min_pt.x = std::min(pt.x, min_pt.x);
      min_pt.y = std::min(pt.y, min_pt.y);
      min_pt.z = std::min(pt.z, min_pt.z);

      max_pt.x = std::max(pt.x, max_pt.x);
      max_pt.y = std::max(pt.y, max_pt.y);
      max_pt.z = std::max(pt.z, max_pt.z);
      
      verts_[vindex] = pt;
      vindex++;
    }
    else if(commandstr == "normal:") {
      Dir3D dir;
      sscanf(arg, "normal: %f %f %f", &dir.x, &dir.y, &dir.z);
      norms_[nindex] = dir;
      nindex++;
    }
    else if(commandstr == "triangle:") {
      if(max_vert < 0) {
        printf("ERROR: NUMBER OF VERTICES NOT SPECIFIED. SKIPPING\n");
        continue;
      }
      int p1, p2, p3;
      Triangle *t = new Triangle;
      sscanf(arg, "triangle: %d %d %d", &p1, &p2, &p3);
      t -> p1_ = verts_[p1];
      t -> p2_ = verts_[p2];
      t -> p3_ = verts_[p3];
      
      int area = 0;
      // get normal, per face normal
      int mag = (area < 0) ? -1 : 1;
      Dir3D n = mag * cross(t->p3_ - t->p1_, t->p2_ - t->p1_).normalized();
      t -> n1_ = n;
      t -> n2_ = n;
      t -> n3_ = n;
      t -> mat_ = cur_mat;
      tris_.push_back(t);
    }
    else if(commandstr == "normal_triangle:") {
      if(max_vert < 0 || max_norm < 0) {
        printf("ERROR: NUMBER OF VERTICES/NORMALS NOT SPECIFIED. SKIPPING\n");
        continue;
      }
      int p1, p2, p3, n1, n2, n3;
      Triangle *t = new Triangle;
      sscanf(arg, "normal_triangle: %d %d %d %d %d %d", &p1, &p2, &p3, &n1, &n2, &n3);
      t -> p1_ = verts_[p1];
      t -> p2_ = verts_[p2];
      t -> p3_ = verts_[p3];
      
      t -> n1_ = norms_[n1];
      t -> n2_ = norms_[n2];
      t -> n3_ = norms_[n3];
      t -> mat_ = cur_mat;
      tris_.push_back(t);
    }
  }
  right_ = cross(up_, forward_);
  up_ = cross(forward_, right_);
  forward_ = forward_.normalized();
  right_ = right_.normalized();
  up_ = up_.normalized();

  // now do the scene/BVH stuff
  Dimension scene_d;
  scene_d.min_x = min_pt.x;
  scene_d.min_y = min_pt.y;
  scene_d.min_z = min_pt.z;
  
  scene_d.max_x = max_pt.x;
  scene_d.max_y = max_pt.y;
  scene_d.max_z = max_pt.z;
  if(verts_ != nullptr)
    delete [] verts_;
  if(norms_ != nullptr)
    delete [] norms_;
  clock_t e = std::clock();
  printf("File Load Time: %.4fs\n", (float)(e - s)/ CLOCKS_PER_SEC);
  s = std::clock();  
  printf("File load complete. Beginning BVH construction\n");
  scene = new BVH(scene_d, tris_, balls_);
  e = std::clock();
  printf("BVH Init Time: %.4fs\n", (float)(e - s)/ CLOCKS_PER_SEC);
  fclose(in_file);
  img_ = new Image(width_, height_);
}

void RayTracer::RayTrace() {
    Image* result = img_;
    float half_width = width_ * .5;
    float half_height = height_ * .5;
    float d = half_height / tanf(HalfAngleFOV_ * (M_PI / 180.0));
    int width(width_), height(height_);
    int size = width * height;
    #pragma omp parallel for schedule(guided)
    for(int y = 0; y < height; y++) {
      #pragma omp parallel
      for(int x = 0; x < width; x++) {
        // generate ray
        // TODO - try u = -half_width * (x+.5)/width + (.5-x/width)*half_width
        // and compare the results
        //printf("\rRay %d of %d", x + y*width, size);
        //fflush(stdout);
        float u = half_width - (width * ((x + 0.5)/width));
        float v = half_height - (height * ((y + 0.5)/height));
        Point3D ray_point = eye_ - d*forward_ + u*right_+ v*up_;
        
        Dir3D ray_dir = (ray_point - eye_);
        
        // Get color
        result->setPixel(x, y, RayRecurse(eye_, ray_dir, 0));
      }
    }
    const char * name = output_name_.c_str();
    result->write(name);
    //delete result;
}

Color RayTracer::RayRecurse(Point3D p, Dir3D d, int depth) {
  if(depth > max_depth_) {
    // It is black not background for a reason
    return Color(0,0,0);
  }
  // shoot ray from p in direction dir
  Color color = background_clr_;
  Hit hit_info;
  bool hit = SceneIntersect(p, d, hit_info);
  if(hit) {
    color = Color(0,0,0);
    Point3D r_pt = hit_info.pos + .0001 * hit_info.norm;
    for(int i = 0; i < lights_.size(); i++) {
      color = color + LightPoint(lights_[i], r_pt, hit_info.m, hit_info.norm);
    }
    Color ambient = amb_light_ * hit_info.m.a_;
    color = color + ambient;
    // now we do the recursion
    // reflect
    if(hit_info.m.s_.r + hit_info.m.s_.g + hit_info.m.s_.b > 0.0) {
      Dir3D reflect = d - (2*dot(hit_info.norm, d)*hit_info.norm);
      Color ret = RayRecurse(r_pt, reflect.normalized(), depth+1);
      color = color + (ret * hit_info.m.s_);
    }
    // refraction
    if(hit_info.m.t_.r + hit_info.m.t_.g + hit_info.m.t_.b > 0.0) {
      float dir = dot(d, hit_info.norm);
      Dir3D r_in;
      if(dir > 0) {
        // leaving
        r_in = Refract(d.normalized(), -1 * hit_info.norm, hit_info.m.ior_, 1.0);
      }
      if (dir < 0) {
        // entering
        r_in = Refract(d.normalized(), hit_info.norm, 1.0, hit_info.m.ior_);
      }
      // wiggle in
      Point3D in_s = hit_info.pos + .005 * r_in;
      // find the out point
      // recurse!
      Color refrac_clr = RayRecurse(in_s, r_in, depth+1);
      Color amt = hit_info.m.t_;
      color = color + Color(refrac_clr.r*amt.r, refrac_clr.g*amt.g, refrac_clr.b*amt.b);
    }
  }
  return color;
}

Dir3D RayTracer::Refract(Dir3D d, Dir3D n, float n_i, float n_r) {
  float m = n_i / n_r;
  float dot_prod = dot(d,n);
  //float sign = (dot_prod < 0.0) ? -1.0 : (dot_prod > 0.0) ? 1.0 : 0.0;
  Dir3D t = m * (d - n*dot_prod) - n*sqrt(1.0 - m*m*(1.0 - dot_prod*dot_prod));
  return t.normalized();
}

// Apply lighting to the point P using Phong Illumination,
// as well as reflection and refraction
Color RayTracer::LightPoint(Light l, Point3D p, Material m, Dir3D n) {
  Color clr(0,0,0);
  // get ambient light
    float time;
    float i = 1.0;  
    Dir3D to_light;
    Dir3D to_eye;
    // get light vectors
    switch(l.type_) {
      case point: {
      to_light = l.pos_ - p;
      to_eye = eye_ - p;
      i = 1.0 / (dot(to_light, to_light) + .0001);
      time = to_light.magnitudeSqr();
      break;
    }
      case dir: {
      to_light = (l.dir1_ * -1.0).normalized();
      to_eye = eye_ - p;
      time = INFINITY;
      break;
    }
      case spot: {
      // get the angle from the light and point
      to_light = l.pos_ - p;
      time = to_light.magnitudeSqr();
      float angle = acosf(dot(-1*to_light.normalized(), (l.dir1_).normalized()));
      if(angle < l.angle1_) {
        // behave like point light
        to_eye = eye_ - p;
        i = 1.0 / (dot(to_light, to_light) + .0001);
      }
      else if(angle > l.angle2_) {
        to_eye = eye_ - p;
        i = 0.0;
      }
      else {
        // linear falloff
        to_eye = eye_ - p;
        i = ((l.angle2_ - angle) / (l.angle2_ - l.angle1_)) * (1.0 / (dot(to_light, to_light) + .0001));
      }
    }
    }

    // now check if we are in shadow
    Hit shadow_hit;
    Point3D wiggle = p + to_light.normalized()*.1;
    bool in_shadow = SceneIntersect(wiggle, to_light.normalized(), shadow_hit);
    //  if in_shadow && time < dist
    
    if(in_shadow && (shadow_hit.pos - p).magnitudeSqr() < time) {
      // we are in the shadow AND we hit an object
      // before reaching the light
      return clr;
    }
    else {
      float dir = dot(n, to_light);
      Dir3D norm = (dir < 0) ? -1.0*n : n;
      Dir3D r  = (2 * dot(to_light, norm) * norm - to_light).normalized();
      float kd = i * std::max((float)0.0, dot(norm, to_light.normalized()));
      if(kd == 0.0) {
        int a = 0;
      }
      float ks = std::max((float)0.0, pow(dot(r, to_eye.normalized()), m.ns_));
      
      Color diffuse = l.clr_* (m.d_ * kd);
      Color specular = m.s_ * ks;
      clr = clr + diffuse + specular;
    }

  return clr;
}

bool RayTracer::SceneIntersectFast(Point3D p, Dir3D d) {
  Hit temp;
  return scene->BVHIntersect(p, d, temp);
  for(int i = 0; i < balls_.size(); i++) {
    if (sphereIntersect(p, d, balls_[i], temp)) {
      return true;
    }
  }
  for(int i = 0; i < tris_.size(); i++) {
    if (TriangleIntersect(p, d, *tris_[i], temp)) {
      return true;
    }
  }
  return false;
}

// Check if the ray intersects the scene. if so, save the results in
// i_pt and t (intersection point and at time t)
bool RayTracer::SceneIntersect(Point3D p, Dir3D d, Hit &hit) {
  bool result = false;
  Hit cur_hit;
  return scene->BVHIntersect(p, d, hit);
  for(int i = 0; i < balls_.size(); i++) {
    // Find the closest ball, not the first ball
    if (sphereIntersect(p, d, balls_[i], cur_hit) && cur_hit.t < hit.t) {
      hit.pos = cur_hit.pos;
      hit.t = cur_hit.t;
      hit.m = balls_[i].mat_;
      hit.norm = (hit.pos - balls_[i].center_).normalized();
      result = true;
    }
  }
  for(int i = 0; i < tris_.size(); i++) {
    if(TriangleIntersect(p, d, *tris_[i], cur_hit) && cur_hit.t < hit.t) {
      hit = cur_hit;
      hit.m = *(tris_[i]->mat_);
      result = true;
    }
  }
  return result;
}

// Check to make sure triangle intersect actually works
void RayTracer::TriangleSanityCheck() {
  printf("-----------BEGINNING SANITY CHECK-----------\n");
  Triangle t;
  t.p1_ = Point3D(0,1,0);
  t.p2_ = Point3D(-1,0,0);
  t.p3_ = Point3D(1,0,0);
  Hit hit;
  if(TriangleIntersect(Point3D(0,0,-1), Dir3D(0,0,1), t, hit)) {
    printf("Triangle intersects!\n");
  }
  else {
    printf("You shouldn't see this!\n");
  }
  if(TriangleIntersect(Point3D(0,0,-1), Dir3D(1,0,0), t, hit)) {
    printf("You shouldn't see this!\n");
  }
  else {
    printf("Triangle doesn't intersect!\n");
  }
  if(TriangleIntersect(Point3D(0,0,-1), Dir3D(1,1,1), t, hit)) {
    printf("Triangle intersects!\n");
  }
  else {
    printf("You shouldn't see this!\n");
  }
  if(TriangleIntersect(Point3D(0,0,-1), Dir3D(0,0,-1), t, hit)) {
    printf("You shouldn't see this!\n");
  }
  else {
    printf("Triangle doesn't intersect!\n");
  }
  
}
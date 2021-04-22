#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "image_lib.h"

#include "PGA_3D.h"
#include "structs.h"
#include "BVH.h"
#include <string>
#include <vector>

class RayTracer {
public:
    RayTracer(Point3D eye=Point3D(0,0,0), Dir3D forward=Dir3D(0,0,-1), Dir3D right=Dir3D(1,0,0), Dir3D up=Dir3D(0,1,0), float half_FOV=45.0);
    ~RayTracer();
    void InitFromFile(std::string input_file_name);
    void RayTrace();
    void TriangleSanityCheck();
    Image* getImg() { return img_; };
private:
    // This is where the result is drawn to
    int width_;
    int height_;
    std::string output_name_;
    // This is the necessary information for the
    // Camera
    float HalfAngleFOV_;
    Point3D eye_;
    Dir3D forward_;
    Dir3D right_;
    Dir3D up_;
    int max_depth_;
    
    Point3D *verts_;
    Dir3D *norms_;
    BVH *scene;
    Image* img_;
    // this is the spheres which may appear in the scene
    std::vector<Sphere> balls_;
    std::vector<Triangle*> tris_;
    std::vector<Material*> mats_;
    std::vector<Light> lights_;
    // this is lighting information
    Color background_clr_;
    Color amb_light_;

    Color RayRecurse(Point3D p, Dir3D d, int depth);
    Color LightPoint(Light l, Point3D p, Material m, Dir3D n);
    bool SceneIntersect(Point3D p, Dir3D d, Hit &hit);
    bool SceneIntersectFast(Point3D p, Dir3D d);
    Dir3D Refract(Dir3D d, Dir3D n, float n_i, float n_r);
};

#endif  // RAYTRACER_H
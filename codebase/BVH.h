#ifndef BVH_H
#define BVH_H

#include "structs.h"
#include <vector>

class BVH {
public:
    BVH(Dimension D, std::vector<Triangle*> &t, std::vector<Sphere> &s);
    ~BVH() { Destruct(scene); };
    bool BVHIntersect(Point3D p, Dir3D d, Hit &hit);
    static void SanityCheck();
    int tree_depth = 1;
private:
    void initBVH(Dimension D, std::vector<Triangle*> &t, std::vector<Sphere> &s);
    void buildBottomUp(std::vector<Triangle*> t, std::vector<Sphere> s);
    void buildRecurse(AABB* node, std::vector<AABB*> leaves);
    void initAABB(AABB *box, int depth, int prev_t);
    void Destruct(AABB *box);
    bool IntersectRecurse(AABB* box, Point3D p, Dir3D d, Hit &hit, int &checks);
    void GetAllIntersections(AABB* box, Point3D p, Dir3D d, std::vector<Hit> &hit);    
    
    int TriSide(float thresh, int axis, Triangle t);
    int SphereSide(float thresh, int axis, Sphere s);
    int PointSide(float thresh, int axis, Point3D p);

    bool contains(AABB *box, Triangle &t);
    bool contains(AABB *box, Sphere &s);
    bool contains(AABB *box, Point3D &p);
    bool RayIntersect(AABB *box, Point3D p, Dir3D d, float& i_t);
    bool BoxIntersect(AABB* box, Point3D p, Dir3D d, Hit& hit);
    AABB *scene;
    int max_elems = 15;
    int max_depth = 50;
    
    void printBVH(AABB* box, int depth);

};

#endif

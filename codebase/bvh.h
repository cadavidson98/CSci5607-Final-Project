#ifndef bvh_h
#define bvh_h

#include <vector>

#include "structs.h"

using namespace std;

struct node {
    Dimension AABB_;
    node * l_child_ = nullptr;
    node * r_child_ = nullptr;
    vector<int> tri_offsets_;
};

struct compressed_node {
    Dimension AABB_;
    int l_child_ = -1;
    int r_child_ = -1;
    vector<int> tri_offsets_;
};

struct triangle_info {
  Dimension AABB_;
  Point3D centroid_;
  int tri_offset_;
};

class bvh {
  public:
    bvh() {};
    bvh(vector<TriangleGL*> tris);
    bool intersect(Point3D p, Dir3D d, Hit& hit);
    nodeGL * getCompact(int & num_nodes);
  private:
    vector<TriangleGL*> triangles_;
    vector<compressed_node> bvh_nodes_;
    int max_tri = 4;

    Dimension getExtent(vector<node*> nodes);
    Dimension getExtent(vector<triangle_info> tris);
    Dimension getExtent(vector<Point3D> pts);
    float surfaceArea(Dimension AABB);
    vector<triangle_info> boundTriangles();
    bool splitSAH(vector<triangle_info> in_tris, vector<triangle_info> &bin_1, vector<triangle_info> &bin_2);
    bool splitMidpoint(vector<triangle_info> in_tris, vector<triangle_info> &bin_1, vector<triangle_info> &bin_2);
    void buildRecurse(int node_offset, vector<triangle_info> tris);
    bool intersectRecurse(int node_offset, Point3D p, Dir3D d, Hit& hit, float &best_time);
    bool intersectIterative(Point3D p, Dir3D d, Hit& hit);
    void validate(node * cur_node);
};

#endif  // bvh_h
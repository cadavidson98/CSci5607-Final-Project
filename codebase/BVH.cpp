#include "BVH.h"
#include "collision.h"
#include <string>
#include <algorithm>

void BVH::SanityCheck() {
    // Build and example BVH and make sure all the elements are accounted for
    Dimension d;
    d.min_x = 0;
    d.min_y = 0;
    d.min_z = 0;
    d.max_x = 100;
    d.max_y = 100;
    d.max_z = 3;
    std::vector<Triangle*> t;
    std::vector<Sphere> s;
    // generate 100 spheres
    for(int y = 0; y < 10; y++) {
        for(int x = 0; x < 5; x++) {
            Sphere sp;
            sp.center_ = Point3D(x*10 + 5, y*10 + 5, 2);
            sp.radius_ = 1;
            s.push_back(sp);
        }
    }
    printf("Building BVH\n");
    BVH bvh(d, t, s);
    printf("Done building, testing intersection\n");
    Point3D p(5,5,5), i_pt;
    Dir3D dir(0,0,-1), i_n;
    Hit hit_info;
    bool hit = bvh.BVHIntersect(p, dir, hit_info);
    std::string str = (hit) ? "yes" : "no";
    printf("does (0,0,-1) hit anything?\t");
    printf("%s", str.c_str());
    printf("\t at (%.2f,%.2f, %.2f)\n", i_pt.x, i_pt.y, i_pt.z);

    hit = bvh.BVHIntersect(Point3D(75, 5, 2), Dir3D(-1, 0, 0), hit_info);
    str = (hit) ? "yes" : "no";
    printf("does (0,0,1) hit anything?\t");
    printf("%s", str.c_str());
    
}

///---------------------------------------------------------------------------------///
///---------------------------------------------------------------------------------///
///---------------------------------------------------------------------------------///
BVH::BVH(Dimension d, std::vector<Triangle*> &t, std::vector<Sphere> &s) {
    scene = nullptr;
    if(t.size() + s.size() > 0) {
        initBVH(d, t, s);
    }
    //printBVH(scene, 0);
}

void BVH::initBVH(Dimension d, std::vector<Triangle*> &t, std::vector<Sphere> &s) {
    scene = new AABB;
    scene->dim = d;
    scene->s_ = s;
    scene->t_ = t;
    buildBottomUp(t, s);
    //initAABB(scene, 0, -1);
    
}

void BVH::initAABB(AABB* box, int depth, int prev_t) {
    int sum = box->t_.size() + box->s_.size();
    if(sum <= max_elems || depth >= max_depth || prev_t == box->t_.size()) {
        // This box is fine
        //printf("Leaf created at depth %d with %d triangles\n", depth, prev_t);
        tree_depth = (depth > tree_depth) ? depth: tree_depth; 
        return;
    }
    // calculate the dimension of this box
    for(int i=0; i < box->t_.size(); i++) {
        // get min and max bounds of triangle
        Triangle *t = box->t_[i];
        std::pair<float, float> x_bounds = std::minmax({t->p1_.x, t->p2_.x, t->p3_.x});
        std::pair<float, float> y_bounds = std::minmax({t->p1_.y, t->p2_.y, t->p3_.y});
        std::pair<float, float> z_bounds = std::minmax({t->p1_.z, t->p2_.z, t->p3_.z});

        box->dim.min_x = std::min(box->dim.min_x, x_bounds.first);
        box->dim.min_y = std::min(box->dim.min_y, y_bounds.first);
        box->dim.min_z = std::min(box->dim.min_z, z_bounds.first);
        
        box->dim.max_x = std::max(box->dim.max_x, x_bounds.second);
        box->dim.max_y = std::max(box->dim.max_y, y_bounds.second);
        box->dim.max_z = std::max(box->dim.max_z, z_bounds.second);
    }
    // find a split along the major axis
    Dimension d = box->dim; 
    float thresh;
    float arr[3];
    arr[0] = d.max_x - d.min_x;
    arr[1] = d.max_y - d.min_y;
    arr[2] = d.max_z - d.min_z;
    int index = -1;
    float max = -1.0;
    for(int i = 0; i < 3; i++) {
        if(arr[i] > max) {
            index++;
            max = arr[i];
        }
    }
    Dimension d1(d), d2(d);
    if(index == 0) {
        // split on x
        d1.max_x = d.min_x + .5 * arr[0];
        d2.min_x = d.min_x + .5 * arr[0];
        thresh = d.min_x + .5 * arr[0];
    }
    else if (index == 1) {
        // split on y
        d1.max_y = d.min_y + .5 * arr[1];
        d2.min_y = d.min_y + .5 * arr[1];
        thresh = d.min_y + .5 * arr[1];
    }
    else {
        // split on z
        d1.max_z = d.min_z + .5 * arr[2];
        d2.min_z = d.min_z + .5 * arr[2];
        thresh = d.min_z + .5 * arr[2];
    }
    
    box->l_child = new AABB;
    box->r_child = new AABB;
    box->l_child->dim = d1;
    box->r_child->dim = d2;
    
    // Check which side it belongs on
    int num_t = box->t_.size();
    for(auto iter = box->t_.begin(); iter != box->t_.end(); iter++) {
        Triangle *tri = (*iter);
        int side = TriSide(thresh, index, *tri);
        if(side == -1) {
            // add to left child only
            box->l_child->t_.push_back(tri);
        }
        else if( side == 1) {
            // add to right child only
            box->r_child->t_.push_back(tri);
        }
        else if( side == 0) {
            // add to both children
            box->l_child->t_.push_back(tri);
            box->r_child->t_.push_back(tri);
        }
        else {
            printf("ERROR! Element is not left, right, or in both boxes!\n");
        }
        box->t_.erase(iter);
        iter--;
    }
    int num_s = box->s_.size();
    for(auto iter = box->s_.begin(); iter != box->s_.end(); iter++) {
        Sphere sph = *iter;
        int side = SphereSide(thresh, index, sph);
        if(side == -1) {
            // add to left child only
            box->l_child->s_.push_back(sph);
        }
        else if( side == 1) {
            // add to right child only
            box->r_child->s_.push_back(sph);
        }
        else if( side == 0) {
            // add to both children
            box->l_child->s_.push_back(sph);
            box->r_child->s_.push_back(sph);
        }
        else {
            printf("ERROR! Element is not left, right, or in both boxes!\n");
        }
        iter = box->s_.erase(iter);
        // for loop sanity
        iter--;
    }
    // sanity check, is this box empty?
    box->t_.shrink_to_fit();
    box->s_.shrink_to_fit();
    int new_size = box->t_.size() + box->s_.size();
    if(new_size != 0) {
        printf("ERROR: not all elements were placed in a box!\n");
    }
    // recurse!
    initAABB(box->l_child, depth+1, num_t);
    initAABB(box->r_child, depth+1, num_t);
}

void BVH::printBVH(AABB* box, int depth) {
    if(box == nullptr) {
        return;
    }
    if(box -> l_child == nullptr && box->r_child == nullptr) {
        printf("Leaf at depth %d\n\tnumTriangles = %zd\n", depth, box->t_.size());
    }
    else {
        printBVH(box->l_child, depth+1);
        printBVH(box->r_child, depth+1);
    }
}

bool sort(AABB* box1, AABB* box2) {
    //get box center
    float c_x1(box1->dim.max_x - box1->dim.min_x), c_y1(box1->dim.max_y - box1->dim.min_y), c_z1(box1->dim.max_z - box1->dim.min_z);
    float c_x2(box2->dim.max_x - box2->dim.min_x), c_y2(box2->dim.max_y - box2->dim.min_y), c_z2(box2->dim.max_z - box2->dim.min_z);
    return (c_z1 != c_z2) ? c_z1 < c_z2 : (c_y1 != c_y2) ? c_y1 < c_y2 : c_x1 < c_x2;
}

void BVH::buildBottomUp(std::vector<Triangle*> t, std::vector<Sphere> s) {
    // start by building an AABB around every triangle
    int size = t.size();
    std::vector<AABB*> boxes;
    for(int i = 0; i < size; i++) {
        AABB* leaf_node = new AABB;
        leaf_node->t_.push_back(t[i]);
        //calculate the bounds
        Dimension leaf_d;
        std::pair<float, float> x_bounds = std::minmax({t[i]->p1_.x, t[i]->p2_.x, t[i]->p3_.x});
        std::pair<float, float> y_bounds = std::minmax({t[i]->p1_.y, t[i]->p2_.y, t[i]->p3_.y});
        std::pair<float, float> z_bounds = std::minmax({t[i]->p1_.z, t[i]->p2_.z, t[i]->p3_.z});
        
        leaf_d.min_x = x_bounds.first;
        leaf_d.min_y = y_bounds.first;
        leaf_d.min_z = z_bounds.first;
        
        leaf_d.max_x = x_bounds.second;
        leaf_d.max_y = y_bounds.second;
        leaf_d.max_z = z_bounds.second;
        leaf_node -> dim = leaf_d;
        boxes.push_back(leaf_node);
    }
    size = s.size();
    for(int i = 0; i < size; i++) {
        AABB* leaf_node = new AABB;
        //calculate the bounds
        leaf_node->s_.push_back(s[i]);
        leaf_node->dim.min_x = s[i].center_.x - s[i].radius_;
        leaf_node->dim.min_y = s[i].center_.y - s[i].radius_;
        leaf_node->dim.min_z = s[i].center_.z - s[i].radius_;

        leaf_node->dim.max_x = s[i].center_.x + s[i].radius_;
        leaf_node->dim.max_y = s[i].center_.y + s[i].radius_;
        leaf_node->dim.max_z = s[i].center_.z + s[i].radius_;
        boxes.push_back(leaf_node);
    }
    // now that we have all the leaf nodes, we just need to sort
    // get semi major axis
    scene = new AABB;
    buildRecurse(scene, boxes);
}

bool sortX(const AABB* box1, const AABB* box2) {
    return (box1->dim.min_x+box1->dim.max_x) < (box2->dim.min_x+box2->dim.max_x);
}

bool sortY(const AABB* box1, const AABB* box2) {
    return (box1->dim.min_y+box1->dim.max_y) < (box2->dim.min_y+box2->dim.max_y);
}

bool sortZ(const AABB* box1, const AABB* box2) {
    return (box1->dim.min_z+box1->dim.max_z) < (box2->dim.min_z+box2->dim.max_z);
}

void BVH::buildRecurse(AABB* box, std::vector<AABB*> leaves) {
    // calculate the box dimension
    if(leaves.size() == 1) {
        box -> l_child = leaves[0];
        return;
    }
    else if(leaves.size() == 2) {
        box->l_child = leaves[0];
        box->r_child = leaves[1];
        return;
    }
    for(int i = 0; i < leaves.size(); i++) {
        // get the bounds of each box
        box->dim.min_x = std::min(box->dim.min_x, leaves[i]->dim.min_x);
        box->dim.min_y = std::min(box->dim.min_y, leaves[i]->dim.min_y);
        box->dim.min_z = std::min(box->dim.min_z, leaves[i]->dim.min_z);

        box->dim.max_x = std::max(box->dim.max_x, leaves[i]->dim.max_x);
        box->dim.max_y = std::max(box->dim.max_y, leaves[i]->dim.max_y);
        box->dim.max_z = std::max(box->dim.max_z, leaves[i]->dim.max_z);
    }
    // get the semi major axis
    Dimension d = box->dim; 
    float thresh;
    float arr[3];
    arr[0] = d.max_x - d.min_x;
    arr[1] = d.max_y - d.min_y;
    arr[2] = d.max_z - d.min_z;
    int index = -1;
    float max = -1.0;
    for(int i = 0; i < 3; i++) {
        if(arr[i] > max) {
            index++;
            max = arr[i];
        }
    }
    // now split
    if(index == 0) {
        // split on x
        std::sort(leaves.begin(), leaves.end(), sortX);
    }
    else if (index == 1) {
        // split on y
        std::sort(leaves.begin(), leaves.end(), sortY);
    }
    else {
        // split on z
        std::sort(leaves.begin(), leaves.end(), sortZ);
    }

    // now give half the leaves to each child
    int half = leaves.size() / 2;
    box -> l_child = new AABB;
    buildRecurse(box->l_child, std::vector<AABB*>(leaves.begin(), leaves.begin()+half));
    box -> r_child = new AABB;
    buildRecurse(box->r_child, std::vector<AABB*>(leaves.begin()+half, leaves.end()));
}

void BVH::Destruct(AABB* box) {
    // bottom up destruction
    if(box == nullptr) {
        return;
    }
    Destruct(box->l_child);
    Destruct(box->r_child);
    delete box;
}

bool BVH::BVHIntersect(Point3D p, Dir3D d, Hit &hit) {
    // Start by trying to intersect the scene
    // keep intersecting nodes until we hit a leaf node
    // (the child nodes will be null)
    //std::vector<Hit> hits;
    float t = 100000;
    /*GetAllIntersections(scene, p, d, hits);
    for(int i = 0; i < hits.size(); i++) {
        // get first hit
        if(hits[i].t < t) {
            hit = hits[i];
            t = hit.t;
        }
    }*/
    int checks  = 0;
    bool result = IntersectRecurse(scene, p, d, hit, checks);
    return result;
}

void BVH::GetAllIntersections(AABB* box, Point3D p, Dir3D d, std::vector<Hit> &hit) {
    float t_0;
    bool hit_box = RayIntersect(box, p, d, t_0);
    if(!hit_box) {
        return;
    }
    if(box -> l_child == nullptr && box -> r_child == nullptr) {
        Hit hit_info;
        int size = box->t_.size();
        for(int i = 0; i < size; i++) {
            if(TriangleIntersect(p, d, *(box->t_[i]), hit_info)) {
                hit.push_back(hit_info);
            }
        }
        size = box->s_.size();
        for(int i = 0; i < size; i++) {
            if(sphereIntersect(p, d, box->s_[i], hit_info)) {
                hit.push_back(hit_info);
            }
        }
        return;
    }
    // not a leaf, so recurse
    GetAllIntersections(box->l_child, p, d, hit);
    GetAllIntersections(box->r_child, p, d, hit);
}

bool BVH::IntersectRecurse(AABB* box, Point3D p, Dir3D d, Hit &hit, int &checks) {
    // check if the box intersects
    float t_0, t_1(10000), t_2(10000);
    bool hit_box = RayIntersect(box, p ,d, t_0);
    if(!hit_box) {
        return false;
    }
    // check if this is a leaf
    if(box->l_child == nullptr && box->r_child == nullptr) {
        // check all the elements here
        return BoxIntersect(box, p, d, hit);
    }
    // now check the left and right boxes
    Hit hit_l, hit_r;
    bool hit1 = IntersectRecurse(box->l_child, p, d, hit_l, ++checks);
    bool hit2 = IntersectRecurse(box->r_child, p, d, hit_r, ++checks);
    if(hit1 || hit2) {
        hit = (hit_l.t < hit_r.t) ? hit_l : hit_r;
        return true;
    }
    return false;
    /*bool hit_l = RayIntersect(box->l_child, p ,d, t_1);
    bool hit_r = RayIntersect(box->r_child, p, d, t_2);
    // We recurse down the closest box first
    if (t_1 < t_2) {
        // go left first
        bool hit_elem = IntersectRecurse(box->l_child, p ,d, hit, ++checks);
        if(hit_elem) {
            return hit_elem;
        }
        if(hit_r) {
            hit_elem = IntersectRecurse(box->r_child, p, d, hit, ++checks);
        }
        return hit_elem;
    }
    else if(t_1 > t_2){
        // go right first
        bool hit_elem = IntersectRecurse(box->r_child, p ,d, hit, ++checks);
        if(hit_elem) {
            return hit_elem;
        }
        if(hit_l) {
            hit_elem = IntersectRecurse(box->l_child, p, d, hit, ++checks);
        }
        return hit_elem;
    }
    else {
        // we somehow hit both at the same time, so we return whichever is closer
        Hit l_hit, r_hit;
        hit_l = IntersectRecurse(box->l_child, p, d, l_hit, ++checks);
        hit_r = IntersectRecurse(box->r_child, p ,d, r_hit, ++checks);
        if(hit_l && hit_r) {
            hit = (l_hit.t < r_hit.t) ? l_hit : r_hit;
            return hit_l;
        }
        else if(hit_l) {
            hit = l_hit;
            return hit_l;
        }
        else if(hit_r) {
            hit = r_hit;
            return hit_r;
        }
    }
    // we should never get here, but return false to be safe
    return false;*/
}

// Check if the ray intersects the scene. if so, save the results in
// i_pt and t (intersection point and at time t)
bool BVH::BoxIntersect(AABB* box, Point3D p, Dir3D d, Hit &hit) {
  bool result = false;
  Hit cur_hit;
  float i_t(10000);
  int size = box->s_.size();
  for(int i = 0; i < size; i++) {
    // Find the closest ball, not the first ball
    if (sphereIntersect(p, d, box->s_[i], cur_hit) && cur_hit.t < i_t) {
      i_t = cur_hit.t;
      hit.pos = cur_hit.pos;
      hit.t = cur_hit.t;
      hit.m = box->s_[i].mat_;
      hit.norm = (hit.pos - box->s_[i].center_).normalized();
      result = true;
    }
  }
  size = box->t_.size();
  for(int i = 0; i < size; i++) {
    if(TriangleIntersect(p, d, *box->t_[i], cur_hit) && cur_hit.t < i_t) {
      hit.pos = cur_hit.pos;
      i_t = cur_hit.t;
      hit.t = cur_hit.t;
      hit.m = *(box->t_[i])->mat_;
      hit.norm = cur_hit.norm;
      result = true;
    }
  }
  return result;
}

bool BVH::RayIntersect(AABB* box, Point3D p, Dir3D d, float& i_t) {
    // intersection calculated using the slab method
    if(box == nullptr) {
        return false;
    }
    Dimension dim = box -> dim;
    float tmin(-INFINITY), tmax(INFINITY);
    //if(d.x != 0.0) {
        float tx1 = (dim.min_x - p.x)/d.x;
        float tx2 = (dim.max_x - p.x)/d.x;
        tmin = std::max(tmin, std::min(tx1, tx2));
        tmax = std::min(tmax, std::max(tx1, tx2));
    //}
    //if(d.y != 0.0) {
        float ty1 = (dim.min_y - p.y)/d.y;
        float ty2 = (dim.max_y - p.y)/d.y;
        tmin = std::max(tmin, std::min(ty1, ty2));
        tmax = std::min(tmax, std::max(ty1, ty2));
    //}
    //if(d.z != 0.0) {
        float tz1 = (dim.min_z - p.z)/d.z;
        float tz2 = (dim.max_z - p.z)/d.z;
        tmin = std::max(tmin, std::min(tz1, tz2));
        tmax = std::min(tmax, std::max(tz1, tz2));
    //}
    if(tmax > 0 && tmax >= tmin) {
        i_t = (tmin >= 0) ? tmin : tmax;
        return true;
    }
    return false;
}

// Check which side of the partition the element belongs in:
// -1 if left
// 0 if both
// 1 if right
int BVH::TriSide(float thresh, int axis, Triangle t) {
    // Check which sides the triangle pts fall on
    Point3D p = t.p1_;
    int p1_side = PointSide(thresh, axis, p);
    int p2_side = PointSide(thresh, axis, t.p2_);
    int p3_side = PointSide(thresh, axis, t.p3_);
    // if the sum is -3, then all points are on the left
    // if the sum is 3, then all points are on the right
    // otherwise it is on both sides
    int sum = p1_side + p2_side + p3_side;
    if (sum == 3) {
        return 1;
    }
    else if (sum == -3) {
        return -1;
    }
    else {
        return 0;
    }
}

int BVH::SphereSide(float thresh, int axis, Sphere s) {
    int cen_side = PointSide(thresh, axis, s.center_);
    if(cen_side == -1) {
        // check if the radius overlaps the line
        int max_side = PointSide(thresh, axis, s.center_ + Dir3D(s.radius_, s.radius_, s.radius_));
        return (max_side == -1) ? -1 : 0;
    }
    else if(cen_side == 1) {
        int min_side = PointSide(thresh, axis, s.center_ + Dir3D(-s.radius_, -s.radius_, -s.radius_));
        return (min_side == 1) ? 1 : 0;
    }
    return 0; 
}

int BVH::PointSide(float thresh, int axis, Point3D p) {
    float dif;
    if(axis == 0) {
        // x
        dif = p.x - thresh;
    }
    else if(axis == 1) {
        // y
        dif = p.y - thresh;
    }
    else if(axis == 2) {
        // z
        dif = p.z - thresh;
    }

    if(dif < 0.0) {
        return -1;
    }
    // (using fabs for precision)
    else if(fabs(dif) <= .00000001) {
        return 0;
    }
    else {
        return 1;
    }
}

// for a sphere to be within a box, we can wrap it in a box
// TODO - fix this!
bool BVH::contains(AABB* box, Sphere &s) {
    Dimension d = box -> dim;
    float r = s.radius_;

    Point3D min_corner = s.center_ + r*Dir3D(-1,-1,-1);
    Point3D max_corner = s.center_ + r*Dir3D(1, 1, 1);
    Dimension s_box;

    s_box.min_x = min_corner.x;
    s_box.min_y = min_corner.y;
    s_box.min_z = min_corner.z;

    s_box.max_x = max_corner.x;
    s_box.max_y = max_corner.y;
    s_box.max_z = max_corner.z;

    // if min greater than max or max less than min,
    // then no intersect
    if(contains(box, min_corner) || contains(box, max_corner)) {
        return true;
    }
    return false;
}

bool BVH::contains(AABB* box, Point3D &p) {
    Dimension d = box -> dim;
    bool in_x = d.min_x <= p.x && p.x <= d.max_x;
    bool in_y = d.min_y <= p.y && p.y <= d.max_y;
    bool in_z = d.min_z <= p.z && p.z <= d.max_z;
    if(in_x && in_y && in_z) {
        return true;
    }
    return false;
}
// for a triangle to be within a box, all the points must be inside it
bool BVH::contains(AABB* box, Triangle &t) {
    if(box == nullptr) {
        return false;
    }
    return contains(box, t.p1_) || contains(box, t.p2_) || contains(box, t.p3_);
}
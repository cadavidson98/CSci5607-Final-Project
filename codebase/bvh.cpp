#include "bvh.h"

#include <algorithm>
#include <stack>

using namespace std;

/**
 * Create a new bvh, that's it...
 */ 
bvh::bvh(vector<TriangleGL*> tris) {
    // start by making the bvh leaves
    triangles_ = tris;
    vector<triangle_info> leaves = boundTriangles();
    // now construct the bvh
    Dimension scene_dim = getExtent(leaves);
    compressed_node rt;
    bvh_nodes_.push_back(rt);
    buildRecurse(0, leaves);
    //validate(root);
}

/**
 * Create a compact version of the BVH which has all structs properly aligned for the GPU
 * The number of nodes can be found by using the sizeof() method
 */
nodeGL * bvh::getCompact(int &num_nodes) {
    num_nodes = bvh_nodes_.size();
    nodeGL* compact_bvh = new nodeGL[num_nodes];
    for (size_t i = 0; i < num_nodes; ++i) {
        nodeGL a_node;
        // set dimension values
        a_node.dim.max[0] = bvh_nodes_[i].AABB_.max_x;
        a_node.dim.max[1] = bvh_nodes_[i].AABB_.max_y;
        a_node.dim.max[2] = bvh_nodes_[i].AABB_.max_z;
        
        a_node.dim.min[0] = bvh_nodes_[i].AABB_.min_x;
        a_node.dim.min[1] = bvh_nodes_[i].AABB_.min_y;
        a_node.dim.min[2] = bvh_nodes_[i].AABB_.min_z;

        // set the node information
        a_node.l_child_offset = bvh_nodes_[i].l_child_;
        a_node.r_child_offset = bvh_nodes_[i].r_child_;
        if (bvh_nodes_[i].l_child_ == -1 && bvh_nodes_[i].l_child_ == -1) {
            a_node.triangle_offset = bvh_nodes_[i].tri_offsets_[0];
        }
        else {
            a_node.triangle_offset = -1;
        }
        compact_bvh[i] = a_node;
    }

    return compact_bvh;
}

void bvh::validate(node * cur_node) {
    if(cur_node == nullptr) return;
    if(cur_node ->l_child_ == nullptr && cur_node -> r_child_ == nullptr) {
        // this is a leaf
        for(size_t i = 0; i < cur_node->tri_offsets_.size(); ++i) {
            TriangleGL* tri = triangles_[cur_node->tri_offsets_[i]];
            printf("\tTriangle (%.2f, %2f, %.2f), (%.2f, %2f, %.2f), (%.2f, %2f, %.2f)\n",
            tri->p1[0], tri->p1[1], tri->p1[2],
            tri->p2[0], tri->p2[1], tri->p2[2],
            tri->p3[0], tri->p3[1], tri->p3[2]);
        }
        return;
    }
    validate(cur_node ->l_child_);
    validate(cur_node->r_child_);
}

/**
 * Construct a leaf node containing 1 singular triangle
 * This is trivially done by calculating the extent of the triangle
 */ 
vector<triangle_info> bvh::boundTriangles() {
    size_t num_tri = triangles_.size();
    vector<triangle_info> tri_nodes(num_tri);
    for(size_t i = 0; i < num_tri; ++i) {
        tri_nodes[i].tri_offset_ = i;
        
        Dimension tri_bnd;
        pair<float, float> x_bnds = minmax({triangles_[i]->p1[0], triangles_[i]->p2[0], triangles_[i]->p3[0]});
        pair<float, float> y_bnds = minmax({triangles_[i]->p1[1], triangles_[i]->p2[1], triangles_[i]->p3[1]});
        pair<float, float> z_bnds = minmax({triangles_[i]->p1[2], triangles_[i]->p2[2], triangles_[i]->p3[2]});
        
        tri_nodes[i].AABB_.min_x = x_bnds.first;
        tri_nodes[i].AABB_.max_x = x_bnds.second;
        
        tri_nodes[i].AABB_.min_y = y_bnds.first;
        tri_nodes[i].AABB_.max_y = y_bnds.second;
        
        tri_nodes[i].AABB_.min_z = z_bnds.first;
        tri_nodes[i].AABB_.max_z = z_bnds.second;

        tri_nodes[i].centroid_ = Point3D(.5*x_bnds.second+.5*x_bnds.first, .5*y_bnds.second+.5*y_bnds.first, .5*z_bnds.second+.5*z_bnds.first);
    }
    return tri_nodes;
}

/**
 * Determine the smallest boundary that encapsulates all the
 * triangles
 */ 
Dimension bvh::getExtent(vector<triangle_info> tris) {
    Dimension extent;
    size_t num_tris = tris.size();
    for(size_t i = 0; i < num_tris; ++i) {
        extent.max_x = max(extent.max_x, tris[i].AABB_.max_x);
        extent.max_y = max(extent.max_y, tris[i].AABB_.max_y);
        extent.max_z = max(extent.max_z, tris[i].AABB_.max_z);

        extent.min_x = min(extent.min_x, tris[i].AABB_.min_x);
        extent.min_y = min(extent.min_y, tris[i].AABB_.min_y);
        extent.min_z = min(extent.min_z, tris[i].AABB_.min_z);
    }
    return extent;
}

/**
 * Determine the smallest boundary that encapsulates all the
 * nodes.
 */ 
Dimension bvh::getExtent(vector<node*> nodes) {
    Dimension extent;
    size_t num_nodes = nodes.size();
    for(size_t i = 0; i < num_nodes; ++i) {
        extent.max_x = max(extent.max_x, nodes[i]->AABB_.max_x);
        extent.max_y = max(extent.max_y, nodes[i]->AABB_.max_y);
        extent.max_z = max(extent.max_z, nodes[i]->AABB_.max_z);

        extent.min_x = min(extent.min_x, nodes[i]->AABB_.min_x);
        extent.min_y = min(extent.min_y, nodes[i]->AABB_.min_y);
        extent.min_z = min(extent.min_z, nodes[i]->AABB_.min_z);
    }
    return extent;
}

/**
 * Determine the smallest boundary that encapsulates all the
 * points
 */ 
Dimension bvh::getExtent(vector<Point3D> pts) {
    Dimension extent;
    size_t num_pts = pts.size();
    for(size_t i = 0; i < num_pts; ++i) {
        extent.max_x = max(extent.max_x, pts[i].x);
        extent.max_y = max(extent.max_y, pts[i].y);
        extent.max_z = max(extent.max_z, pts[i].z);

        extent.min_x = min(extent.min_x, pts[i].x);
        extent.min_y = min(extent.min_y, pts[i].y);
        extent.min_z = min(extent.min_z, pts[i].z);
    }
    return extent;
}

/**
 * Calculate the surace area of the axis aligned bounding box 
 */
float bvh::surfaceArea(Dimension AABB) {
    return 
    2 * ((AABB.max_x - AABB.min_x) * (AABB.max_y - AABB.min_y)) +
    2 * ((AABB.max_x - AABB.min_x) * (AABB.max_z - AABB.min_z)) +
    2 * ((AABB.max_z - AABB.min_z) * (AABB.max_y - AABB.min_y));
}

/**
 * Construct a bvh by using the Surface Area Heuristic to subdivide
 * the leaf nodes provided 
 */
void bvh::buildRecurse(int node_offset, vector<triangle_info> tris) {
    if(tris.size() == 2) {
        compressed_node l_child, r_child;
        l_child.AABB_ = tris[0].AABB_;
        l_child.tri_offsets_.push_back(tris[0].tri_offset_);

        r_child.AABB_ = tris[1].AABB_;
        r_child.tri_offsets_.push_back(tris[1].tri_offset_);

        bvh_nodes_.push_back(l_child);
        bvh_nodes_.push_back(r_child);

        bvh_nodes_[node_offset].AABB_ = getExtent(tris);
        bvh_nodes_[node_offset].l_child_ = node_offset+1;
        bvh_nodes_[node_offset].r_child_ = node_offset+2;
        return;
    }
    if(tris.size() == 1) {
        compressed_node l_child;
        l_child.AABB_ = tris[0].AABB_;
        l_child.tri_offsets_.push_back(tris[0].tri_offset_);

        bvh_nodes_.push_back(l_child);

        bvh_nodes_[node_offset].AABB_ = getExtent(tris);
        bvh_nodes_[node_offset].l_child_ = node_offset+1;
        return;
    }
    bvh_nodes_[node_offset].AABB_ = getExtent(tris);
    vector<triangle_info> bin_1(0), bin_2(0);
    if(!splitMidpoint(tris, bin_1, bin_2)) {
        bvh_nodes_[node_offset].l_child_ = 
        bvh_nodes_[node_offset].r_child_ = -1;
        bvh_nodes_[node_offset].AABB_ = getExtent(tris);
        for(size_t i = 0; i < tris.size(); ++i) {
            bvh_nodes_[node_offset].tri_offsets_.push_back(tris[i].tri_offset_);
        }
        //cout << "created a leaf with " << tris.size() << " triangles" << endl;
        return;
    }
    // create left child
    compressed_node new_node;
    bvh_nodes_.push_back(new_node);
    int l_offset = node_offset + 1;
    bvh_nodes_[node_offset].l_child_ = l_offset;
    buildRecurse(l_offset, bin_1);
    // create right child
    bvh_nodes_.push_back(new_node);
    int r_offset = bvh_nodes_.size() - 1;
    bvh_nodes_[node_offset].r_child_ = r_offset;
    buildRecurse(r_offset, bin_2);
}

struct Bucket {
    Dimension bounds_;
    vector<triangle_info> tris_;
};

bool bvh::splitSAH(vector<triangle_info> in_tris, vector<triangle_info> &bin_1,
                   vector<triangle_info> &bin_2) {
    // bound the node
    Dimension scene_bnds = getExtent(in_tris);
    vector<Point3D> centroid_pts;
    for(size_t i = 0; i < in_tris.size(); ++i) {
        centroid_pts.push_back(in_tris[i].centroid_);
    }
    float min_sah = INFINITY;
    Dimension centroid_bnds = getExtent(centroid_pts);
    float range_x = centroid_bnds[0].second - centroid_bnds[0].first;
    float range_y = centroid_bnds[1].second - centroid_bnds[1].first;
    float range_z = centroid_bnds[2].second - centroid_bnds[2].first;
    float ranges[3] = { range_x, range_y, range_z };
    int axis = 0;
    // parition along semi-major axis, but using centroids to mitigate
    // large triangle issues...
    for(int i = 1; i < 3; i++) {
        if(ranges[i] > ranges[axis]) {
            axis = i;
        }
    }
    if(ranges[axis] < .0001) {
        return false;
    }
    float parent_sa = surfaceArea(scene_bnds);
    
    // try binning the triangles
    int num_bins = 16;
    Bucket bins[16];
        
    for(size_t i = 0; i < in_tris.size(); ++i) {
        // find out what bin this triangle belongs in
        // this can be done by manually finding the threshhold value
        // or just normalizing its value
        float range = ranges[axis];
        float centroid_pos = (axis == 0) ? in_tris[i].centroid_.x : (axis == 1) ? in_tris[i].centroid_.y : in_tris[i].centroid_.z;
        float normalized_position = (centroid_pos - centroid_bnds[axis].first) / range;
        int bucket_index = min(int(num_bins * normalized_position), num_bins - 1);
        bins[bucket_index].tris_.push_back(in_tris[i]);
    }

    // get the size of each bin
    for(int i = 0; i < num_bins; ++i) {
            bins[i].bounds_ = getExtent(bins[i].tris_);
    }
    float costs[15];
    // now get all the possible SAH values
    for(int i = 0; i < num_bins - 1; ++i) {
        int j = 0;
        Dimension bound1, bound2;
        int num_tris1(0), num_tris2(0);
        // this is one half split
        for(; j <= i; ++j) {
            bound1 = bound1.Union(bins[j].bounds_);
            num_tris1 += bins[j].tris_.size();
        }
        // this is the other half split
        for(; j < num_bins; ++j) {
            bound2 = bound2.Union(bins[j].bounds_);
            num_tris2 += bins[j].tris_.size();
        }
        // finally, calculate the SAH!
        costs[i] = .125 + (num_tris1 * surfaceArea(bound1) + num_tris2 * surfaceArea(bound2)) / parent_sa;
    }
    // now that we have all the costs, we can find the best split
    float min_cost = costs[0];
    int min_split = 0;
    for(int i = 0; i < num_bins - 1; ++i) {
        if(costs[i] < min_cost) {
            min_split = i;
        }
    }
    // Finally, determine whether or not we need to create a leaf node
    float parent_cost = .125 + in_tris.size();
    if(min_cost < parent_cost || in_tris.size() > max_tri) {
        // we need to subdivide
        float cen_min = centroid_bnds[axis].first;
        auto mid = partition(in_tris.begin(), in_tris.end(), [=](const triangle_info & tri) {
            float centroid_pos = (axis == 0) ? tri.centroid_.x : (axis == 1) ? tri.centroid_.y : tri.centroid_.z;
            float normalized_position = (centroid_pos - cen_min) / ranges[axis];
            int bucket_index = min(int(num_bins * normalized_position), num_bins - 1);
            return bucket_index <= min_split;
        });
        bin_1 = vector<triangle_info>(in_tris.begin(), mid);
        bin_2 = vector<triangle_info>(mid, in_tris.end());
        // we can continue to subdivide;
        return true;
    } else {
        // this needs to become a leaf node
        return false;
    }
}   

bool bvh::splitMidpoint(vector<triangle_info> in_tris, vector<triangle_info> &bin_1,
                        vector<triangle_info> &bin_2) {
    vector<Point3D> centroids;
    for (size_t i = 0; i < in_tris.size(); ++i) {
        centroids.push_back(in_tris[i].centroid_);
    }
    Dimension d = getExtent(centroids); 
    float thresh;
    float arr[3];
    arr[0] = abs(d.max_x - d.min_x);
    arr[1] = abs(d.max_y - d.min_y);
    arr[2] = abs(d.max_z - d.min_z);
    int index = 0;
    float max = arr[0];
    for(int i = 1; i < 3; ++i) {
        if(arr[i] > max) {
            index = i;
            max = arr[i];
        }
    }
    // now split
    if(index == 0) {
        // split on x
        std::sort(in_tris.begin(), in_tris.end(), [](const triangle_info& tri_1, const triangle_info &tri_2) {
            return tri_1.centroid_.x < tri_2.centroid_.x;
        });
    }
    else if (index == 1) {
        // split on y
        std::sort(in_tris.begin(), in_tris.end(), [](const triangle_info& tri_1, const triangle_info &tri_2) {
            return tri_1.centroid_.y < tri_2.centroid_.y;
        });
    }
    else {
        // split on z
        std::sort(in_tris.begin(), in_tris.end(), [](const triangle_info& tri_1, const triangle_info &tri_2) {
            return tri_1.centroid_.z < tri_2.centroid_.z;
        });
    }

    // now give half the leaves to each child
    int half = in_tris.size() / 2;
    bin_1 = vector<triangle_info>(in_tris.begin(), in_tris.begin()+half);
    bin_2 = vector<triangle_info>(in_tris.begin()+half, in_tris.end());
    return true;
}
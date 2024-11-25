#ifndef BVH_H
#define BVH_H

#include <vector>
#include <memory>
#include "sphere.h"
#include "cylinder.h"
#include "triangle.h"
#include "ray.h"
#include "shader_result.h"

struct AABB {
    std::vector<float> min_point{std::numeric_limits<float>::max(), 
                                std::numeric_limits<float>::max(), 
                                std::numeric_limits<float>::max()};
    std::vector<float> max_point{-std::numeric_limits<float>::max(), 
                                -std::numeric_limits<float>::max(), 
                                -std::numeric_limits<float>::max()};
    
    float surfaceArea() const {
        float dx = max_point[0] - min_point[0];
        float dy = max_point[1] - min_point[1];
        float dz = max_point[2] - min_point[2];
        return 2.0f * (dx * dy + dy * dz + dz * dx);
    }

    bool intersect(const Ray& ray, float& t_min) const;
    void expandToInclude(const AABB& other);
    void expandToInclude(const std::vector<float>& point);
};

class BVHNode {
public:
    AABB bounds;
    std::unique_ptr<BVHNode> left;
    std::unique_ptr<BVHNode> right;
    std::vector<Sphere> spheres;
    std::vector<Cylinder> cylinders;
    std::vector<Triangle> triangles;

    bool isLeaf() const {
        return !left && !right;
    }
};

struct BVHTraversal {
    float t;  // Distance to intersection
    BVHNode* node;
    
    BVHTraversal(float _t, BVHNode* _node) : t(_t), node(_node) {}
};

class BVH {
public:
    BVH(const std::vector<Sphere>& spheres,
        const std::vector<Cylinder>& cylinders,
        const std::vector<Triangle>& triangles);
        
    ShaderResult intersect(const Ray& ray) const;
    ShaderResult intersect(const Ray& ray, const std::vector<float>& backgroundcolor) const;
    bool isOccluded(const Ray& ray, float maxDistance) const;

private:
   std::unique_ptr<BVHNode> root;
    
    // Add these constants
    static constexpr int MIN_OBJECTS_PER_LEAF = 1;
    static constexpr int MAX_OBJECTS_PER_LEAF = 4;  // Reduced from 12
    static constexpr float TRAVERSAL_COST = 1.0f;
    static constexpr float INTERSECTION_COST = 2.0f;

    float calculateNodeCost(const AABB& bounds, int numObjects) const {
        float surface_area = bounds.surfaceArea();
        return surface_area * numObjects * INTERSECTION_COST;
    }

    std::unique_ptr<BVHNode> buildNode(std::vector<Sphere>& spheres,
                                      std::vector<Cylinder>& cylinders,
                                      std::vector<Triangle>& triangles);
                                      
    
};

#endif
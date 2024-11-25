#include "bvh.h"
#include "vector_utils.h"
#include <algorithm>
#include <queue>
#include <chrono>
#include <iostream>

bool AABB::intersect(const Ray& ray, float& t_min) const {
    float tmin = 0.0f;
    float tmax = std::numeric_limits<float>::infinity();

    for (int i = 0; i < 3; i++) {
        if (ray.direction[i] != 0.0f) {
            float t1 = (min_point[i] - ray.origin[i]) / ray.direction[i];
            float t2 = (max_point[i] - ray.origin[i]) / ray.direction[i];

            if (t1 > t2) std::swap(t1, t2);
            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);
            
            if (tmax <= tmin) return false;
        } else if (ray.origin[i] < min_point[i] || ray.origin[i] > max_point[i]) {
            return false;
        }
    }
    t_min = tmin;
    return true;
}

void AABB::expandToInclude(const AABB& other) {
    for (int i = 0; i < 3; i++) {
        min_point[i] = std::min(min_point[i], other.min_point[i]);
        max_point[i] = std::max(max_point[i], other.max_point[i]);
    }
}

void AABB::expandToInclude(const std::vector<float>& point) {
    for (int i = 0; i < 3; i++) {
        min_point[i] = std::min(min_point[i], point[i]);
        max_point[i] = std::max(max_point[i], point[i]);
    }
}

BVH::BVH(const std::vector<Sphere>& spheres, 
         const std::vector<Cylinder>& cylinders, 
         const std::vector<Triangle>& triangles) {
    std::vector<Sphere> spheres_copy = spheres;
    std::vector<Cylinder> cylinders_copy = cylinders;
    std::vector<Triangle> triangles_copy = triangles;
    root = buildNode(spheres_copy, cylinders_copy, triangles_copy);
}

std::unique_ptr<BVHNode> BVH::buildNode(std::vector<Sphere>& spheres,
                                       std::vector<Cylinder>& cylinders,
                                       std::vector<Triangle>& triangles) {
    auto node = std::make_unique<BVHNode>();
    
    // Calculate bounds
    if (!spheres.empty() || !cylinders.empty() || !triangles.empty()) {
        node->bounds = AABB();  // Initialize with infinity
        
        // Add sphere bounds
        for (const auto& sphere : spheres) {
            for (int i = 0; i < 3; i++) {
                node->bounds.min_point[i] = std::min(node->bounds.min_point[i], sphere.center[i] - sphere.radius);
                node->bounds.max_point[i] = std::max(node->bounds.max_point[i], sphere.center[i] + sphere.radius);
            }
        }
        
        // Add cylinder bounds
        for (const auto& cylinder : cylinders) {
            std::vector<float> top = {
                cylinder.center[0] + cylinder.height * cylinder.axis[0],
                cylinder.center[1] + cylinder.height * cylinder.axis[1],
                cylinder.center[2] + cylinder.height * cylinder.axis[2]
            };
            for (int i = 0; i < 3; i++) {
                node->bounds.min_point[i] = std::min(node->bounds.min_point[i], std::min(cylinder.center[i], top[i]) - cylinder.radius);
                node->bounds.max_point[i] = std::max(node->bounds.max_point[i], std::max(cylinder.center[i], top[i]) + cylinder.radius);
            }
        }
        
        // Add triangle bounds
        for (const auto& triangle : triangles) {
            for (int i = 0; i < 3; i++) {
                node->bounds.min_point[i] = std::min(node->bounds.min_point[i], std::min({triangle.v0[i], triangle.v1[i], triangle.v2[i]}));
                node->bounds.max_point[i] = std::max(node->bounds.max_point[i], std::max({triangle.v0[i], triangle.v1[i], triangle.v2[i]}));
            }
        }
    }
    
    // Leaf node criteria
    size_t total_objects = spheres.size() + cylinders.size() + triangles.size();
    if (total_objects <= 4) {  // Small leaf size for better splitting
        node->spheres = std::move(spheres);
        node->cylinders = std::move(cylinders);
        node->triangles = std::move(triangles);
        return node;
    }
    
    // Find the longest axis
    int axis = 0;
    float extent = node->bounds.max_point[0] - node->bounds.min_point[0];
    for (int i = 1; i < 3; i++) {
        float dim = node->bounds.max_point[i] - node->bounds.min_point[i];
        if (dim > extent) {
            extent = dim;
            axis = i;
        }
    }
    
    // Split point
    float split = node->bounds.min_point[axis] + extent * 0.5f;
    
    // Partition objects
    std::vector<Sphere> left_spheres, right_spheres;
    std::vector<Cylinder> left_cylinders, right_cylinders;
    std::vector<Triangle> left_triangles, right_triangles;
    
    // Split spheres
    for (auto& sphere : spheres) {
        if (sphere.center[axis] < split) {
            left_spheres.push_back(sphere);
        } else {
            right_spheres.push_back(sphere);
        }
    }
    
    // Split cylinders
    for (auto& cylinder : cylinders) {
        if (cylinder.center[axis] < split) {
            left_cylinders.push_back(cylinder);
        } else {
            right_cylinders.push_back(cylinder);
        }
    }
    
    // Split triangles
    for (auto& triangle : triangles) {
        float center = (triangle.v0[axis] + triangle.v1[axis] + triangle.v2[axis]) / 3.0f;
        if (center < split) {
            left_triangles.push_back(triangle);
        } else {
            right_triangles.push_back(triangle);
        }
    }
    
    // Handle degenerate splits
    // Handle degenerate splits
    if (left_spheres.empty() && left_cylinders.empty() && left_triangles.empty()) {
        // Move half of each type to the left
        size_t mid_spheres = spheres.size() / 2;
        size_t mid_cylinders = cylinders.size() / 2;
        size_t mid_triangles = triangles.size() / 2;
        
        left_spheres = std::vector<Sphere>(spheres.begin(), spheres.begin() + mid_spheres);
        right_spheres = std::vector<Sphere>(spheres.begin() + mid_spheres, spheres.end());
        
        left_cylinders = std::vector<Cylinder>(cylinders.begin(), cylinders.begin() + mid_cylinders);
        right_cylinders = std::vector<Cylinder>(cylinders.begin() + mid_cylinders, cylinders.end());
        
        left_triangles = std::vector<Triangle>(triangles.begin(), triangles.begin() + mid_triangles);
        right_triangles = std::vector<Triangle>(triangles.begin() + mid_triangles, triangles.end());
        
    } else if (right_spheres.empty() && right_cylinders.empty() && right_triangles.empty()) {
        // Move half of each type to the right
        size_t mid_spheres = left_spheres.size() / 2;
        size_t mid_cylinders = left_cylinders.size() / 2;
        size_t mid_triangles = left_triangles.size() / 2;
        
        right_spheres = std::vector<Sphere>(left_spheres.begin() + mid_spheres, left_spheres.end());
        left_spheres = std::vector<Sphere>(left_spheres.begin(), left_spheres.begin() + mid_spheres);
        
        right_cylinders = std::vector<Cylinder>(left_cylinders.begin() + mid_cylinders, left_cylinders.end());
        left_cylinders = std::vector<Cylinder>(left_cylinders.begin(), left_cylinders.begin() + mid_cylinders);
        
        right_triangles = std::vector<Triangle>(left_triangles.begin() + mid_triangles, left_triangles.end());
        left_triangles = std::vector<Triangle>(left_triangles.begin(), left_triangles.begin() + mid_triangles);
    }
    
    // Recursively build children
    node->left = buildNode(left_spheres, left_cylinders, left_triangles);
    node->right = buildNode(right_spheres, right_cylinders, right_triangles);
    
    return node;
}
// Define a comparator for the priority queue
struct BVHTraversalComparator {
    bool operator()(const BVHTraversal& a, const BVHTraversal& b) {
        return a.t > b.t;  // Min-heap based on distance
    }
};

ShaderResult BVH::intersect(const Ray& ray) const{
    std::vector<float> backgroundcolor = {0,0,0};
    return intersect(ray, backgroundcolor);
}

ShaderResult BVH::intersect(const Ray& ray, const std::vector<float>& backgroundcolor) const {

    ShaderResult result;
    result.color = backgroundcolor;
    result.intersected = false;
    
    struct StackEntry {
        BVHNode* node;
        float t_min;
    };
    
    StackEntry stack[64];  // Fixed-size stack
    int stack_ptr = 0;
    
    float t_min;
    if (root && root->bounds.intersect(ray, t_min)) {
        stack[stack_ptr++] = {root.get(), t_min};
    }
    
    float closest_t = std::numeric_limits<float>::max();
    
    while (stack_ptr > 0) {
        auto [node, t] = stack[--stack_ptr];
        
        if (t >= closest_t) continue;
        
        if (node->isLeaf()) {
            // Test all objects in leaf
            for (const auto& sphere : node->spheres) {
                float t;
                if (sphere.intersectSphere(ray, t) && t < closest_t) {
                    closest_t = t;
                    result.intersected = true;
                    result.intersected_material = sphere.material;
                    
                    result.intersection_point = {
                        ray.origin[0] + t * ray.direction[0],
                        ray.origin[1] + t * ray.direction[1],
                        ray.origin[2] + t * ray.direction[2]
                    };
                    
                    result.normal = {
                        result.intersection_point[0] - sphere.center[0],
                        result.intersection_point[1] - sphere.center[1],
                        result.intersection_point[2] - sphere.center[2]
                    };
                    normalize(result.normal);
                    
                    result.uv_coordinates = sphere.getUV(result.intersection_point);
                }
            }
            
            for (const auto& cylinder : node->cylinders) {
                float t;
                if (cylinder.intersectCylinder(ray, t) && t < closest_t) {
                    closest_t = t;
                    result.intersected = true;
                    result.intersected_material = cylinder.material;

                    result.intersection_point = {
                        ray.origin[0] + t * ray.direction[0],
                        ray.origin[1] + t * ray.direction[1],
                        ray.origin[2] + t * ray.direction[2]
                    };

                    std::vector<float> pc = {
                        result.intersection_point[0] - cylinder.center[0],
                        result.intersection_point[1] - cylinder.center[1],
                        result.intersection_point[2] - cylinder.center[2]
                    };
                    
                    
                    float dot = pc[0] * cylinder.axis[0] + pc[1] * cylinder.axis[1] + pc[2] * cylinder.axis[2];

                    result.normal = {
                        pc[0] - dot * cylinder.axis[0],
                        pc[1] - dot * cylinder.axis[1],
                        pc[2] - dot * cylinder.axis[2]
                    };
                    normalize(result.normal);

                    result.uv_coordinates = cylinder.getUV(result.intersection_point);
                }
            }

            for (const auto& triangle : node->triangles) {
                float t;
                if (triangle.intersectTriangle(ray, t) && t < closest_t) {
                    closest_t = t;
                    result.intersected = true;
                    result.intersected_material = triangle.material;

                    result.intersection_point = {
                        ray.origin[0] + t * ray.direction[0],
                        ray.origin[1] + t * ray.direction[1],
                        ray.origin[2] + t * ray.direction[2]
                    };

                    std::vector<float> edge1 = {triangle.v1[0] - triangle.v0[0],
                                                triangle.v1[1] - triangle.v0[1],
                                                triangle.v1[2] - triangle.v0[2]};

                    std::vector<float> edge2 = {triangle.v2[0] - triangle.v0[0],
                                                triangle.v2[1] - triangle.v0[1],
                                                triangle.v2[2] - triangle.v0[2]};

                    result.normal = {edge1[1] * edge2[2] - edge1[2] * edge2[1],
                                     edge1[2] * edge2[0] - edge1[0] * edge2[2],
                                     edge1[0] * edge2[1] - edge1[1] * edge2[0]};
                    normalize(result.normal);

                    result.uv_coordinates = triangle.getUV(result.intersection_point);
                }
            }
        } else {
            float t_left, t_right;
            bool hit_left = node->left->bounds.intersect(ray, t_left);
            bool hit_right = node->right->bounds.intersect(ray, t_right);
            
            if (hit_left && hit_right) {
                if (t_left > t_right) {
                    std::swap(node->left, node->right);
                    std::swap(t_left, t_right);
                }
                stack[stack_ptr++] = {node->right.get(), t_right};
                stack[stack_ptr++] = {node->left.get(), t_left};
            } else if (hit_left) {
                stack[stack_ptr++] = {node->left.get(), t_left};
            } else if (hit_right) {
                stack[stack_ptr++] = {node->right.get(), t_right};
            }
        }
    }

    return result;
}

bool BVH::isOccluded(const Ray& ray, float maxDistance) const {

    struct StackEntry {
        BVHNode* node;
        float t_min;
    };
    
    StackEntry stack[64];
    int stack_ptr = 0;

    float t_min;
    if (root && root->bounds.intersect(ray, t_min)) {
        stack[stack_ptr++] = {root.get(), t_min};
    }
    
    while (stack_ptr > 0) {
        auto [node, t] = stack[--stack_ptr];
        
        if (t >= maxDistance) continue;
        
        if (node->isLeaf()) {
            // Quick intersection tests for occlusion
            for (const auto& sphere : node->spheres) {
                float t;
                if (sphere.intersectSphere(ray, t) && t < maxDistance) {
                    return true;  // Found occlusion
                }
            }
            
            for (const auto& cylinder : node->cylinders) {
                float t;
                if (cylinder.intersectCylinder(ray, t) && t < maxDistance) {
                    return true;  // Found occlusion, exit early
                }
            }
            for (const auto& triangle : node->triangles) {
                float t;
                if (triangle.intersectTriangle(ray, t) && t < maxDistance) {
                    return true;  // Found occlusion, exit early
                }
            }
        } else {
            float t_left, t_right;
            bool hit_left = node->left->bounds.intersect(ray, t_left);
            bool hit_right = node->right->bounds.intersect(ray, t_right);
            
            if (hit_left && hit_right) {
                if (t_left > t_right) {
                    std::swap(node->left, node->right);
                    std::swap(t_left, t_right);
                }
                if (t_right < maxDistance) stack[stack_ptr++] = {node->right.get(), t_right};
                if (t_left < maxDistance) stack[stack_ptr++] = {node->left.get(), t_left};
            } else if (hit_left && t_left < maxDistance) {
                stack[stack_ptr++] = {node->left.get(), t_left};
            } else if (hit_right && t_right < maxDistance) {
                stack[stack_ptr++] = {node->right.get(), t_right};
            }
        }
    }

    return false;  // No occlusion found
}
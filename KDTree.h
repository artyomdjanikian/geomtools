#pragma once
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>

using Eigen::Vector3d;
using std::vector;
using std::unique_ptr;
using std::make_unique;

// KDTree for storing 3d points, supports insert point, remove point, radius search

class KDTree3D {
private:
    struct Node {
        Vector3d point;
        int axis;
        unique_ptr<Node> left;
        unique_ptr<Node> right;

        Node(const Vector3d& pt, int ax) : point(pt), axis(ax), left(nullptr), right(nullptr) {}
    };

    unique_ptr<Node> root;

    // Helper to build KD-Tree
    unique_ptr<Node> buildRecursive(vector<Vector3d>& pts, int depth) {
        if (pts.empty()) return nullptr;

        int axis = depth % 3;
        auto comparator = [axis](const Vector3d& a, const Vector3d& b) {
            return a[axis] < b[axis];
        };

        size_t medianIdx = pts.size() / 2;
        std::nth_element(pts.begin(), pts.begin() + medianIdx, pts.end(), comparator);
        Vector3d medianPoint = pts[medianIdx];

        vector<Vector3d> leftPts(pts.begin(), pts.begin() + medianIdx);
        vector<Vector3d> rightPts(pts.begin() + medianIdx + 1, pts.end());

        auto node = make_unique<Node>(medianPoint, axis);
        node->left = buildRecursive(leftPts, depth + 1);
        node->right = buildRecursive(rightPts, depth + 1);
        return node;
    }

    // Helper for insertion
    void insertRecursive(unique_ptr<Node>& node, const Vector3d& point, int depth) {
        if (!node) {
            node = make_unique<Node>(point, depth % 3);
            return;
        }

        int axis = node->axis;
        if (point[axis] < node->point[axis]) {
            insertRecursive(node->left, point, depth + 1);
        } else {
            insertRecursive(node->right, point, depth + 1);
        }
    }

    // Helper to find min in subtree
    Node* findMin(Node* node, int axis, int depth) {
        if (!node) return nullptr;

        int currentAxis = node->axis;
        if (currentAxis == axis) {
            if (node->left) return findMin(node->left.get(), axis, depth + 1);
            return node;
        }

        Node* leftMin = findMin(node->left.get(), axis, depth + 1);
        Node* rightMin = findMin(node->right.get(), axis, depth + 1);
        Node* minNode = node;

        if (leftMin && leftMin->point[axis] < minNode->point[axis]) minNode = leftMin;
        if (rightMin && rightMin->point[axis] < minNode->point[axis]) minNode = rightMin;

        return minNode;
    }

    // Helper for deletion
    unique_ptr<Node> deleteRecursive(unique_ptr<Node> node, const Vector3d& point, int depth) {
        if (!node) return nullptr;

        int axis = node->axis;
        if ((node->point - point).norm() < 1e-8) {
            if (node->right) {
                Node* minNode = findMin(node->right.get(), axis, depth + 1);
                node->point = minNode->point;
                node->right = deleteRecursive(std::move(node->right), minNode->point, depth + 1);
            } else if (node->left) {
                Node* minNode = findMin(node->left.get(), axis, depth + 1);
                node->point = minNode->point;
                node->right = deleteRecursive(std::move(node->left), minNode->point, depth + 1);
                node->left = nullptr;
            } else {
                return nullptr;
            }
            return node;
        }

        if (point[axis] < node->point[axis]) {
            node->left = deleteRecursive(std::move(node->left), point, depth + 1);
        } else {
            node->right = deleteRecursive(std::move(node->right), point, depth + 1);
        }
        return node;
    }

    // Helper for range search
    void searchRecursive(Node* node, const Vector3d& target, double eps, vector<Vector3d>& result) const {
        if (!node) return;

        if ((node->point - target).norm() <= eps) {
            result.push_back(node->point);
        }

        int axis = node->axis;
        double diff = target[axis] - node->point[axis];

        if (diff <= eps) searchRecursive(node->left.get(), target, eps, result);
        if (diff >= -eps) searchRecursive(node->right.get(), target, eps, result);
    }

public:
    KDTree3D() : root(nullptr) {}

    void build(const vector<Vector3d>& points) {
        vector<Vector3d> pts = points;
        root = buildRecursive(pts, 0);
    }

    void insert(const Vector3d& point) {
        insertRecursive(root, point, 0);
    }

    void remove(const Vector3d& point) {
        root = deleteRecursive(std::move(root), point, 0);
    }

    vector<Vector3d> radiusSearch(const Vector3d& query, double eps) const {
        vector<Vector3d> result;
        searchRecursive(root.get(), query, eps, result);
        return result;
    }
};

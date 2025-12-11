#pragma once
#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <tuple>
#include <eigen3/Eigen/Core>

// 3D integer cell coordinate
struct CellCoord {
    int x, y, z;

    bool operator==(const CellCoord& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// Custom hash for CellCoord for unordered_map
namespace std {
    template <>
    struct hash<CellCoord> {
        size_t operator()(const CellCoord& coord) const {
            size_t hx = hash<int>()(coord.x);
            size_t hy = hash<int>()(coord.y);
            size_t hz = hash<int>()(coord.z);
            return ((hx ^ (hy << 1)) >> 1) ^ (hz << 1);
        }
    };
}

// spatial hash grid for storing 3d points, supports radiusSearch, updatePoints
class SpatialHashGrid
{
public:
    explicit SpatialHashGrid(double cellSize_, const std::vector<Eigen::Vector3d> &points_) : cellSize(cellSize_), points(points_) {
        for(size_t iPnt = 0; iPnt < points.size(); iPnt++)
            insertOrUpdate(iPnt);
    }

    void updatePoints() {

        // very few updates, write shouldUpdate, if so insert iPnt into a vector and after the loop do all updates

        nUpdates = 0;
        for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {
            insertOrUpdate(iPnt);
        }
//        printf("%lu points, %d updates\n", points.size(), nUpdates);
    }

    size_t searchRadius(const Eigen::Vector3d &position, double radius, const std::function<bool(size_t)> &foundPoint) const
    {
        size_t nPnts = 0;
        int r = static_cast<int>(std::ceil(radius / cellSize));
        CellCoord center = getCellCoord(position);

        for (int dx = -r; dx <= r; ++dx)
        {
            for (int dy = -r; dy <= r; ++dy)
            {
                for (int dz = -r; dz <= r; ++dz)
                {
                    CellCoord neighbor{center.x + dx, center.y + dy, center.z + dz};
                    auto it = grid.find(neighbor);
                    if (it != grid.end())
                    {
                        for (auto iPnt : it->second)
                        {
                            if ((points[iPnt] - position).squaredNorm() <= radius * radius)
                            {
                                nPnts++;
                                if (foundPoint(iPnt))
                                    return nPnts;
                            }
                        }
                    }
                }
            }
        }
        return nPnts;
    }

    // TODO : push_back, size, operator [] for points vector access
    // operator [] should check if point is moved to a different CellCoord and update grid accordingly


private:
    int nUpdates = 0;

    void insertOrUpdate(size_t iPnt)
    {
        auto position = points[iPnt];
        CellCoord newCell = getCellCoord(position);
        auto it = pointCellMap.find(iPnt);

        if (it != pointCellMap.end())
        {
            CellCoord &oldCell = it->second;
            if (!(newCell == oldCell))
            {
                nUpdates++;
                removeFromCell(iPnt, oldCell);
                addToCell(iPnt, newCell);
                it->second = newCell;
            }
        }
        else
        {
            addToCell(iPnt, newCell);
            pointCellMap[iPnt] = newCell;
        }
    }

    void remove(size_t iPnt)
    {
        auto it = pointCellMap.find(iPnt);
        if (it != pointCellMap.end())
        {
            removeFromCell(iPnt, it->second);
            pointCellMap.erase(it);
        }
    }


private:
    double cellSize;
    const std::vector<Eigen::Vector3d> &points;

    std::unordered_map<CellCoord, std::vector<size_t>> grid;
    std::unordered_map<size_t, CellCoord> pointCellMap;

    CellCoord getCellCoord(const Eigen::Vector3d &pos) const
    {
        return { // cellSize should be large enough to avoid overflow   
            static_cast<int>(std::floor(pos.x() / cellSize)),
            static_cast<int>(std::floor(pos.y() / cellSize)),
            static_cast<int>(std::floor(pos.z() / cellSize))};
    }

    void addToCell(size_t iPnt, const CellCoord &cell)
    {
        grid[cell].push_back(iPnt);
    }

    void removeFromCell(size_t iPnt, const CellCoord &cell)
    {
        auto &vec = grid[cell];
        vec.erase(std::remove(vec.begin(), vec.end(), iPnt), vec.end());
        if (vec.empty())
        {
            grid.erase(cell);
        }
    }
};

#if 0

int main() {
    SpatialHashGrid grid(1.0f); // cell size = 1.0 unit

    Eigen::Vector3d p1(1, 0.1f, 0.2f, 0.3f);
    Eigen::Vector3d p2(2, 1.1f, 0.9f, 0.5f);
    Eigen::Vector3d p3(3, 3.0f, 3.0f, 3.0f);

    grid.insertOrUpdate(&p1);
    grid.insertOrUpdate(&p2);
    grid.insertOrUpdate(&p3);

    double queryX = 0.0f, queryY = 0.0f, queryZ = 0.0f, radius = 1.5f;
    auto nearby = grid.queryNearby(queryX, queryY, queryZ, radius);

    std::cout << "Nearby points:" << std::endl;
    for (auto p : nearby) {
        std::cout << "ID " << p->id << " at (" << p->x << ", " << p->y << ", " << p->z << ")" << std::endl;
    }

    return 0;
}

🛠️ Optional Enhancements:

    Use memory pools for performance in real-time simulations.

    Implement multithreaded insert/query for large simulations.

    Add support for multiple point types or tags (e.g., by layer).

    Let me know if you want GPU-accelerated versions or integration with physics engines (like Bullet or PhysX).
#endif

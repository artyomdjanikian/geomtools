#pragma once

#include "tools3d.h"
#include <vector>
#include <array>
#include <limits>


// octree implementation, stores triangles
// access through cOctreeIter

class cOctreeNode {

friend class cOctree;
friend class cOctreeIter;
public:

private:
    int m_id = -1;
    std::vector<uint32_t> m_facets;
    std::array<cOctreeNode *, 8> * m_subnodes = nullptr;
};

class cOctreeIter {
friend class cOctree;

private:
    cOctreeIter(cOctreeNode * pRoot, Bounds3d bounds) {
        sOctreeNodePath path;
        path.bounds = bounds;
        path.minCornerIndex = 0;
        path.edgeLen = (1<<19);
        path.pNode = pRoot;
        
        m_path.push_back(path);
    }

public:

    // access
    bool isValid() const { return (m_path.size() && m_path.back().pNode != nullptr); }
    bool isLeaf() const { return m_path.back().pNode->m_subnodes == nullptr; }
    bool hasChild(uint8_t iChild) const { return (*m_path.back().pNode->m_subnodes)[iChild] != nullptr; };
    std::vector<uint32_t> &getFacets() const { return m_path.back().pNode->m_facets; }

    // navigation
    size_t getDepth() const {return m_path.size();}

    bool down(uint8_t iChild) {
//        printf(">> %lu depth, down %d\n", m_path.size(), iChild);

        Eigen::Vector3d cornerPnt = getCorner(iChild);
        Eigen::Vector3d centerPnt = (m_path.back().bounds.minBound +  m_path.back().bounds.maxBound)* 0.5;

        Bounds3d bounds;
        bounds += cornerPnt;
        bounds += centerPnt;

        uint64_t edgeLen = m_path.back().edgeLen >> 1;
        if(edgeLen == 0) {
//            printf("    edge len %lu\n", edgeLen);
            return false;
        }

//        printf("%lu, %lu\n", m_path.back().minCornerIndex, m_path.back().edgeLen );

        uint64_t cornerIndex = getCornerIndex(iChild, m_path.back().minCornerIndex, edgeLen);
 
//        printf("    %lu, %lu\n", cornerIndex, edgeLen);

        sOctreeNodePath path;
        
        path.bounds = bounds;
        path.minCornerIndex = cornerIndex;
        path.edgeLen = edgeLen;
        path.pNode = isLeaf() ? nullptr : (*m_path.back().pNode->m_subnodes)[iChild];
        path.iChildIndex = 0;

//        printf("     %lu, %lu\n", m_path.back().pNode, m_path.back().pNode->m_subnodes);

        m_path.push_back(path);

//        printf("     %d depth, downed %d\n", m_path.size(), iChild);
              
//        printf("     %lu\n", m_path.back().pNode);
        return true;
    }

    bool up() {
//        printf("     %d depth, up\n", m_path.size());

        if(m_path.size()) {
            m_path.pop_back();
            return true;
        }

        return false;
    }

    void skip() {
        if(m_path.size())
            m_path.back().iChildIndex = 8;
    }

    const cOctreeIter & operator++() {
//        printf("++, %d depth\n", m_path.size());

        while(getDepth()) {
            if(isLeaf() || getChildIndex() == 8) {
                up();
                incrChildIndex();
            }
            else {
                if(hasChild(getChildIndex())) {
                    down(getChildIndex());
                    break;
                }
                else
                    incrChildIndex();
            }
        }

//        printf("++, %d new depth\n", m_path.size());

        return *this;
    };

    // geometry

    const Bounds3d & getBounds() const {
        if(!m_path.size())
            throw std::logic_error(" !!! invalid iterator !!!\n");

        return m_path.back().bounds;
    }

    std::array<Eigen::Vector3d, 4> getFaceCoords(uint8_t iFace) const {
        auto faceCorners = getFace(iFace);
        std::array<Eigen::Vector3d, 4> face;
        
        for(size_t iC = 0; iC < 4; iC++)
            face[iC] = getCorner(faceCorners[iC]);

        return face;
    }

    std::array<uint64_t, 4> getFaceIndices(uint8_t iFace) const {
        auto faceCorners = getFace(iFace);
        std::array<uint64_t, 4> face;
        
        for(size_t iC = 0; iC < 4; iC++)
            face[iC] = getCornerIndex(faceCorners[iC], m_path.back().minCornerIndex, m_path.back().edgeLen);

        return face;
    }

private:

    Eigen::Vector3d getCorner(uint8_t ii) const {

        uint8_t iX = ii&1;
        uint8_t iY = ii&2 >> 1;
        uint8_t iZ = ii&4 >> 2;

        const auto &bounds = m_path.back().bounds;

        Eigen::Vector3d cornerPnt = bounds.minBound;

        for(uint8_t iC = 0; iC < 3; iC++) {         
            uint8_t maskC = 1 << iC;

            if(ii & maskC)
                cornerPnt[iC] = bounds.maxBound[iC];
        }

        return cornerPnt;
    }

    uint64_t getCornerIndex(uint8_t index, uint64_t minCornerIndex, uint64_t edgeLen) const {

        // This is how (i, j, k) coordinates of tree node are stored as 64 bigetCornerIndext unsigned integer:
        /// 4 bits reserved, 20 bits for K, 20 bits for J, 20 bits for I

//        printf("corner %lu, edge len %lu, index %d\n", minCornerIndex, edgeLen, index);

        uint64_t indexIJK[3];

        uint64_t maskC = (1 << 20) - 1;
        uint64_t bitShift = 0;

        for(int iC = 0; iC < 3; iC++) {
            indexIJK[iC] = (minCornerIndex & maskC) >> bitShift;

            bitShift += 20;
            maskC <<= 20;
        }

        for(uint8_t iC = 0; iC < 3; iC++) {         
            uint8_t maskC = 1 << iC;

//            printf("    %d & %d = %d\n", index, maskC, index & maskC);

            if(index & maskC)
                indexIJK[iC] += edgeLen;
        }

//        printf("    >> %lu, %lu, %lu <<\n", indexIJK[0], indexIJK[1], indexIJK[2]);

        uint64_t cornerIndex = 0;

        for(uint8_t iC = 0; iC < 3; iC++)
            cornerIndex += indexIJK[iC] << (20*iC);

//        printf("corner index %lu\n", cornerIndex);

        return cornerIndex;
    }

    std::array<uint64_t, 4> getFace(uint8_t iFace) const {

        std::array<uint64_t, 4> face;

        switch(iFace) {
            case 0: 
                face = {0, 4, 6, 2}; //{2, 6, 4, 0};
                break;
            case 1:
                face = {1, 5, 4, 0}; //{0, 4, 5, 1};
                break;
            case 2:
                face = {2, 3, 1, 0}; //{0, 1, 3, 2};
                break;
            case 3:
                face = {3, 7, 5, 1}; //{1, 5, 7, 3};
                break;
            case 4:
                face = {2, 6, 7, 3}; //{3, 7, 6, 2};
                break;
            case 5:
                face = {7, 6, 4, 5}; //{5, 4, 6, 7};
                break;
            default:
            printf("   !!! wrong face index !!!\n");
        }

        return face;
    }
 
    cOctreeNode * getNode() {
        if(m_path.size())
            return m_path.back().pNode;

        return nullptr;
    }

    int getChildIndex() const {
        if(m_path.size())
            return m_path.back().iChildIndex;

        return -1;
    }

    void incrChildIndex() {
        if(m_path.size()) {
            m_path.back().iChildIndex++;
        }
    }

    private:

    struct sOctreeNodePath {
        Bounds3d bounds;
        uint64_t minCornerIndex = 0;
        uint64_t edgeLen = 0;

        cOctreeNode * pNode = nullptr;
        uint8_t iChildIndex = 0;
    };

    std::vector<sOctreeNodePath> m_path;
};

class cOctree 
{
public:

    void AddTriangle(Triangle3d triangle) {m_triangles.push_back(triangle);
                                           if(m_triangles.size() > std::numeric_limits<uint32_t>::max())
                                            throw std::logic_error(" !!! maximum number of triangles exceeded !!!"); }

    void Start() {
        for(const auto &triangle : m_triangles) {
            m_bounds += triangle.Pnt(0);
            m_bounds += triangle.Pnt(1);
            m_bounds += triangle.Pnt(2);
        }

        Eigen::Vector3d span = m_bounds.maxBound - m_bounds.minBound;
        m_len = std::max(span[0], span[1]);
        m_len = std::max(m_len, span[2]);

        Eigen::Vector3d diag(m_len, m_len, m_len);
        diag += m_bounds.minBound;
        m_bounds += diag;

        m_pRoot = new cOctreeNode();
        m_pRoot->m_id = m_nextId++;
        m_pRoot->m_subnodes = nullptr;

        for(uint32_t iF = 0; iF < m_triangles.size(); iF++)
            m_pRoot->m_facets.push_back(iF);
    }

    cOctreeIter begin() { return cOctreeIter(m_pRoot, m_bounds); }

    bool Subdivide(cOctreeIter &iter) {
//        printf("subdivide %d at depth %lu\n", iter.isValid(), iter.getDepth());

        cOctreeNode * pNode = iter.getNode();

        // get array of leafBounds;
        auto &facets = pNode->m_facets;

//        printf("subdivide %d facets\n", facets.size());
//        printf("   %d\n", iter.m_path.back().iChildIndex); 

        if(facets.size() == 0)
            return false;

        // allocate array of 8 pointers
        pNode->m_subnodes = new std::array<cOctreeNode*, 8>;
        for(int iN = 0; iN < 8; iN++)
            (*pNode->m_subnodes)[iN] = nullptr;

        std::array<Bounds3d, 8> leafBounds;
        std::array<std::vector<uint32_t>, 8> leafFacets;

        for(int iOct = 0; iOct < 8; iOct++) {
//            printf("down %d\n", iOct);
            iter.down(iOct);
            leafBounds[iOct] = iter.getBounds();
            iter.up();
        }

        for(int iLeaf = 0; iLeaf < 8; iLeaf++) {

            // leafBounds[iLeaf] to double boxcenter[3], double boxhalfsize[3]
//            printf(" %d, %.6f\n", iLeaf, (leafBounds[iLeaf].maxBound-leafBounds[iLeaf].minBound).norm());

            double boxCenter[3];
            double boxHalfsize[3];
            for(int iC = 0; iC < 3; iC++) {
                boxCenter[iC] = (leafBounds[iLeaf].minBound[iC]+leafBounds[iLeaf].maxBound[iC])*0.5;
                boxHalfsize[iC] = (leafBounds[iLeaf].maxBound[iC]-leafBounds[iLeaf].minBound[iC])*0.5;
            }

            for(int iFacet = 0; iFacet < static_cast<int>(facets.size()); iFacet++) {

                double triVerts[3][3];
                for(int iP = 0; iP < 3; iP++)
                    for(int iC = 0; iC < 3; iC++) {
                    // facet to double triverts[3][3]
                        triVerts[iP][iC] = m_triangles[facets[iFacet]].Pnt(iP)[iC];
                }

                int triXbox = triBoxOverlap(boxCenter, boxHalfsize, triVerts);
 
                if(triXbox)
                    leafFacets[iLeaf].push_back(facets[iFacet]);                    
            }
        }

        // allocate 8 nodes
        // link children to node

        int nAllocated = 0;
        for(int iLeaf = 0; iLeaf < 8; iLeaf++) {
//            printf("   %d. %d facets\n", iLeaf, leafFacets[iLeaf].size());
            if(leafFacets[iLeaf].size()) {
                // allocate leaf node
                nAllocated++;

                (*pNode->m_subnodes)[iLeaf] = new cOctreeNode();
                (*pNode->m_subnodes)[iLeaf]->m_id = m_nextId++;
                (*pNode->m_subnodes)[iLeaf]->m_subnodes = nullptr;
                (*pNode->m_subnodes)[iLeaf]->m_facets = std::move(leafFacets[iLeaf]);

//                printf("*** allocate %d, %d facets\n", (*pNode->m_subnodes)[iLeaf]->m_id, (*pNode->m_subnodes)[iLeaf]->m_facets.size());
            }
        }

        if(!nAllocated) {
//            printf("delete\n");
            delete pNode->m_subnodes;
            pNode->m_subnodes = nullptr;
            return false;
        }

        facets.clear();

        return true;
    }

public:

    std::vector<Triangle3d> m_triangles;

    int           m_nextId = 0;
    cOctreeNode * m_pRoot;
    Bounds3d      m_bounds;
    double        m_len;
};

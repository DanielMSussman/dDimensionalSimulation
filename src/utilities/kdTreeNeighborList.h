#ifndef kdTreeNeighborList_H
#define kdTreeNeighborList_H

#include "hyperrectangularCellList.h"
#include "baseNeighborList.h"
#include "kdTreeInterface.h"

/*! \file kdTreeNeighborList.h */
//!take a set of positions, sort those positions according to a cellList, and create data structures of possible neighbors of each particle
class kdTreeNeighborList : public baseNeighborList
    {
    public:
        //!basic constructor has a box and a range. pbc ONLY WORKS for hypercubes
        kdTreeNeighborList(scalar range, BoxPtr _box, double epsilon, bool pbc = true)
            {
            saveDistanceData = true;
            kdTree.radius = range;
            kdTree.epsilon = epsilon;
            dVec bds;
            _box->getBoxDims(bds);
            kdTree.boxLength  = bds.x[0];
            kdTree.periodicBoundary = pbc;
            Box = _box;
            }

        //!the interface class to CGAL
        kdTreeInterface kdTree;
    protected:

        BoxPtr Box;
        //!Save the displacement and distances associated with neihgbors?
        bool saveDistanceData;
        //! compute via GPU
        virtual void computeGPU(GPUArray<dVec> &points)
            {
            std::cout << "GPU routine not available for kdTree -- defaulting to CPU method" << std::endl;
            computeCPU(points);
            }
        //! compute via CPU
        virtual void computeCPU(GPUArray<dVec> &points);
    };

#endif


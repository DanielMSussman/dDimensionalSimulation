#ifndef neighborList_H
#define neighborList_H

#include "hyperrectangularCellList.h"
#include "baseNeighborList.h"
/*! \file neighborList.h */
//!take a set of positions, sort those positions according to a cellList, and create data structures of possible neighbors of each particle
class neighborList : public baseNeighborList
    {
    public:
        //!basic constructor has a box and a range
        neighborList(scalar range, BoxPtr _box, int subGridReduction = 1);

        //!The cell list that will help out
        shared_ptr<hyperrectangularCellList> cellList;
        //!Enforce GPU operation
        virtual void setGPU(bool _useGPU=true){
            useGPU = _useGPU;
            cellList->setGPU(useGPU);
            };
    protected:

        //!Save the displacement and distances associated with neihgbors?
        bool saveDistanceData;
        //!first index is Nmax, second is whether to recompute
        GPUArray<int> assist;
        //! compute via GPU
        virtual void computeGPU(GPUArray<dVec> &points);
        //! compute via CPU
        virtual void computeCPU(GPUArray<dVec> &points);
        //!Initialization and helper without using the GPU
        void resetNeighborsCPU(int size, int _nmax);
        //!Initialization and helper
        void resetNeighborsGPU(int size,int _nmax);

    };

#endif

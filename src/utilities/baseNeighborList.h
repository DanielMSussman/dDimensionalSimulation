#ifndef baseneighborList_H
#define baseneighborList_H

#include "kernelTuner.h"
/*! \file baseneighborList.h */
//!take a set of positions, sort those positions according to a cellList, and create data structures of possible neighbors of each particle
class baseNeighborList
    {
    public:
        //!basic constructor has a box and a range
        baseNeighborList(){};

        //!computethe neighborlist of the set of points passed in
        void computeNeighborLists(GPUArray<dVec> &points)
            {
            if(useGPU)
                {
                computeGPU(points);
                }
            else
                {
                computeCPU(points);
                }
            };

        void printNeighborInfo(GPUArray<dVec> &points, int idx)
            {
            ArrayHandle<dVec> p(points);
            cout << "particlePosition:" << endl;
            printdVec(p.data[idx]);
            ArrayHandle<unsigned int> npp(neighborsPerParticle);
            ArrayHandle<int> ns(particleIndices);
            cout <<"indices: " << endl;
            for (int i=0; i< npp.data[idx]; ++i)
                {
                int nIdx = ns.data[neighborIndexer(i,idx)];
                cout << nIdx << " ";
                printdVec(p.data[nIdx]);
                }
            cout << endl;
            };
        //!Enforce GPU operation
        virtual void setGPU(bool _useGPU=true){
            useGPU = _useGPU;
            };
        //!whether the updater does its work on the GPU or not
        bool useGPU;
        //!The Box used to compute periodic distances
        BoxPtr Box;

        void setBox(BoxPtr _bx){Box=_bx;};

        //!indexes the neighbors of each particle
        Index2D neighborIndexer;

        //! An array containing the number of elements in each neighborhood
        GPUArray<unsigned int> neighborsPerParticle;
        //!An array containing the indices of neighbors of each particle. So, partilceIndices[neighborIndexer(nn,pp)] gives the index of the nth particle in the neighborhood of particle pp
        GPUArray<int> particleIndices;
        //!An array saving the displacement data associated with each neighbor pair. distances[neighborIndexer(nn,pp)]
        GPUArray<dVec> neighborVectors;
        //!An array saving the distance data associated with each neighbor pair. distances[neighborIndexer(nn,pp)]
        GPUArray<scalar> neighborDistances;

        //!An internal counter
        int computations;

        //!maximum range that neighbors need to be kept at
        scalar maxRange;
        //! the maximum number of particles found in any neighborhood
        int Nmax;
        //!kernelTuner object
        shared_ptr<kernelTuner> nlistTuner;
    protected:

        //!Save the displacement and distances associated with neihgbors?
        bool saveDistanceData;
        //! compute via GPU
        virtual void computeGPU(GPUArray<dVec> &points){};
        //! compute via CPU
        virtual void computeCPU(GPUArray<dVec> &points){};
    };
#endif

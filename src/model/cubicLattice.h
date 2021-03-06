#ifndef cubicLattice_H
#define cubicLattice_H

#include "simpleModel.h"
#include "indexer.h"
#include "latticeBoundaries.h"

/*! \file cubicLattic.h
\brief puts degrees of freedom on a cubic lattice... probably for spin-like models
*/

//!define a type of simple model which places all degrees of freedom (which are still d-dimensional) on a cubic lattice with nearest neighbor interactions
class cubicLattice : public simpleModel
    {
    public:
        //!The base constructor takes the number of lattice sites along the cubic edge
        cubicLattice(int l, bool _slice = false,bool _useGPU = false);

        //!A rectilinear set of lattice sits
        cubicLattice(int lx, int ly, int lz, bool _slice = false,bool _useGPU = false);

        //!move the degrees of freedom
        virtual void moveParticles(GPUArray<dVec> &displacements,scalar scale = 1.);

        //!initialize each d.o.f. to be a unit spin on the sphere
        void setSpinsRandomly(noiseSource &noise);

        //! return the integer corresponding to the given site, along with the indices of the six nearest neighbors
        int getNeighbors(int target, vector<int> &neighbors, int &neighs);
        //!decide to slice sites
        void sliceIndices(bool _s=true){sliceSites = _s;};
        //!given a triple, determine what
        int latticeSiteToLinearIndex(const int3 &target);
        //!indexer for lattice sites
        Index3D latticeIndex;

        //!return the mean spin
        virtual dVec averagePosition()
            {
            dVec ans(0.0);
            ArrayHandle<dVec> spins(positions);
            ArrayHandle<int> t(types);
            int nSites=0;
            for(int i = 0; i < N; ++i)
                if(t.data[i] <= 0)
                    {
                    nSites+=1;
                    ans += spins.data[i];
                    };
            ans = (1.0/nSites)*ans;
            return ans;
        };

        //!list of the non-bulk objects in the simulations
        GPUArray<boundaryObject> boundaries;

        //!assign a collection of lattice sites to a new boundaryObject
        void createBoundaryObject(vector<int> &latticeSites, boundaryType _type, scalar Param1, scalar Param2);

    protected:
        //!a utility function for initialization
        void initializeNSites();
        //! should we use a memory-efficient slicing scheme?
        bool sliceSites;

        //!lattice sites per edge
        int L;

        //!normalize vector length when moving spins?
        bool normalizeSpins;
    };
#endif

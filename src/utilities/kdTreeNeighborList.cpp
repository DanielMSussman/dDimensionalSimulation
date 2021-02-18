#include "kdTreeNeighborList.h"

void kdTreeNeighborList::computeCPU(GPUArray<dVec> &points)
    {
    int Npoints = points.getNumElements();
    ArrayHandle<dVec> h_p(points,access_location::host,access_mode::read);
    kdTree.findNeighbors(Npoints,h_p.data);
    
    if(neighborsPerParticle.getNumElements() != Npoints)
        neighborsPerParticle.resize(Npoints);
    ArrayHandle<unsigned int> h_npp(neighborsPerParticle,access_location::host,access_mode::overwrite);
    for (int i = 0; i < Npoints; ++i)
        h_npp.data[i] = kdTree.allNeighbors[i].size();

    Nmax = kdTree.maximumNeighborNum;
    neighborIndexer = Index2D(Nmax,Npoints);
    if(particleIndices.getNumElements() != neighborIndexer.getNumElements())
        {
        particleIndices.resize(neighborIndexer.getNumElements());
        if(saveDistanceData)
            {
            neighborVectors.resize(neighborIndexer.getNumElements());
            neighborDistances.resize(neighborIndexer.getNumElements());
            };
        };

    ArrayHandle<int> h_idx(particleIndices,access_location::host,access_mode::overwrite);
    ArrayHandle<dVec> h_vec(neighborVectors,access_location::host,access_mode::overwrite);
    ArrayHandle<scalar> h_dist(neighborDistances,access_location::host,access_mode::overwrite);
    for (int i = 0; i < Npoints; ++i)
        {
        int nNeigh = kdTree.allNeighbors[i].size();
        dVec dist;
        for (int j = 0; j < nNeigh; ++j)
            {
            int nidx= kdTree.allNeighbors[i][j];
            h_idx.data[neighborIndexer(j,i)] = nidx;
            if(saveDistanceData)
                {
                Box->minDist(h_p.data[i],h_p.data[nidx],dist);
                h_vec.data[neighborIndexer(j,i)] = dist;
                h_dist.data[neighborIndexer(j,i)] = norm(dist);
                };
            }
        }

    };

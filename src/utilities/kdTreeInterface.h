#ifndef kdTreeInterface_h
#define kdTreeInterface_h

//! Access the d-dimensional range searching of CGAL
/*! A class for interfacing with the CGAL library's d-dimensional spatial searching functionality.
In particular, we make use of (potentially fuzzy) spherical range searching in the context of k-d tree neighbor searching
*/
class kdTreeInterface
    {
    public:
        kdTreeInterface(double _r=1.0, double _e=0.0, double _l=1.0, bool _pbc = false)
            {radius = _r; epsilon = _e; boxLength=_l, periodicBoundary = _pbc;};

        void findNeighbors(int N, dVec *points);

        double radius;
        double epsilon;
        double boxLength;
        bool periodicBoundary;
        std::vector< std::vector<int> > allNeighbors;
        bool binaryVecIterate(std::vector<int> &v);
    };
#endif

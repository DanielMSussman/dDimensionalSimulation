#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <utility>
#include "kdTreeInterface.h"

typedef CGAL::Epick_d<CGAL::Dimension_tag<DIMENSION> > K;
typedef K::Point_d Point_d;
typedef boost::tuple<Point_d,int>   Point_and_int;
typedef CGAL::Search_traits_d<K,CGAL::Dimension_tag<DIMENSION> >  Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int, CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
                                    Traits_base>    Traits;
typedef CGAL::Random_points_in_cube_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

/*! \file neighborList.cpp */

//!tool to iterate through a sequence of binary combinations. used to brute-force periodic neighbors
bool kdTreeInterface::binaryVecIterate(std::vector<int> &v)
    {
    int max = v.size();if(max <1) return false;
    v[0]+=1;
    int dd = 0;
    while(v[dd] > 1)
        {
        v[dd] = 0;
        dd+=1;
        if(dd == max) return false;
        v[dd] +=1;
        };
    return true;
    }

void kdTreeInterface::findNeighbors(int N, dVec *points)
    {
    std::vector<int> indices(N);
    std::vector<Point_d> pts(N);
    for (int ii = 0; ii < N; ++ii)
        {
        indices[ii]=ii;
        Point_d tempPoint(points[ii].x + 0,points[ii].x + DIMENSION);
        pts[ii] = tempPoint;
        }
    Tree tree(boost::make_zip_iterator(boost::make_tuple( pts.begin(),indices.begin() )),
              boost::make_zip_iterator(boost::make_tuple( pts.end(),indices.end() ) )
              );
    tree.build();

    //find neighbors
    Point_and_int periodicPoint;
    allNeighbors.resize(N);
    maximumNeighborNum = 0;
    for (int ii = 0; ii < N; ++ii)
        {
        std::vector<Point_and_int> pointNeighs; pointNeighs.reserve(30);
        Fuzzy_sphere sphereSearch(pts[ii],radius,epsilon);
        tree.search(back_inserter(pointNeighs),sphereSearch);

        //search periodic images if necessary
        if(periodicBoundary)
            {
            std::vector<int> nearbyBoundaries;
            for (int dd = 0; dd < DIMENSION;++dd)
                {
                if(boxLength-pts[ii][dd]<radius+epsilon)
                    nearbyBoundaries.push_back(dd+1);
                if(pts[ii][dd]<radius+epsilon)
                    nearbyBoundaries.push_back(-(dd+1));
                }

            std::vector<int> periodicCombinations(nearbyBoundaries.size(),0);
            while(binaryVecIterate(periodicCombinations))
                {
                double originalPoint[DIMENSION];
                std::copy(std::begin(points[ii].x), std::end(points[ii].x),std::begin(originalPoint));
                Point_and_int pbcPoint = pts[ii];
                for (int jj = 0; jj < nearbyBoundaries.size(); ++jj)
                    {
                    int coordinateIndex = abs(nearbyBoundaries[jj])-1;
                    double  sgn = (nearbyBoundaries[jj] > 0) ? 1. : -1.;
                    originalPoint[coordinateIndex] += periodicCombinations[jj]*sgn*boxLength;

                    }
                Point_d tempPoint(originalPoint + 0,originalPoint+ DIMENSION);
                boost::get<0>(pbcPoint) = tempPoint;
                Fuzzy_sphere sphereSearch2(pbcPoint,radius,epsilon);
                tree.search(back_inserter(pointNeighs),sphereSearch2);
                };
            };

        //insert all neighbor indicies into data structure
        allNeighbors[ii].clear();
        int nNeighs = pointNeighs.size();
        if(nNeighs > maximumNeighborNum) maximumNeighborNum = nNeighs;
        for (int jj = 0; jj < nNeighs; ++jj)
            {
            allNeighbors[ii].push_back(boost::get<1>(pointNeighs[jj]));
            }
        };
    };

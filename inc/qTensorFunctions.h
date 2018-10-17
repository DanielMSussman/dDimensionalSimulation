#ifndef QTENSORFUNCTIONS_H
#define QTENSORFUNCTIONS_H

#include "std_include.h"

/*!
The Q-tensor has five independent components, which will get passed around in dVec structures...
a dVec of q[0,1,2,3,4] corresponds to the symmetric traceless tensor laid out as
    (q[0]    q[1]        q[2]    )
Q = (q[1]    q[3]        q[4]    )
    (q[2]    q[4]   -(q[0]+q[3]) )

This file implements handy manipulations and functions of the Q-tensor as laid out this way
 */


#ifdef NVCC
#define HOSTDEVICE __host__ __device__ inline
#else
#define HOSTDEVICE inline __attribute__((always_inline))
#endif

/*! \file qTensorFunctions.h */

/** @defgroup Functions functions
 * @{
 \brief Utility functions that can be called from host or device
 */

//! return a qtensor given a director and a value of S0
HOSTDEVICE void qTensorFromDirector(scalar3 n, scalar S0, dVec &q)
    {
    //first, ensure that n is a unit vector
    scalar scale = sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
    if(scale == 0) // this shouldn't be called. ifit is, by default make \hat{[}n} the z-direction
        {
        q[0] = -0.5*S0; q[3] = -0.5*S0;
        }
    else if (scale != 1)
        {
        n.x /=scale;n.y /=scale;n.z /=scale;
        }
    q[0] = 1.5*(n.x*n.x-1.0/3.0)*S0;
    q[1] = 1.5*n.x*n.y;
    q[2] = 1.5*n.x*n.z;
    q[3] = 1.5*(n.y*n.y-1.0/3.0)*S0;
    q[4] = 1.5*n.y*n.z;
    };

//!Tr(Q^2) = Q_{kl}Q_{lk}
HOSTDEVICE scalar TrQ2(dVec &q)
    {
    return 2.0*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4] + q[0]*q[3]);
    };

//!Tr(Q^3)
HOSTDEVICE scalar TrQ3(dVec &q)
    {
    return -3.0*((-(q[1]*q[1]) + q[3]*q[3] + q[4]*q[4])*q[0] + q[0]*q[0]*q[3] - q[1]*q[1]*q[3] + q[2]*q[2]*q[3] - 2*q[1]*q[2]*q[4]);
    };

//!(Tr(Q^2))^2
HOSTDEVICE scalar TrQ2Squared(dVec &q)
    {
    scalar Q2 = TrQ2(q);
    return Q2*Q2;
    };

//!derivative of Tr(Q^2) w/r/t q[0] .. q[4]
HOSTDEVICE dVec derivativeTrQ2(dVec &q)
    {
    dVec ans;
    ans[0] = 2*(2*q[0] + q[3]);
    ans[1] = 4*q[1];
    ans[2] = 4*q[2];
    ans[3] = 2*(q[0] + 2*q[3]);
    ans[4] = 4*q[4];
    return ans;
    };

//!derivative of Tr(Q^3) w/r/t q[0] .. q[4]
HOSTDEVICE dVec derivativeTrQ3(dVec &q)
    {
    dVec ans;
    ans[0] = -3.0*(-q[1]*q[1] + q[3]*q[3] + q[4]*q[4] + 2*q[0]*q[3]);
    ans[1] = 6*(q[0]*q[1] + q[1]*q[3] + q[2]*q[4]);
    ans[2] = -6*q[2]*q[3] + 6*q[1]*q[4];
    ans[3] = -3*(q[0]*q[0] - q[1]*q[1] + q[2]*q[2] + 2*q[0]*q[3]);
    ans[4] = 6*q[1]*q[2] - 6*q[0]*q[4];
    return ans;
    };

//!derivative of (Tr(Q^2))^2 w/r/t q[0] .. q[4]
HOSTDEVICE dVec derivativeTrQ2Squared(dVec &q)
    {
    dVec ans;
    ans[0] = 8*(2*q[0] + q[3])*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4] + q[0]*q[3]);
    ans[1] = 16*q[1]*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4] + q[0]*q[3]);
    ans[2] = 16*q[2]*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4] + q[0]*q[3]);
    ans[3] = 8*(q[0] + 2*q[3])*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4] + q[0]*q[3]);
    ans[4] = 16*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4] + q[0]*q[3])*q[4];
    return ans;
    };

//!Q_{jk}Q_{ki}
HOSTDEVICE dVec QjkQki(dVec &q)
    {
    dVec ans;
    ans[0] = q[0]*q[0] + q[1]*q[1] + q[2]*q[2];
    ans[1] = q[0]*q[1] + q[1]*q[3] + q[2]*q[4];
    ans[2] = q[0]*q[2] - q[2]*(q[0]+q[3]) + q[1]*q[4];
    ans[3] = q[1]*q[1] + q[3]*q[3] + q[4]*q[4];
    ans[4] = q[1]*q[2] - q[4]*(q[0]+q[3]) + q[3]*q[4];
    return ans;
    };

//! determinant of a qt matrix
HOSTDEVICE scalar determinantOfQ(dVec &q)
    {
    return q[3]*(q[1]*q[1] - q[0]*q[0] - q[2]*q[2]) + 2*q[1]*q[2]*q[4] +
           q[0]*(q[1]*q[1] - q[3]*q[3] - q[4]*q[4]);
    };

//!Get the eigenvalues of a real symmetric traceless 3x3 matrix
HOSTDEVICE void eigenvaluesOfQ(dVec &q,scalar &a,scalar &b,scalar &c)
    {
    scalar p1 = q[1]*q[1] + q[2]*q[2]+ q[4]*q[4];
    if(p1 ==0) // diagonal matrix case
        {
        a=q[0];
        b=q[3];
        c=-q[0]-q[3];
        return;
        }
    //since q is traceless, some of these expressions are simpler than expected
    scalar p2 = 2*(q[0]*q[0]+q[3]*q[3] + q[0]*q[3]) + 2*p1;
    scalar p = sqrt(p2/6.0);
    dVec B = (1.0/p) * q;
    scalar r = 0.5*determinantOfQ(B);
    scalar phi;
    if(r <= -1)
        phi = PI/3.0;
    else if (r >= 1)
        phi = 0.0;
    else
        phi = acos(r)/3.0;

    a= 2.0*p*cos(phi);
    c= 2.0*p*cos(phi+2.0*PI/3.0);
    b=-a-c;
    }

/** @} */ //end of group declaration
#undef HOSTDEVICE
#endif

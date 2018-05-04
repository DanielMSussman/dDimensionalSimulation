#ifndef dDimensionalVectorTypes_h
#define dDimensionalVectorTypes_h
/*! \file dDimensionalVectorTypes.h 
defines dVec class (d-dimensional array of scalars)
defines iVec class (d-dimensional array of ints)
*/

//!dVec is an array whose length matches the dimension of the system
class dVec
    {
    public:
        HOSTDEVICE dVec(){};
        HOSTDEVICE dVec(const scalar value)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                x[dd] = value;
            };
        HOSTDEVICE dVec(const dVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                x[dd] = other.x[dd];
            };

        scalar x[DIMENSION];

        //mutating operators
        dVec& operator=(const dVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                this->x[dd] = other.x[dd];
            return *this;
            }
        dVec& operator-=(const dVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                this->x[dd] -= other.x[dd];
            return *this;
            }
        dVec& operator+=(const dVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                this->x[dd] += other.x[dd];
            return *this;
            }
    };

//!Less than operator for dVecs just sorts by the x-coordinate
HOSTDEVICE bool operator<(const dVec &a, const dVec &b)
    {
    return a.x[0]<b.x[0];
    }

//!Equality operator tests for.... equality of all elements
HOSTDEVICE bool operator==(const dVec &a, const dVec &b)
    {
    for (int dd = 0; dd <DIMENSION; ++dd)
        if(a.x[dd]!= b.x[dd]) return false;
    return true;
    }

//!return a dVec with all elements equal to one number
HOSTDEVICE dVec make_dVec(scalar value)
    {
    dVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = value;
    return ans;
    }

//!component-wise addition of two dVecs
HOSTDEVICE dVec operator+(const dVec &a, const dVec &b)
    {
    dVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a.x[dd]+b.x[dd];
    return ans;
    }

//!component-wise subtraction of two dVecs
HOSTDEVICE dVec operator-(const dVec &a, const dVec &b)
    {
    dVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a.x[dd]-b.x[dd];
    return ans;
    }

//!component-wise multiplication of two dVecs
HOSTDEVICE dVec operator*(const dVec &a, const dVec &b)
    {
    dVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a.x[dd]*b.x[dd];
    return ans;
    }

//!multiplication of dVec by scalar
HOSTDEVICE dVec operator*(const scalar &a, const dVec &b)
    {
    dVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a*b.x[dd];
    return ans;
    }

//!multiplication of dVec by scalar
HOSTDEVICE dVec operator*(const dVec &b, const scalar &a)
    {
    dVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a*b.x[dd];
    return ans;
    }

//!print a dVec to screen
inline __attribute__((always_inline)) void printdVec(dVec a)
    {
    cout <<"{";
    for (int dd = 0; dd < DIMENSION; ++dd)
        if(dd != DIMENSION-1)
            cout << a.x[dd] <<", ";
        else
            cout << a.x[dd];

    cout << "}" << endl;
    };


//!iVec is an array of ints whose length matches the dimension of the system
class iVec
    {
    public:
        HOSTDEVICE iVec(){};
        HOSTDEVICE iVec(const int value)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                x[dd] = value;
            };
        HOSTDEVICE iVec(const iVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                x[dd] = other.x[dd];
            };

        int x[DIMENSION];

        //mutating operators
        iVec& operator=(const iVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                this->x[dd] = other.x[dd];
            return *this;
            }
        iVec& operator-=(const iVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                this->x[dd] -= other.x[dd];
            return *this;
            }
        iVec& operator+=(const iVec &other)
            {
            for (int dd = 0; dd < DIMENSION; ++dd)
                this->x[dd] += other.x[dd];
            return *this;
            }
    };

//!Less than operator for dVecs just sorts by the x-coordinate
HOSTDEVICE bool operator<(const iVec &a, const iVec &b)
    {
    return a.x[0]<b.x[0];
    }

//!Equality operator tests for.... equality of all elements
HOSTDEVICE bool operator==(const iVec &a, const iVec &b)
    {
    for (int dd = 0; dd <DIMENSION; ++dd)
        if(a.x[dd]!= b.x[dd]) return false;
    return true;
    }

//!return a iVec with all elements equal to one number
HOSTDEVICE iVec make_dVec(int value)
    {
    iVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = value;
    return ans;
    }

//!component-wise addition of two iVecs
HOSTDEVICE iVec operator+(const iVec &a, const iVec &b)
    {
    iVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a.x[dd]+b.x[dd];
    return ans;
    }

//!component-wise subtraction of two iVecs
HOSTDEVICE iVec operator-(const iVec &a, const iVec &b)
    {
    iVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a.x[dd]-b.x[dd];
    return ans;
    }

//!component-wise multiplication of two iVecs
HOSTDEVICE iVec operator*(const iVec &a, const iVec &b)
    {
    iVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a.x[dd]*b.x[dd];
    return ans;
    }

//!multiplication of iVec by int
HOSTDEVICE iVec operator*(const int &a, const iVec &b)
    {
    iVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a*b.x[dd];
    return ans;
    }

//!multiplication of iVec by int
HOSTDEVICE iVec operator*(const iVec &b, const int &a)
    {
    iVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans.x[dd] = a*b.x[dd];
    return ans;
    }

//!modular addition of iVec (elementwise)
HOSTDEVICE iVec modularAddition(const iVec &i1, const iVec &i2, const iVec &max)
    {
    iVec ans;
    for (int dd = 0; dd < DIMENSION; ++dd)
        {
        ans.x[dd] = (i1.x[dd]+i2.x[dd])%max.x[dd];
        if(ans.x[dd] <0) ans.x[dd] += max.x[dd];
        };
    return ans;
    };

//!print a iVec to screen
inline __attribute__((always_inline)) void printiVec(iVec a)
    {
    cout <<"{";
    for (int dd = 0; dd < DIMENSION; ++dd)
        if(dd != DIMENSION-1)
            cout << a.x[dd] <<", ";
        else
            cout << a.x[dd];

    cout << "}" << endl;
    };

//! iterate through an iVec... on the first call, pass in (it = min except it.x[0] = min.x[0]-1
HOSTDEVICE bool iVecIterate(iVec &it, const iVec &min, const iVec &max)
    {
        it.x[0] +=1;
        int dd = 0;
        while(it.x[dd] >= max.x[dd]+1)
            {
            it.x[dd] = min.x[dd];
            dd +=1;
            if (dd == DIMENSION) return false;
            it.x[dd] += 1;
            };
        return true;
    };
#endif

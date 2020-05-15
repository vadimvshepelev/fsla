/*
-----------------------------------------------------------------------------
Adapted from Ogre3D
-----------------------------------------------------------------------------
*/

#ifndef __Vector4_H__
#define __Vector4_H__

#include <iostream>
#include <vector>
#include <assert.h>

using namespace std;

/** 4-dimensional homogeneous vector.
*/
class Vector4
{
public:

	double x, y, z, w;

    static const Vector4 ZERO;

    inline Vector4()
    {
    }

    inline Vector4( const double fX, const double fY, const double fZ, const double fW )
        : x( fX ), y( fY ), z( fZ ), w( fW)
    {
    }

    inline explicit Vector4( const double scalar )
        : x( scalar )
        , y( scalar )
        , z( scalar )
        , w( scalar )
    {
    }

	inline Vector4(vector<double> V) : x(V[0]), y(V[1]), z(V[2]), w(V[3]) {}

	inline double operator [] ( const size_t i ) const
    {
        assert( i < 4 );

        return *(&x+i);
    }

	inline double& operator [] ( const size_t i )
    {
        assert( i < 4 );

        return *(&x+i);
    }

    inline Vector4& operator = ( const Vector4& rkVector )
    {
        x = rkVector.x;
        y = rkVector.y;
        z = rkVector.z;
        w = rkVector.w;

        return *this;
    }

	inline Vector4& operator = ( const double fScalar)
	{
		x = fScalar;
		y = fScalar;
		z = fScalar;
		w = fScalar;
		return *this;
	}

    inline bool operator == ( const Vector4& rkVector ) const
    {
        return ( x == rkVector.x &&
            y == rkVector.y &&
            z == rkVector.z &&
            w == rkVector.w );
    }

    inline bool operator != ( const Vector4& rkVector ) const
    {
        return ( x != rkVector.x ||
            y != rkVector.y ||
            z != rkVector.z ||
            w != rkVector.w );
    }

    // arithmetic operations
    inline Vector4 operator + ( const Vector4& rkVector ) const
    {
        return Vector4(
            x + rkVector.x,
            y + rkVector.y,
            z + rkVector.z,
            w + rkVector.w);
    }

    inline Vector4 operator - ( const Vector4& rkVector ) const
    {
        return Vector4(
            x - rkVector.x,
            y - rkVector.y,
            z - rkVector.z,
            w - rkVector.w);
    }

    inline Vector4 operator * ( const double fScalar ) const
    {
        return Vector4(
            x * fScalar,
            y * fScalar,
            z * fScalar,
            w * fScalar);
    }

    inline Vector4 operator * ( const Vector4& rhs) const
    {
        return Vector4(
            rhs.x * x,
            rhs.y * y,
            rhs.z * z,
            rhs.w * w);
    }

    inline Vector4 operator / ( const double fScalar ) const
    {
        assert( fScalar != 0.0 );

        double fInv = 1.0 / fScalar;

        return Vector4(
            x * fInv,
            y * fInv,
            z * fInv,
            w * fInv);
    }

    inline Vector4 operator / ( const Vector4& rhs) const
    {
        return Vector4(
            x / rhs.x,
            y / rhs.y,
            z / rhs.z,
            w / rhs.w);
    }

    inline const Vector4& operator + () const
    {
        return *this;
    }

    inline Vector4 operator - () const
    {
        return Vector4(-x, -y, -z, -w);
    }

    inline friend Vector4 operator * ( const double fScalar, const Vector4& rkVector )
    {
        return Vector4(
            fScalar * rkVector.x,
            fScalar * rkVector.y,
            fScalar * rkVector.z,
            fScalar * rkVector.w);
    }

    inline friend Vector4 operator / ( const double fScalar, const Vector4& rkVector )
    {
        return Vector4(
            fScalar / rkVector.x,
            fScalar / rkVector.y,
            fScalar / rkVector.z,
            fScalar / rkVector.w);
    }

    inline friend Vector4 operator + (const Vector4& lhs, const double rhs)
    {
        return Vector4(
            lhs.x + rhs,
            lhs.y + rhs,
            lhs.z + rhs,
            lhs.w + rhs);
    }

    inline friend Vector4 operator + (const double lhs, const Vector4& rhs)
    {
        return Vector4(
            lhs + rhs.x,
            lhs + rhs.y,
            lhs + rhs.z,
            lhs + rhs.w);
    }

    inline friend Vector4 operator - (const Vector4& lhs, double rhs)
    {
        return Vector4(
            lhs.x - rhs,
            lhs.y - rhs,
            lhs.z - rhs,
            lhs.w - rhs);
    }

    inline friend Vector4 operator - (const double lhs, const Vector4& rhs)
    {
        return Vector4(
            lhs - rhs.x,
            lhs - rhs.y,
            lhs - rhs.z,
            lhs - rhs.w);
    }

    // arithmetic updates
    inline Vector4& operator += ( const Vector4& rkVector )
    {
        x += rkVector.x;
        y += rkVector.y;
        z += rkVector.z;
        w += rkVector.w;

        return *this;
    }

    inline Vector4& operator -= ( const Vector4& rkVector )
    {
        x -= rkVector.x;
        y -= rkVector.y;
        z -= rkVector.z;
        w -= rkVector.w;

        return *this;
    }

    inline Vector4& operator *= ( const double fScalar )
    {
        x *= fScalar;
        y *= fScalar;
        z *= fScalar;
        w *= fScalar;
        return *this;
    }

    inline Vector4& operator += ( const double fScalar )
    {
        x += fScalar;
        y += fScalar;
        z += fScalar;
        w += fScalar;
        return *this;
    }

    inline Vector4& operator -= ( const double fScalar )
    {
        x -= fScalar;
        y -= fScalar;
        z -= fScalar;
        w -= fScalar;
        return *this;
    }

    inline Vector4& operator *= ( const Vector4& rkVector )
    {
        x *= rkVector.x;
        y *= rkVector.y;
        z *= rkVector.z;
        w *= rkVector.w;

        return *this;
    }

    inline Vector4& operator /= ( const double fScalar )
    {
        assert( fScalar != 0.0 );

        double fInv = 1.0 / fScalar;

        x *= fInv;
        y *= fInv;
        z *= fInv;
        w *= fInv;

        return *this;
    }

    inline Vector4& operator /= ( const Vector4& rkVector )
    {
        x /= rkVector.x;
        y /= rkVector.y;
        z /= rkVector.z;
        w /= rkVector.w;

        return *this;
    }

    /** Calculates the dot (scalar) product of this vector with another.
    */
    inline double dotProduct(const Vector4& vec) const
    {
        return x * vec.x + y * vec.y + z * vec.z + w * vec.w;
	}

    /** Function for writing to a stream.
    */
	inline friend std::ostream& operator <<
		( std::ostream& o, const Vector4& v )
    {
        o << "Vector4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
        return o;
    }
};

#endif

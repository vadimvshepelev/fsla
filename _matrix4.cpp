/*
-----------------------------------------------------------------------------
Adapted from Ogre3D
-----------------------------------------------------------------------------
*/


#include "_matrix4.h"


const Matrix4 Matrix4::ZERO(
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0 );


const Matrix4 Matrix4::IDENTITY(
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1 );


inline static double
    MINOR(const Matrix4& m, const size_t r0, const size_t r1, const size_t r2, 
							const size_t c0, const size_t c1, const size_t c2)
{
    return m[r0][c0] * (m[r1][c1] * m[r2][c2] - m[r2][c1] * m[r1][c2]) -
        m[r0][c1] * (m[r1][c0] * m[r2][c2] - m[r2][c0] * m[r1][c2]) +
        m[r0][c2] * (m[r1][c0] * m[r2][c1] - m[r2][c0] * m[r1][c1]);
}


Matrix4 Matrix4::adjoint() const
{
    return Matrix4( MINOR(*this, 1, 2, 3, 1, 2, 3),
        -MINOR(*this, 0, 2, 3, 1, 2, 3),
        MINOR(*this, 0, 1, 3, 1, 2, 3),
        -MINOR(*this, 0, 1, 2, 1, 2, 3),

        -MINOR(*this, 1, 2, 3, 0, 2, 3),
        MINOR(*this, 0, 2, 3, 0, 2, 3),
        -MINOR(*this, 0, 1, 3, 0, 2, 3),
        MINOR(*this, 0, 1, 2, 0, 2, 3),

        MINOR(*this, 1, 2, 3, 0, 1, 3),
        -MINOR(*this, 0, 2, 3, 0, 1, 3),
        MINOR(*this, 0, 1, 3, 0, 1, 3),
        -MINOR(*this, 0, 1, 2, 0, 1, 3),

        -MINOR(*this, 1, 2, 3, 0, 1, 2),
        MINOR(*this, 0, 2, 3, 0, 1, 2),
        -MINOR(*this, 0, 1, 3, 0, 1, 2),
        MINOR(*this, 0, 1, 2, 0, 1, 2));
}


double Matrix4::determinant() const
{
    return m[0][0] * MINOR(*this, 1, 2, 3, 1, 2, 3) -
        m[0][1] * MINOR(*this, 1, 2, 3, 0, 2, 3) +
        m[0][2] * MINOR(*this, 1, 2, 3, 0, 1, 3) -
        m[0][3] * MINOR(*this, 1, 2, 3, 0, 1, 2);
}


Matrix4 Matrix4::inverse() const
{
    double m00 = m[0][0], m01 = m[0][1], m02 = m[0][2], m03 = m[0][3];
    double m10 = m[1][0], m11 = m[1][1], m12 = m[1][2], m13 = m[1][3];
    double m20 = m[2][0], m21 = m[2][1], m22 = m[2][2], m23 = m[2][3];
    double m30 = m[3][0], m31 = m[3][1], m32 = m[3][2], m33 = m[3][3];

    double v0 = m20 * m31 - m21 * m30;
    double v1 = m20 * m32 - m22 * m30;
    double v2 = m20 * m33 - m23 * m30;
    double v3 = m21 * m32 - m22 * m31;
    double v4 = m21 * m33 - m23 * m31;
    double v5 = m22 * m33 - m23 * m32;

    double t00 = + (v5 * m11 - v4 * m12 + v3 * m13);
    double t10 = - (v5 * m10 - v2 * m12 + v1 * m13);
    double t20 = + (v4 * m10 - v2 * m11 + v0 * m13);
    double t30 = - (v3 * m10 - v1 * m11 + v0 * m12);

    double invDet = 1 / (t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03);

    double d00 = t00 * invDet;
    double d10 = t10 * invDet;
    double d20 = t20 * invDet;
    double d30 = t30 * invDet;

    double d01 = - (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
    double d11 = + (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
    double d21 = - (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
    double d31 = + (v3 * m00 - v1 * m01 + v0 * m02) * invDet;

    v0 = m10 * m31 - m11 * m30;
    v1 = m10 * m32 - m12 * m30;
    v2 = m10 * m33 - m13 * m30;
    v3 = m11 * m32 - m12 * m31;
    v4 = m11 * m33 - m13 * m31;
    v5 = m12 * m33 - m13 * m32;

    double d02 = + (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
    double d12 = - (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
    double d22 = + (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
    double d32 = - (v3 * m00 - v1 * m01 + v0 * m02) * invDet;

    v0 = m21 * m10 - m20 * m11;
    v1 = m22 * m10 - m20 * m12;
    v2 = m23 * m10 - m20 * m13;
    v3 = m22 * m11 - m21 * m12;
    v4 = m23 * m11 - m21 * m13;
    v5 = m23 * m12 - m22 * m13;

    double d03 = - (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
    double d13 = + (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
    double d23 = - (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
    double d33 = + (v3 * m00 - v1 * m01 + v0 * m02) * invDet;

    return Matrix4(
        d00, d01, d02, d03,
        d10, d11, d12, d13,
        d20, d21, d22, d23,
        d30, d31, d32, d33);
}

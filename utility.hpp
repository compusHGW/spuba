#ifndef PROJECT_UTITLITY_HPP
#define PROJECT_UTITLITY_HPP

#include <cmath>

inline void norm(double *a) {
    double temp = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= temp;
    a[1] /= temp;
    a[2] /= temp;
    //printf("%f %f %f\n",a[0],a[1],a[2]);
}

inline double _len(const double *a) {
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

inline double _dot(const double *a, const double *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double Angle(const double *a, const double *b) {
    return acos(_dot(a, b) / sqrt(_dot(a, a)) * sqrt(_dot(b, b)));
}

inline void _vec(const double *a, const double *b, double *res) {
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}

inline void _add(const double *a, const double *b, double *res) {
    res[0] = a[0] + b[0];
    res[1] = a[1] + b[1];
    res[2] = a[2] + b[2];
}


inline void PolarToCartesian(const double *polar, double *cartesian) {
    cartesian[0] = polar[2] * sin(polar[0]) * cos(polar[1]);
    cartesian[1] = polar[2] * sin(polar[0]) * sin(polar[1]);
    cartesian[2] = polar[2] * cos(polar[0]);
}

inline void MultMatrix(double matrix[][3], double mult[][3]) {
    register double temp[3];
    temp[0] = matrix[0][0] * mult[0][0] + matrix[0][1] * mult[1][0] + matrix[0][2] * mult[2][0];
    temp[1] = matrix[0][0] * mult[0][1] + matrix[0][1] * mult[1][1] + matrix[0][2] * mult[2][1];
    temp[2] = matrix[0][0] * mult[0][2] + matrix[0][1] * mult[1][2] + matrix[0][2] * mult[2][2];
    matrix[0][0] = temp[0];
    matrix[0][1] = temp[1];
    matrix[0][2] = temp[2];

    temp[0] = matrix[1][0] * mult[0][0] + matrix[1][1] * mult[1][0] + matrix[1][2] * mult[2][0];
    temp[1] = matrix[1][0] * mult[0][1] + matrix[1][1] * mult[1][1] + matrix[1][2] * mult[2][1];
    temp[2] = matrix[1][0] * mult[0][2] + matrix[1][1] * mult[1][2] + matrix[1][2] * mult[2][2];
    matrix[1][0] = temp[0];
    matrix[1][1] = temp[1];
    matrix[1][2] = temp[2];

    temp[0] = matrix[2][0] * mult[0][0] + matrix[2][1] * mult[1][0] + matrix[2][2] * mult[2][0];
    temp[1] = matrix[2][0] * mult[0][1] + matrix[2][1] * mult[1][1] + matrix[2][2] * mult[2][1];
    temp[2] = matrix[2][0] * mult[0][2] + matrix[2][1] * mult[1][2] + matrix[2][2] * mult[2][2];
    matrix[2][0] = temp[0];
    matrix[2][1] = temp[1];
    matrix[2][2] = temp[2];
}

inline void RotateZ(double matrix[][3], double matrixOut[][3], double alpha) {
    register double s, c;
    sincos(alpha, &s, &c);
    matrix[0][0] = c;
    matrix[0][1] = -s;
    matrix[0][2] = 0.0;
    matrix[1][0] = s;
    matrix[1][1] = c;
    matrix[1][2] = 0.0;
    matrix[2][0] = 0.0;
    matrix[2][1] = 0.0;
    matrix[2][2] = 1.0;
    MultMatrix(matrixOut, matrix);
}

inline void RotateX(double matrix[][3], double matrixOut[][3], double alpha) {
    register double s, c;
    sincos(alpha, &s, &c);
    matrix[0][0] = 1.0;
    matrix[0][1] = 0.0;
    matrix[0][2] = 0.0;
    matrix[1][0] = 0.0;
    matrix[1][1] = c;
    matrix[1][2] = -s;
    matrix[2][0] = 0.0;
    matrix[2][1] = s;
    matrix[2][2] = c;
    MultMatrix(matrixOut, matrix);
}

inline void RotateY(double matrix[][3], double matrixOut[][3], double alpha) {
    register double s, c;
    sincos(alpha, &s, &c);
    matrix[0][0] = c;
    matrix[0][1] = 0.0;
    matrix[0][2] = s;
    matrix[1][0] = 0.0;
    matrix[1][1] = 1.0;
    matrix[1][2] = 0.0;
    matrix[2][0] = -s;
    matrix[2][1] = 0.0;
    matrix[2][2] = c;
    MultMatrix(matrixOut, matrix);
}



inline void CartesianToPolar(const double *cartesian, double *polar) {
    //(theta,phi,r)
    polar[2] = sqrt(cartesian[0] * cartesian[0] + cartesian[1] * cartesian[1] + cartesian[2] * cartesian[2]);
    polar[0] = acos(cartesian[2] / polar[2]);
    polar[1] = atan2(cartesian[1], cartesian[0]);
}

inline void CartesianToCylinder(const double *cartesian, double *cylin) {
    //(r,phi,z)
    cylin[0] = sqrt(cartesian[0] * cartesian[0] + cartesian[1] * cartesian[1]);
    cylin[2] = cartesian[2];
    if ((cartesian[0] == 0) & (cartesian[1] == 0)) cylin[1] = 0.0;
    else if (cartesian[0] >= 0) cylin[1] = asin(cartesian[1] / cylin[0]);
    else cylin[1] = -asin(cartesian[1] / cylin[0]) + M_PI;
}

inline void set_identity(double matrix[3][3]) {
    matrix[0][0] = 1.0;
    matrix[0][1] = 0.0;
    matrix[0][2] = 0.0;
    matrix[1][0] = 0.0;
    matrix[1][1] = 1.0;
    matrix[1][2] = 0.0;
    matrix[2][0] = 0.0;
    matrix[2][1] = 0.0;
    matrix[2][2] = 1.0;
}

#endif //PROJECT_UTITLITY_HPP

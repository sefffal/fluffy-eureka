/*
 * vectors.c
 *
 * Library for fast 3D vector operations.x
 *
 * Author : Will Thompson
 * Date   : 2016-09-11
 *
 * Implementation File
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vectors.h"  /* Include the header file for this library from the current director */

/* Vector3
 * Add two Vector3's components together, and return the result as a new Vector3.
 */
Vector3 vadd(Vector3 a, Vector3 b) {
    Vector3 out;
    out.x = a.x + b.x;
    out.y = a.y + b.y;
    out.z = a.z + b.z;
    assert(out.x == a.x + b.x && out.y == a.y + b.y && out.z == a.z + b.z);
    return out;
}

/* Vector3
 * Multiply two Vector3's components together, and return the result as a new Vector3.
 */
Vector3 vmul(Vector3 a, Vector3 b) {
    Vector3 out;
    out.x = a.x * b.x;
    out.y = a.y * b.y;
    out.z = a.z * b.z;
    assert(out.x == a.x * b.x && out.y == a.y * b.y && out.z == a.z * b.z);    
    return out;
}

/* Vector3
 * Calculate the magnitude/length of a Vector3 and return the result.
 */
double vmag(Vector3 a) {
    double out = sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2));
    assert(out >= 0.0);
    return out;
}
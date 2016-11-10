/*
 * vectors.c
 *
 * Library for fast 3D vector operations.x
 *
 * Author : Will Thompson
 * Date   : 2016-09-11
 *
 * Header file.
 */

/* This jiggery-pokery is to prevent the contents from 
 * being included multiple times if it is <included>
 * more than once (very possible)
 */
#ifndef VECTORS_H
# define VECTORS_H


/* Vector3
 * Struct representing a 3D vector of some physical quantity.
 */
typedef struct {
    double x;
    double y;
    double z;
} Vector3;

Vector3 vadd(Vector3, Vector3);
Vector3 vmul(Vector3, Vector3);
double vmag(Vector3);

/* End jiggery-pokery */
#endif
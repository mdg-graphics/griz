/* $Id$ */
/*
 * geometric.h - Vectors, transformation matrices, transforms.
 *
 *      modified by Donald J. Dovey
 *      Lawrence Livermore National Laboratory
 *      Jan 28 1992
 *
 */


/**
 ** sipp - SImple Polygon Processor
 **
 **  A general 3d graphic package
 **
 **  Copyright Jonas Yngvesson  (jonas-y@isy.liu.se) 1988/89/90/91
 **            Inge Wallin      (ingwa@isy.liu.se)         1990/91
 **
 ** This program is free software; you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 1, or any later version.
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU General Public License for more details.
 ** You can receive a copy of the GNU General Public License from the
 ** Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 **/


#ifndef GEOMETRIC_H
#define GEOMETRIC_H

/*****************************************************************
 * TAG( X_AXIS Y_AXIS Z_AXIS )
 *
 * Some definitions to keep the code readable.
 */
#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2


/*****************************************************************
 * TAG( VEC_SET VEC_DOT VEC_LENGTH VEC_COPY VEC_NEGATE VEC_ADD )
 * TAG( VEC_SUB VEC_SCALE VEC_CROSS VEC_ADDS VEC_AVER )
 *
 * Vector operations.
 * Capitalized types denote Vectors and other aggregates.  Lower
 * case denote scalars.
 */

#define VEC_SET(A, xx, yy, zz)   { (A)[0] = (xx); \
                                   (A)[1] = (yy); \
                                   (A)[2] = (zz); }
#define VEC_DOT(A, B)      ((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2])
#define VEC_DOT_2D(A, B)   ((A)[0]*(B)[0] + (A)[1]*(B)[1])
#define VEC_LENGTH(A)      ((float)sqrt((double)VEC_DOT(A, A)))
#define VEC_LENGTH_2D(A)   ((float)sqrt((double)VEC_DOT_2D(A, A)))
#define VEC_COPY(B, A)     { (B)[0] = (A)[0]; \
                             (B)[1] = (A)[1]; \
                             (B)[2] = (A)[2]; }
#define VEC_NEGATE(A)      { (A)[0] = -(A)[0]; \
                             (A)[1] = -(A)[1]; \
                             (A)[2] = -(A)[2]; }
#define VEC_ADD(C, A, B)   { (C)[0] = (A)[0] + (B)[0]; \
                             (C)[1] = (A)[1] + (B)[1]; \
                             (C)[2] = (A)[2] + (B)[2]; }
#define VEC_SUB(C, A, B)   { (C)[0] = (A)[0] - (B)[0]; \
                             (C)[1] = (A)[1] - (B)[1]; \
                             (C)[2] = (A)[2] - (B)[2]; }
#define VEC_SUB_2D(C, A, B) { (C)[0] = (A)[0] - (B)[0]; \
                              (C)[1] = (A)[1] - (B)[1]; }
#define VEC_SCALE(C, s, A)  { (C)[0] = (s) * (A)[0]; \
                              (C)[1] = (s) * (A)[1]; \
                              (C)[2] = (s) * (A)[2]; }
#define VEC_CROSS(C, A, B)   { (C)[0] = (A)[1]*(B)[2] - (A)[2]*(B)[1]; \
                               (C)[1] = (A)[2]*(B)[0] - (A)[0]*(B)[2]; \
                               (C)[2] = (A)[0]*(B)[1] - (A)[1]*(B)[0]; }
#define VEC_ADDS(C, a, A, B)	 { (C)[0] = (a)*(A)[0] + (B)[0]; \
				   (C)[1] = (a)*(A)[1] + (B)[1]; \
				   (C)[2] = (a)*(A)[2] + (B)[2]; }
#define VEC_AVER(C, a, A, b, B)	 { (C)[0] = (a)*(A)[0] + (b)*(B)[0]; \
				   (C)[1] = (a)*(A)[1] + (b)*(B)[1]; \
			 	   (C)[2] = (a)*(A)[2] + (b)*(B)[2]; }


/*****************************************************************
 * TAG( L_INTERP )
 *
 * Linear interpolation between two points.  Notice that t = 0
 * gives point A and t = 1 gives point B.
 */
#define L_INTERP(C, t, A, B)	 { (C)[0] = (1.0-(t))*(A)[0] + (t)*(B)[0]; \
				   (C)[1] = (1.0-(t))*(A)[1] + (t)*(B)[1]; \
			 	   (C)[2] = (1.0-(t))*(A)[2] + (t)*(B)[2]; }


/*****************************************************************
 * TAG( vec_norm )
 *
 * Normalize a vector.
 */
extern void vec_norm( /* Vector */ );


/*****************************************************************
 * TAG( norm_three_pts )
 *
 * Compute a normal vector to a plane from three points in the
 * plane.  norm = cross( P3 - P2, P1 - P2 ).
 */
extern void norm_three_pts( /* norm, pt1, pt2, pt3 */ );


/*****************************************************************
 * TAG( sqr_dist_seg_to_line )
 *
 * Returns the square of the distance between a line segment and a line.
 * The line segment is specified by its two endpoints; the line is
 * specified by a point and a direction vector.
 */
extern float
sqr_dist_seg_to_line( /* seg_pt1, seg_pt2, line_pt, line_dir */ );


/*****************************************************************
 * TAG( near_pt_on_line )
 *
 * Returns the point on a line nearest to the specified point.
 * The line is specified by a point on the line and a direction
 * vector.  The result is returned in near_pt.
 */
extern void near_pt_on_line( /* pt, line_pt, line_dir, near_pt */ );

/*****************************************************************
 * TAG( intersect_line_plane )
 *
 * Find the intersection of a line and a plane.  The routine
 * returns 0 if they do not intersect.  The intersection point is
 * returned in the last argument.
 */
extern int
intersect_line_plane( /* line_pt, line_dir, pl_pt, pl_norm, isect_pt */ );


/*****************************************************************
 * TAG( intersect_line_tri )
 *
 * Test whether a line intersects a triangular polygon.  Returns
 * 1 if the line intersects and 0 if not.  The intersection point
 * is returned in the last parameter.
 */
extern int intersect_ray_tri( /* verts, line_pt, line_dir, isect_pt */ );


/*****************************************************************
 * TAG( intersect_line_quad )
 *
 * Test whether a line intersects a quadrilateral polygon.  Returns
 * 1 if the line intersects and 0 if not.  The intersection point
 * is returned in the last parameter.
 */
extern int intersect_line_quad( /* verts, line_pt, line_dir, isect_pt */ );


/*****************************************************************
 * TAG( area_of_triangle )
 *
 * Returns the area of a triangular polygon.
 */
extern float area_of_triangle( /* verts */ );


/*****************************************************************
 * TAG( area_of_quad )
 *
 * Returns the area of a quadrilateral polygon.
 */
extern float area_of_quad( /* verts */ );


/*
 * Matrix operations.
 *
 * Define a homogenous transformation matrix. The first row (vector) 
 * is the new X axis, i.e. the X axis in the transformed coordinate 
 * system. The second row is the new Y axis, and so on. The last row
 * is the translation, for a transformed point.
 *
 * The reason we surround the rows with a struct is that we
 * don't want to say (Transf_mat *) &foo[0] instead of &foo when
 * sending an address to a matrix as a parameter to a function.
 * Alas, arrays are not first class objects in C.
 */

/*****************************************************************
 * TAG( Transf_mat )
 *
 */
typedef struct {
    float   mat[4][4];
} Transf_mat;


/*****************************************************************
 * TAG( ident_matrix )
 *
 */
extern Transf_mat ident_matrix;


/*****************************************************************
 * TAG( mat_copy ) 
 *
 * A = B
 */
extern void mat_copy( /* A, B */ );


/*
 * Function declarations for the functions in geometric.c
 */
extern Transf_mat  * transf_mat_create(/* Matrix * */);
extern void          mat_translate(/* Matrix *, float, float, float */);
extern void          mat_rot(/* Matrix *, float angle, int axis */);
extern void          mat_rotate_x(/* Matrix *, float */);
extern void          mat_rotate_y(/* Matrix *, float */);
extern void          mat_rotate_z(/* Matrix *, float */);
extern void          mat_rotate(/* Matrix *, Vector, Vector, float */);
extern void          mat_rotate_axis_vec(/* Matrix *, int, Vector, Bool */);
extern void          mat_rotate_vec_vec(/* Matrix *, Vector, Vector */);
extern void          invert_rot_mat(/* Matrix *, Matrix * */);
extern void          mat_scale(/* Matrix *, float, float, float */);
extern void          mat_mirror_plane(/* Matrix *, Vector, Vector */);
extern void          mat_mul(/* Matrix *, Matrix *, Matrix * */);
extern void          mat_to_angles(/* Matrix *, Vector */);
extern void          point_transform(/* Vector, Vector, Matrix * */);
extern void          mat_to_array( /* Matrix *, float * */);

#endif GEOMETRIC_H

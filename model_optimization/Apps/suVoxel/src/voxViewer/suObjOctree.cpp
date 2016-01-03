#include "suObjOctree.h"
/* Normal Vectors */

float CallbackDraw::cube_n[18] =
{
	0.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f,
	-1.0f, 0.0f, 0.0f,
	0.0f, -1.0f, 0.0f,
	0.0f, 0.0f, -1.0f
};

/* Vertex indices, as quads, before triangulation */
GLubyte CallbackDraw::cube_vi[24] =
{
	7, 4, 5, 6,
	3, 7, 6, 2,
	2, 6, 5, 1,
	5, 4, 0, 1,
	3, 0, 4, 7,
	2, 1, 0, 3
};
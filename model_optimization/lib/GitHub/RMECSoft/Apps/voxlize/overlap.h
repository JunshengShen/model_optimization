#ifndef overlap_h_
#define overlap_h_

#include <cmath>

#define OVERLAP_X   0
#define OVERLAP_Y   1
#define OVERLAP_Z   2

int TriBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3]);

#endif //

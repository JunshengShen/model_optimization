#ifndef PREDICATES_H
#define PREDICATES_H
//#define SINGLE

#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif 


REAL
exactinit(void); // call this before anything else

REAL
orient2d(REAL *pa,
         REAL *pb,
		 REAL *pc);

REAL
orient3d(REAL *pa,
         REAL *pb,
         REAL *pc,
         REAL *pd);

REAL
incircle(REAL *pa,
         REAL *pb,
         REAL *pc,
         REAL *pd);

REAL
insphere(REAL *pa,
         REAL *pb,
         REAL *pc,
         REAL *pd,
         REAL *pe);

#endif

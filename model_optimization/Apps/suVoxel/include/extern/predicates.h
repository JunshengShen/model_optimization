#ifndef PREDICATES_H
#define PREDICATES_H
#ifdef SINGLE
#define REAL REAL
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

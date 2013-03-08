#ifndef MacroDefinitions
#define MacroDefinitions

#define SIGN(x)            ((x >= 0.0) ? 1.0 : -1.0)

#define SIGN2(a,b)         ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MAX(a, b)          ((a) > (b) ? (a) : (b))

#define MIN(a, b)          ((a) > (b) ? (b) : (a))

#define SWAP(x, y, tmp)    (tmp) = (x); (x) = (y); (y) = (tmp)

#define DIST(x, y, z)      (sqrt(((x)*(x) + (y)*(y) + (z)*(z))))

#define SQ_DIST(x, y, z)   ((x)*(x) + (y)*(y) + (z)*(z))

#define NOTEQ(x, y)        (fabs((x) - (y)) > MACH_PREC_ZERO)

#define DNEQ(x, y)         ((fabs(x) > 1.0 && fabs(y) > 1.0 && fabs((x) - (y))/fabs(0.5*((x) + (y))) > 1.0e-12) || fabs((x) - (y)) > 1.0e-12)

#endif

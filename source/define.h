#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Ninterpolate 16
//#define RKFn 13			//RKF78
//#define RKFn 7			//DP54
//#define RKFn 6			//RKF45
#define RKFn 4			//RK4

#define Nperturbers 27



//#define def_dayUnit 0.01720209895
 //#define def_dayUnit 0.01720209894846
//#define def_k2 1.0

 //#define def_c 10065.3201686
//#define def_c 10065.320121


#define def_dayUnit 1.0
#define def_k2 0.01720209895 * 0.01720209895
#define def_c 173.144			//AU / day

#define def_AU 149597870700.0           //AU in m

#define def_atol 1.0e-16
#define def_rtol 1.0e-16
#define def_fac 0.84			//safety factor for RKF integrator
#define def_facmax 1.5 
#define def_facmin 0.8


#define def_OldShuffle 0                //set this to 1 when an old cuda version is used which doesn't have shfl_sync operations

//constant memory
__constant__ double a_c[20 * 20];       //20 is considered here to be large enough (>RKFn)
__constant__ double b_c[20];
__constant__ double bb_c[20];
__constant__ double c_c[20];


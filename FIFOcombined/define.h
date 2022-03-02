#define dayUnit 0.01720209895
//#define dayUnit 0.01720209894846
//#define dayUnit 1.0

//#define def_c 10065.3201686
#define def_c 10065.320121
#define def_AU 149597870700.0           //AU in m

//constant memory
__constant__ double a_c[20 * 20];       //20 is considered here to be large enough (>RKFn)
__constant__ double b_c[20];
__constant__ double c_c[20];


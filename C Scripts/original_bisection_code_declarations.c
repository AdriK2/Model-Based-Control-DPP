#include <math.h>
#include <stdio.h>

#define Im1 InputSignal(0,0)


#define Vout OutputSignal(0,0)
//#define y OutputSignal(0,1)
//#define y1 OutputSignal(0,2)

//#define y2 OutputSignal(0,3)
//#define y3 OutputSignal(0,4)
//#define y4 OutputSignal(0,5)
//#define y5 OutputSignal(0,6)
float T0 = 25.0;
float Vinit = 10.0;
float Isc=9.91;                           //short-circuit current

float Voc=39.71/3;                          //open circuit voltage

float a=1;                             //diode ideality factor
//int Ns=60;                               //series connected cells
int Ns=60/3;                               //series connected cells

float Rs=0.0907;                       //series resistance

float Rp=276.43;                          //parallel resistance
float Rsh=276.43;                          //parallel resistance

//long double k=1.38*10e-23;                      //boltzmann's constant
long double k=1.38e-23;                      //boltzmann's constant
//long double q=1.602*10e-19;                     //electron charge
long double q=1.602e-19;                     //electron charge
int i;
double V1;
double Vout0;
//Iteration criterias
int n=30;                               //number of iterations

double tol=1e-10;

// variables for pv emulation adjustment
static int p, j;
static float Kv;

// K_V = V_emulator / V_ideal at each emulator MPP, ascending current.
//   Row 0 = PV1, row 1 = PV2, row 2 = PV3.  
  
#define NPV 3
#define NKV 6

static const float Ical[NPV][NKV] = {
    {1.4180f, 2.3775f, 3.3265f, 4.2755f, 5.2245f, 6.1930f},    /* PV1 */
    {1.4280f, 2.3775f, 3.3365f, 4.2755f, 5.2345f, 6.1830f},    /* PV2 */
    {1.4185f, 2.3780f, 3.3265f, 4.2850f, 5.2440f, 6.1825f}};   /* PV3 */

static const float Kcal[NPV][NKV] = {
    {1.010296f, 1.012674f, 1.022171f, 1.031162f, 1.039855f, 1.044842f},   /* PV1 */
    {1.002211f, 1.011025f, 1.016896f, 1.028843f, 1.035435f, 1.044105f},   /* PV2 */
    {0.997793f, 1.002549f, 1.012957f, 1.020614f, 1.029509f, 1.040295f}};  /* PV3 */
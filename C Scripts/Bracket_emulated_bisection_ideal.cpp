#include <math.h>

#define Im1  InputSignal(0,0)
#define T0   InputSignal(0,1)
/* InputSignal(0,2) (old Vinit) is no longer used */
#define Vout OutputSignal(0,0)

/* Physical constants (exact SI, same as MATLAB) */
#define Q_E 1.602176634e-19
#define K_B 1.380649e-23

/* PV parameters: same values as the MATLAB model */
static const double Isc = 9.91;          /* short-circuit current [A] */
static const double Voc = 39.71/3.0;     /* open-circuit voltage [V] */
static const double A   = 1.0;           /* diode ideality factor */
static const double Ns  = 20.0;          /* series cells (60/3) */
static const double Rs  = 0.0907;        /* series resistance [ohm] */
static const double Rsh = 276.4340;      /* shunt resistance [ohm] */

static double a1, Is1, VocT;             /* temperature-dependent */
static double V_prev;                    /* warm start: last solution */

static void pv_params(double T0C)
{
    double T1 = T0C + 273.15;
    double B  = -0.28*Voc/100.0;         /* Voc temp. coeff. [V/K] */
    VocT = Voc + B*(T1 - 298.15);
    a1   = Q_E/(Ns*K_B*T1*A);
    Is1  = (Isc + Isc*Rs/Rsh - VocT/Rsh) /
           (exp(a1*VocT) - exp(a1*Isc*Rs));
}

/* MPP condition f(V) = I + V*dI/dV, and its partial derivative wrt V */
static double mpp_f(double V, double I, double *df)
{
    double arg = a1*(V + I*Rs), e, N, D, dN, dD;
    if (arg > 700.0) arg = 700.0;
    e  = exp(arg);
    N  = -Is1*a1*e - 1.0/Rsh;
    D  = 1.0 + Rs/Rsh + a1*Is1*Rs*e;
    dN = -Is1*a1*a1*e;
    dD =  a1*a1*Is1*Rs*e;
    *df = N/D + V*(dN*D - N*dD)/(D*D);
    return I + V*N/D;
}

/* Ideal MPP voltage at current Im: safeguarded Newton on [0, VocT] */
static double vmpp_ideal(double Im)
{
    double I = Im, lo = 0.0, hi = VocT, V = V_prev, f, df, Vn;
    int it;

    if (I < 0.0) I = 0.0;
    if (I > Isc) I = Isc;
    if (!(V > lo && V < hi)) V = 0.8*VocT;

    for (it = 0; it < 20; it++) {
        f = mpp_f(V, I, &df);
        if (f > 0.0) lo = V; else hi = V;
        Vn = V - f/df;
        if (!(Vn > lo && Vn < hi)) Vn = 0.5*(lo + hi);
        if (fabs(Vn - V) < 1e-6) { V = Vn; break; }
        V = Vn;
    }
    V_prev = V;
    return V;
}


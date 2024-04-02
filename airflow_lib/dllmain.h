#pragma once

#include <math.h>

#define pi 3.14159265358979323846
#define e 2.71828182845904523536
#define Ts 293.15      // стандартная температура в кельвинах (ГОСТ 30319.0 п.4.3)
#define Ps 101325      // стандартное давление в Паскалях (ГОСТ 30319.0 п.4.3)
#define Pls 1.20445    // плотность воздуха при н.у., кг / м3 (ГОСТ 30319.1-96 таблица 1)



double __cdecl get_C(double, double, double);

double __cdecl get_Ksh(double, double, double, double);

double __cdecl get_qm(double, double, double, double, double, double, double, double, double, double, double, double);

double __cdecl get_true_qm(double, double, double, double, double, double, double, double, double, double, double, double);

double __cdecl get_A(int, double);

double __cdecl get_y(double, double, double, double, double, bool mark = false);

double __cdecl get_Re(double, double, double);

extern "C" __declspec(dllexport) double __cdecl calc_airflow_mu(double, double, double, double, double);
// dllmain.cpp : Определяет точку входа для приложения DLL.
#include "pch.h"
#include "dllmain.h"



BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
        break;
    case DLL_THREAD_ATTACH:
        break;
    case DLL_THREAD_DETACH:
        break;
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}



double get_C(double Re, double B, double D) {
    double A = pow((19000. * B / Re), 0.8);
    double L1 = 0, L2 = 0; // для углового способа отбора давлений
    double M1 = 2 * L2 / (1 - B);
    double M2 = 0;
    if (D < 0.07112) {
        M2 = 0.011 * (0.75 - B) * (2.8 - (D / 0.0254));
    }
    return 0.5961 + 0.0261 * pow(B, 2) - 0.216 * pow(B, 8) + 0.000521 * pow(((B * 1e6) / Re), 0.7) + (
        0.0188 + 0.0063 * A) * pow(B, 3.5) * pow((1e6 / Re), 0.3) + (
            0.043 + 0.08 * pow(e, (-10 * L1)) - 0.123 * pow(e, (-7 * L1))) * (1 - 0.11 * A) * (
                pow(B, 4) / (1 - pow(B, 4))) - 0.031 * (M1 - 0.8 * pow(M1, 1.1)) * pow(B, 1.3) + M2;
}

double get_A(int i, double Re) {
    double tab[3][3][4] = {
        { // index 1
            {8.87, -3.7114, 0.41841, 0},
            {6.7307, -5.5844, 0.732485, 0},
            {-10.244, 5.7094, 0.76477, 0}
        },
        { // index 2
            {27.23, -11.458, 1.6117, -0.07567},
            {-25.928, 12.426, -2.09397, 0.106143},
            {1.7622, -3.8765, 1.05567, -0.076764}
        },
        { // index 3
            {16.5416, -6.60709, 0.88147, -0.039226},
            {322.594, -132.2, 17.795, -0.799765},
            {-92.029, 37.935, -5.1885, 0.23583}
        }
    };

    double sum = 0;
    int tab_index = (Re <= 1e5) ? 0 : (Re <= 3e6) ? 1 : 2;

    for (int k = 0; k < 4; ++k) {
        sum += tab[tab_index][i][k] * pow(fabs(log10(Re)), k);
    }

    return sum;
}

double get_y(double Re, double Ra, double Ra_max, double Ra_min, double D, bool  mark) {
    double Rsh = 0.0, kD = 0.0, kR = 5.035 / Re;

    if (!mark) {
        Rsh = pi * Ra;
    }
    else if (Ra > Ra_max) {
        Rsh = pi * Ra_max;
    }
    else {
        Rsh = pi * Ra_min;
    }

    kD = 0.26954 * Rsh / D;

    return pow((1.74 - 2 * log10((2 * Rsh / D) - ((37.36 * log10(kD - kR * log10(kD + 3.3333 * kR))) / Re))), -2);
}

double get_Ksh(double Re, double Ra, double B, double D) {
    double leftRa_max = 0, leftRa_min = 0, Ra_max = 0, Ra_min = 0;
    // Ra_max   
    if (Re <= 1e4) {
        leftRa_max = 0.718866 * pow(B, (-3.887)) + 0.364;
    }
    else if (B < 0.65) {
        leftRa_max = get_A(0, Re) * pow(B, get_A(1, Re)) + get_A(2, Re);
    }
    else {
        leftRa_max = get_A(0, Re) * pow(0.65, get_A(1, Re)) + get_A(2, Re);
    }

    if (leftRa_max >= 15) {
        Ra_max = 15e-4 * D;
    }
    else {
        Ra_max = leftRa_max * 1e-4 * D;
    }

    // Ra_min
    if (B < 0.65) {
        leftRa_min = 7.1592 - 12.387 * B - (20.118 - 3.469 * B) * log10(Re) + (0.1382 - 0.23762 * B) * pow(log10(Re), 2);
    }
    else {
        leftRa_min = -0.892352 + 0.24308 * log10(Re) - 0.0162562 * pow(log10(Re), 2);
    }

    if (leftRa_min <= 0 or Re < 3e6) {
        Ra_min = 0;
    }
    else {
        Ra_min = 1e-4 * D * leftRa_min;
    }

    if (Ra >= Ra_min and Ra <= Ra_max) {
        return 1;
    }
    else {
        double y0 = get_y(Re, Ra, Ra_max, Ra_min, D), y1 = get_y(Re, Ra, Ra_max, Ra_min, D, true);

        return 1 + 5.22 * pow(B, 3.5) * (y0 - y1);
    }


}

double get_qm(double Re_cur, double d, double Ksu, double E, double B, double D, double Ra, double Kp, double eps, double dP, double T, double P) {
    double C = get_C(Re_cur, B, D); // (ГОСТ 8.586.2 п.5.3.2.1)
    double Ksh = get_Ksh(Re_cur, Ra, B, D); // (ГОСТ 8.586.2 п.5.3.2.2)
    return 0.25 * pi * pow(d, 2) * C * E * Ksh * Kp * eps * pow((2 * dP * Pls * ((P * Ts) / (Ps * T))), 0.5);
}

double get_true_qm(double u, double d_20, double Ksu, double E, double B, double D, double Ra, double Kp, double eps, double dP, double T, double P) {
    double Re_first = 1e6;
    double qm_first = get_qm(Re_first, d_20, Ksu, E, B, D, Ra, Kp, eps, dP, T, P); // (ГОСТ 8.568.5 п.5.2.3)
    double Re_cur = get_Re(qm_first, D, u); // (ГОСТ 8.568.5 п.5.2.5)
    double qm_cur = get_qm(Re_cur = Re_cur, d_20, Ksu, E, B, D, Ra, Kp, eps, dP, T, P);
    while (100 * fabs(qm_cur - qm_first) / qm_first >= 0.001) {
        qm_first = qm_cur;
        Re_cur = get_Re(qm_first, D, u);
        qm_cur = get_qm(Re_cur = Re_cur, d_20, Ksu, E, B, D, Ra, Kp, eps, dP, T, P);
    }

    return round(3600 * qm_cur); // перевод в кг/ч 
}

double get_Re(double qm, double D, double u) {
    return (4 * qm) / (pi * D * u);
};

__declspec(dllexport) double __cdecl calc_airflow_mu(double p_izm, double dp_izm, double t_izm, double d_20, double D_20) {
    // Исходные данные

    double Ra = 0.045e-3; // Среднее арифметическое отклонение профиля шероховатости ИТ(стальные трубы с незначительным налетом ржавчины, ГОСТ 8.586.2 Приложение Д), м
    double r_s = 0.00004;  // Начальный радиус входной кромки диафрагмы, м
    double ty = 1;  // межконтрольный интервал СУ, лет
    double F = 0;  // относительная влажность газа %
    double dP = dp_izm * 98070;  // перепад на диафрагме, Па
    double P = p_izm * 98070;  // абсолютное давление на диафрагме, кгс / см2
    double t = t_izm;  // температура газа, С

    // Расчет

    // Коэффициент, учитывающий изменение диаметра отверстия диафрагмы вызванное отклонением температуры газа от 20
    // Постоянные для материала СУ: сталь марки 12Х18Н9Т (а0=15.6; а1=8.3; а2=-6.5)
    double a_tsu = 1e-6 * (15.6 + 8.3 * (t / 1000) - 6.5 * pow((t / 1000), 2)); // (ГОСТ 8.586.1 Приложение Г)
    double Ksu = 1 + a_tsu * (t - 20); // (ГОСТ 8.586.1 п.5.5)

    // Диаметр отверстия диафрагмы при рабочей температуре
    double d = d_20 * Ksu; // (ГОСТ 8.586.1 п.5.5)

    // Коэффициент, учитывающий изменение диаметра ИТ, вызванное отклонением температуры газа от 20
    // Постоянные для материала ИТ: сталь марки 20 (а0=11.1; а1=7.7; а2=-3.4)
    double a_tt = 1e-6 * (11.1 + 7.7 * (t / 1000) - 3.4 * pow((t / 1000), 2)); // (ГОСТ 8.586.1 Приложение Г)
    double Kt = 1 + a_tt * (t - 20); // (ГОСТ 8.586.1 п.5.5)

    // Диаметр ИТ при рабочей температуре
    double D = D_20 * Kt; // (ГОСТ 8.586.1 п.5.5)

    // Относительный диаметр отверстия диафрагмы
    double B = d / D; // (ГОСТ 8.586.1 п.3.2.12)

    // Коэффициент скорости входа
    double E = 1 / pow((1 - pow(B, 4)), 0.5); // (ГОСТ 8.586.1 п.3.3.9)

    // Поправочный коэффициент учитывающий притупление входной кромки диафрагмы (ГОСТ 8.586.2 п. 5.3.2.4)
    double ark = 0.195e-3; // (ГОСТ 8.586.2 п. 5.3.2.4 описание формулы 5.14)
    double rk = ark - (3. / ty) * (ark - r_s) * (1 - pow(e, -ty / 3)); // (ГОСТ 8.586.2 п. 5.3.2.4 формула 5.15)
    double Kp = 0.9826 + pow((rk / d + 0.0007773), 0.6); // (ГОСТ 8.586.2 п. 5.3.2.4 формула 5.16)




    double T = t + 273.15;  // Термодинамическая температура газа, К (ГОСТ 8.586.5 п.6.3.1)

    double u = 18.27 * 1e-6 * ((291.15 + 120) / (T + 120)) * pow((T / 291.15), (3 / 2.));

    // показатель адиабаты воздуха при измеренной температуре (апроксимация табличных значений полиномом 4й степени)
    double k = -1E-13 * pow(t, 4) + 4E-10 * pow(t, 3) - 3E-07 * pow(t, 2) - 9E-06 * t + 1.4011;

    // Коэффициент расширения (ГОСТ 8.586.2 п.5.3.2.2)
    double eps = 1 - (0.351 + 0.256 * pow(B, 4) + 0.93 * pow(B, 8)) * (1 - pow((1 - dP / P), (1 / k)));

    //######################################################
    double Q = get_true_qm(u, d, Ksu, E, B, D, Ra, Kp, eps, dP, T, P);

    // отсечение ошибок измерения датчиком
    if (Q <= 0) {
        return 0;
    }
    else {
        return Q;
    }
}


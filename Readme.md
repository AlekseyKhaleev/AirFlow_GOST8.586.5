﻿# Расчет расхода воздуха по ГОСТ 8.586.5

## Описание

Данная библиотека на C++ предназначена для вычисления расхода воздуха через диафрагму с угловым способом отбора давления в соответствии с ГОСТ 8.568.5. Расчеты проводятся для воздуха с использованием стандартных сужающих устройств (диафрагм), сделанных из стали марки 12Х18Н9Т, установленных в измерительных трубах из стали марки 20.

## Входные параметры функции

Функция расчета принимает следующие параметры:

- `p_izm` - Давление на диафрагме (в кгс/см²).
- `dp_izm` - Перепад давления на диафрагме (в кгс/см²).
- `t_izm` - Температура перед диафрагмой (в градусах Цельсия).
- `d_20` - Диаметр сужающего устройства (внутренний диаметр диафрагмы) (в метрах).
- `D_20` - Диаметр измерительной трубы (в метрах).

## Исходные данные

В функции используются следующие постоянные значения:

- `Ra = 0.045e-3` - Среднее арифметическое отклонение профиля шероховатости измерительной трубы (ГОСТ 8.586.2 Приложение Д), м.
- `r_s = 0.00004` - Начальный радиус входной кромки диафрагмы, м.
- `ty = 1` - Межконтрольный интервал средства учета (СУ), лет.
- `F = 0` - Относительная влажность газа (%).
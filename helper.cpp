#include "helper.h"

void DataProcess(double &mean, double &stdev, double *Time, int cyctimes)
{
    double temp;
    double sum = 0;
    for (int i = 0; i < cyctimes; i++)
    {
        sum = sum + Time[i];
    }
    mean = sum / cyctimes;
    double temp_sum = 0;
    for (int i = 0; i < cyctimes; i++)
    {
        temp = mean - Time[i];
        temp = temp * temp;
        temp_sum = temp_sum + temp;
    }
    stdev = sqrt(temp_sum / cyctimes);
    stdev = stdev / mean;
}

ZZ PRF_ZZ(int prfkey, ZZ mmod)
{
    ZZ res;
    SetSeed(ZZ(prfkey));
    RandomBnd(res, mmod);
    return res;
}


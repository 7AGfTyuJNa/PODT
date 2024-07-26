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

ZZ PRF_ZZ(const int &prfkey, const ZZ &mmod)
{
    ZZ res;
    SetSeed(ZZ(prfkey));
    RandomBnd(res, mmod);
    return res;
}


void GenerateMatrix(int m, int n, std::vector<std::vector<int>>& delta)
{
    srand(time(0)); // Initialize random seed
    
    for (int i = 0; i < m; ++i)
    {
        // Initialize the row with 0s
        std::vector<int> row(n, 0);
        
        // Randomly select an index to set to 1
        int randomIndex = rand() % n;
        row[randomIndex] = 1;
        
        // Assign the row to the matrix
        delta[i] = row;
    }
}
#pragma once
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

using namespace std;
using namespace NTL;

void DataProcess(double &mean, double &stdev, double *Time, int cyctimes);
ZZ PRF_ZZ(const int &prfkey, const ZZ &mmod);
void GenerateMatrix(int m, int n, std::vector<std::vector<int>>& delta);
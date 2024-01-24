#pragma once
#include "HSS.h"

typedef struct
{
    PK pk;

    int t; // bit length
    int m; // decision nodes
    int k; // leaf nodes
    int h; // tree depth
} Para;

void KeyGen(Para &param, EK &ek0, EK &ek1);
void KeyFree(Para &param, EK &ek0, EK &ek1);

void ClientEnc(vector<Mat<ZZ>> &Ix, Para param, vec_ZZ x);
void ProviderEnc(vector<Mat<ZZ>> &Ix, Mat<ZZ> &Iv, Para param, vec_ZZ y, vec_ZZ v);
void HSSCMP(ZZ &c_b, int b, Para param, EK ekb, Mat<ZZ> Ix, Mat<ZZ> Iy, int &prf_key, Vec<ZZ> M1b);
void ClassificationGen(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, Para param, EK ekb, Vec<ZZ> cmp_res, Mat<ZZ> Iv, Vec<ZZ> M1, 
int &prf_key);
void DTEvaluation(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, Para param, EK ekb, 
vector<Mat<ZZ>> Ix, vector<Mat<ZZ>> Iy, Mat<ZZ> Iv);
void ClDecryption(ZZ &res, Para param, Vec<ZZ> &pc_0, Vec<ZZ> &vv_0, Vec<ZZ> &pc_1, Vec<ZZ> &vv_1);

void ODTE_TIME_TEST(int depth, int msgbit, int cyctimes);
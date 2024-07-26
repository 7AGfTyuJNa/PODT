#pragma once
#include "HSS.h"

typedef struct
{
    PK pk;

    int t; // bit length
    int m; // decision nodes
    int k; // leaf nodes
    int h; // tree depth
    int n; // attributes
} Para;

void KeyGen(Para &param, EK &ek0, EK &ek1);
void KeyFree(Para &param, EK &ek0, EK &ek1);

void FeatureSelection1(Mat<vec_ZZ> &Ix, const Para &param, const vec_ZZ &x);
void FeatureSelection2(Mat<vec_ZZ> &Ix, const Para &param, const vec_ZZ &x, const Mat<vec_ZZ> &Idelta);
void ProviderEnc(Mat<vec_ZZ> &Idelta, Mat<vec_ZZ> &Iy, Mat<ZZ> &Iv, const Para &param, const vec_ZZ &y, const vec_ZZ &v, const vector<vector<int>> &delta);
void HSSCMP(ZZ &c_b, int b, const Para &param, const EK &ekb, const Vec<vec_ZZ> &Ix, const Vec<vec_ZZ> &Iy, int &prf_key, const Vec<ZZ> &M1b);
void ClassificationGen(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, const Para &param, const EK &ekb, const Vec<ZZ> &cmp_res, 
const Mat<ZZ> &Iv, const Vec<ZZ> &M1, int &prf_key);
void DTEvaluation(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, const Para &param, const EK &ekb,
                  const Mat<vec_ZZ> &Ix, const Mat<vec_ZZ> &Iy, const Mat<ZZ> &Iv);
void ClDecryption(ZZ &res, const Para &param, const Vec<ZZ> &pc_0, const Vec<ZZ> &vv_0, const Vec<ZZ> &pc_1, const Vec<ZZ> &vv_1);
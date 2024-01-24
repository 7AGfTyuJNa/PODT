#pragma once
#include "helper.h"

typedef struct
{
    ZZ g1;
    ZZ g2;
    ZZ g3;
} BeaT;

typedef struct
{
    Mat<ZZ> G1;
    vec_ZZ G2;
    vec_ZZ G3;
} BeaMaT;

typedef struct
{
    ZZ P; // Large Prime
    int n; // input length
    int t; // bit length
    int m; // decision nodes
    int k; // leaf nodes
    int h; // tree depth
} TDSCPara;

void ParaGen(TDSCPara &pk, BeaT &BT0, BeaT &BT1, BeaMaT &BMT0, BeaMaT &MBT1);
// Input Preparation
void ClientEnc(vec_ZZ &x0, vec_ZZ &x1, TDSCPara pk, vec_ZZ x);
void ProviderEnc(vec_ZZ &y0, vec_ZZ &y1, vec_ZZ &v0, vec_ZZ &v1,
                 Mat<ZZ> &mat0, Mat<ZZ> &mat1, TDSCPara pk, vec_ZZ y, vec_ZZ v, Mat<ZZ> matt);

void InputSelection(vec_ZZ &res0, vec_ZZ &res1, Mat<ZZ> M0, Mat<ZZ> M1, TDSCPara pk, 
BeaMaT BMT0,
    BeaMaT BMT1,
vec_ZZ x0, vec_ZZ x1);
void SSCMP(ZZ &cmp_res0, ZZ &cmp_res1, TDSCPara pk, BeaT BT0,
    BeaT BT1, ZZ x0, ZZ y0, ZZ x1, ZZ y1);
void ClassificationGen(vec_ZZ &pc0, vec_ZZ &pc1, vec_ZZ &vv0, vec_ZZ &vv1, TDSCPara pk, 
vec_ZZ cmp_res0, vec_ZZ cmp_res1, vec_ZZ v0, vec_ZZ v1);
void DTEvaluation(vec_ZZ &pc0, vec_ZZ &pc1, vec_ZZ &vv0, vec_ZZ &vv1, TDSCPara pk, 
BeaT BT0, BeaT BT1, BeaMaT BMT0, BeaMaT MBT1,
vec_ZZ x0, vec_ZZ x1, vec_ZZ y0, vec_ZZ y1, vec_ZZ v0, vec_ZZ v1,
                 Mat<ZZ> mat0, Mat<ZZ> mat1);
void Decrypt(ZZ &res, vec_ZZ pc0, vec_ZZ pc1, vec_ZZ vv0, vec_ZZ vv1, TDSCPara pk);


void ShareGen(ZZ &res0, ZZ &res1, ZZ modP, ZZ x);
void SSZZMult(ZZ &res0, ZZ &res1, TDSCPara pk, BeaT BT0,
    BeaT BT1, ZZ a0, ZZ b0, ZZ a1, ZZ b1);
void SSMatMult(vec_ZZ &res0, vec_ZZ &res1, TDSCPara pk, Mat<ZZ> M0, vec_ZZ x0, Mat<ZZ> M1, vec_ZZ x1);

void MatAdd(Mat<ZZ> &res, TDSCPara pk, Mat<ZZ> x, Mat<ZZ> y);
void MatSub(Mat<ZZ> &res, TDSCPara pk, Mat<ZZ> x, Mat<ZZ> y);
void MatMul(vec_ZZ &res, TDSCPara pk, Mat<ZZ> x, vec_ZZ y);

void TDSC20_TIME_TEST(int depth, int N_attribute, int msgbit, int cyctimes);

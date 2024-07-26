#pragma once
#include "helper.h"

typedef struct
{
    int l;   // length of ciphertext logN^2 / log B_sk
    ZZ N2;   // N^2
    int k;   // kappa
    ZZ Bmsg; // 2^k
    ZZ Bsk;  // B_msg * B_sk * 2^k = N

    ZZ N;      // B_msg * B_sk * 2^k = N
    bool DJN_OPEN;
    bool Pre_OPEN;
    ZZ hs;     // DJN faster
    Vec<ZZ> D; // Paillier.Enc(d_Bsk)
} PK;

typedef struct
{
    Vec<ZZ> d_Bsk; // d base B_sk
} EK;

void Paillier_Gen(PK &pk, ZZ &d);
void Pailler_Enc(ZZ &ct, const PK &pk, const ZZ &x);
void Pailler_Dec(ZZ &x, const PK &pk, const ZZ &sk, const ZZ &ct);

void HSS_Gen(PK &pk, EK &ek0, EK &ek1);
void HSS_Free(PK &pk, EK &ek0, EK &ek1);

void HSS_Input(Vec<ZZ> &I, const PK &pk, const ZZ &x);
void HSS_M1Gen(Vec<ZZ> &M1, int b, const PK &pk, const EK &ek, int &prf_key);
void HSS_ConvertInput(Vec<ZZ> &Mx, int sigma, const PK &pk, const EK &ek, const Vec<ZZ> &Ix, int &prf_key);
void HSS_Mul(Vec<ZZ> &Mz, const PK &pk, const EK &ek, const Vec<ZZ> &Ix, const Vec<ZZ> &My, int &prf_key);
void HSS_DDLog(ZZ &z, const PK &pk, const ZZ &g);
void HSS_AddMemory(Vec<ZZ> &Mz, const PK &pk, const Vec<ZZ> &Mx, const Vec<ZZ> &My);
void HSS_SubMemory(Vec<ZZ> &Mz, const PK &pk, const Vec<ZZ> &Mx, const Vec<ZZ> &My);

void HSS_AddInput(Vec<ZZ> &Iz, const PK &pk, const Vec<ZZ> &Ix, const Vec<ZZ> &Iy);
void HSS_Evaluate(Vec<ZZ> &y_b_res, int b, const Mat<ZZ> &Ix, const PK &pk, const EK &ekb, int &prf_key, const vector<vector<int>> &F_TEST);
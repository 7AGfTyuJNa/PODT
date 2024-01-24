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
void Pailler_Enc(ZZ &ct, PK pk, ZZ x);
void Pailler_Dec(ZZ &x, PK pk, ZZ sk, ZZ ct);

void HSS_Gen(PK &pk, EK &ek0, EK &ek1);
void HSS_Free(PK &pk, EK &ek0, EK &ek1);

void HSS_Input(Vec<ZZ> &I, PK pk, ZZ x);
void HSS_M1Gen(Vec<ZZ> &M1, int b, PK pk, EK ek, int &prf_key);
void HSS_ConvertInput(Vec<ZZ> &Mx, int sigma, PK pk, EK ek, Vec<ZZ> Ix, int &prf_key);
void HSS_Mul(Vec<ZZ> &Mz, PK pk, EK ek, Vec<ZZ> Ix, Vec<ZZ> My, int &prf_key);
void HSS_DDLog(ZZ &z, PK pk, ZZ g);
void HSS_AddMemory(Vec<ZZ> &Mz, PK pk, Vec<ZZ> Mx, Vec<ZZ> My);
void HSS_SubMemory(Vec<ZZ> &Mz, PK pk, Vec<ZZ> Mx, Vec<ZZ> My);

void HSS_AddInput(Vec<ZZ> &Iz, PK pk, Vec<ZZ> Ix, Vec<ZZ> Iy);
void HSS_Evaluate(Vec<ZZ> &y_b_res, int b, Mat<ZZ> Ix, PK pk, EK ekb, int &prf_key, vector<vector<int>> F_TEST);
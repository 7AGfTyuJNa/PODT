#include "ODTE.h"

void KeyGen(Para &param, EK &ek0, EK &ek1)
{
    HSS_Gen(param.pk, ek0, ek1);
}

void KeyFree(Para &param, EK &ek0, EK &ek1)
{
    HSS_Free(param.pk, ek0, ek1);
}

void ClientEnc(vector<Mat<ZZ>> &Ix, Para param, vec_ZZ x)
{
    Mat<ZZ> tmp;
    tmp.SetDims(param.t, param.pk.l + 1);
    for (int i = 0; i < param.m; ++i)
    {
        for (int j = 0; j < param.t; ++j)
        {
            HSS_Input(tmp[j], param.pk, ZZ(bit(x[i], j)));
        }
        Ix.push_back(tmp);
    }
}

void ProviderEnc(vector<Mat<ZZ>> &Ix, Mat<ZZ> &Iv, Para param, vec_ZZ y, vec_ZZ v)
{
    Mat<ZZ> tmp;
    tmp.SetDims(param.t, param.pk.l + 1);
    for (int i = 0; i < param.m; ++i)
    {
        for (int j = 0; j < param.t; ++j)
        {
            HSS_Input(tmp[j], param.pk, ZZ(bit(y[i], j)));
        }
        Ix.push_back(tmp);
    }

    Iv.SetDims(param.k, param.pk.l + 1);
    for (int j = 0; j < param.k; ++j)
    {
        HSS_Input(Iv[j], param.pk, v[j]);
    }
}

void HSSCMP(ZZ &c_b, int b, Para param, EK ekb, Mat<ZZ> Ix, Mat<ZZ> Iy, int &prf_key, Vec<ZZ> M1b)
{
    Vec<ZZ> Mxb, Mxyb, Mcb;
    HSS_Mul(Mxb, param.pk, ekb, Ix[0], M1b, prf_key);
    HSS_Mul(Mxyb, param.pk, ekb, Iy[0], Mxb, prf_key);
    HSS_SubMemory(Mcb, param.pk, Mxb, Mxyb);

    for (int i = 1; i < param.t; ++i)
    {
        Vec<ZZ> Mcxb, Mcyb, Mcxyb;

        HSS_Mul(Mcxb, param.pk, ekb, Ix[i], Mcb, prf_key);
        HSS_Mul(Mcyb, param.pk, ekb, Iy[i], Mcb, prf_key);
        HSS_Mul(Mcxyb, param.pk, ekb, Ix[i], Mcyb, prf_key);
        HSS_Mul(Mxb, param.pk, ekb, Ix[i], M1b, prf_key);
        HSS_Mul(Mxyb, param.pk, ekb, Iy[i], Mxb, prf_key);

        HSS_SubMemory(Mcb, param.pk, Mcb, Mcxb);
        HSS_SubMemory(Mcb, param.pk, Mcb, Mcyb);
        HSS_AddMemory(Mcb, param.pk, Mcb, Mcxyb);
        HSS_AddMemory(Mcb, param.pk, Mcb, Mcxyb);
        HSS_AddMemory(Mcb, param.pk, Mcb, Mxb);
        HSS_SubMemory(Mcb, param.pk, Mcb, Mxyb);
    }
    c_b = Mcb[0];
}

void ClassificationGen(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, Para param, EK ekb, Vec<ZZ> cmp_res, Mat<ZZ> Iv, Vec<ZZ> M1, int &prf_key)
{
    vec_ZZ vb;
    ZZ pc, v, r_0, r_1;
    int idx;

    for (int i = 0; i < param.k; ++i)
    {
        double tt = GetTime();
        HSS_Mul(vb, param.pk, ekb, Iv[i], M1, prf_key);

        // compute pc_k
        pc = 0;
        idx = param.m + i;

        while (idx)
        {
            bool lr = (idx - 1) & 1; // & 1
            idx = (idx - 1) / 2;
            if (lr)
            {
                AddMod(pc, pc, cmp_res[idx], param.pk.N); // right
            }
            else
            {
                pc = (pc + b - cmp_res[idx]) % param.pk.N; // left
            }
        }
        // random
        r_0 = PRF_ZZ(prf_key++, param.pk.N);
        r_1 = PRF_ZZ(prf_key++, param.pk.N);

        MulMod(r_1, r_1, pc, param.pk.N);
        AddMod(v, r_1, vb[0], param.pk.N);

        MulMod(pc, r_0, pc, param.pk.N);

        pc_b.append(pc);
        vv_b.append(v);
        tt = GetTime() - tt;
        cout << "Each leaf node Eval time: " << tt * 1000 << " ms\n";
    }
}

void DTEvaluation(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, Para param, EK ekb,
                  vector<Mat<ZZ>> Ix, vector<Mat<ZZ>> Iy, Mat<ZZ> Iv)
{
    int prf_key = 0;

    Vec<ZZ> M1, cmp_res;
    HSS_M1Gen(M1, b, param.pk, ekb, prf_key);

    cmp_res.SetLength(param.m);
    for (int i = 0; i < param.m; ++i)
    {
        HSSCMP(cmp_res[i], b, param, ekb, Ix[i], Iy[i], prf_key, M1);
    }

    ClassificationGen(pc_b, vv_b, b, param, ekb, cmp_res, Iv, M1, prf_key);
}

void ClDecryption(ZZ &res, Para param, Vec<ZZ> &pc_0, Vec<ZZ> &vv_0, Vec<ZZ> &pc_1, Vec<ZZ> &vv_1)
{
    for (int i = 0; i < param.k; ++i)
    {
        if (pc_1[i] - pc_0[i] == 0)
        {
            res = vv_1[i] - vv_0[i];
            break;
        }
    }
}
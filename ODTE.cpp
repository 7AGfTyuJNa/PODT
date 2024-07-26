#include "ODTE.h"

void KeyGen(Para &param, EK &ek0, EK &ek1)
{
    HSS_Gen(param.pk, ek0, ek1);
}

void KeyFree(Para &param, EK &ek0, EK &ek1)
{
    HSS_Free(param.pk, ek0, ek1);
}

void FeatureSelection2(Mat<vec_ZZ> &Ix, const Para &param, const vec_ZZ &x, const Mat<vec_ZZ> &Idelta)
{
    int i, j, k;
    Ix.kill();
    Ix.SetDims(param.m, param.t);
    vec_ZZ C1, TEMP;
    HSS_Input(C1, param.pk, ZZ(1));

    for (i = 0; i < param.m; ++i)
    {
        for (j = 0; j < param.t; ++j)
        {
            TEMP = C1;
            for (k = 0; k < param.n; ++k)
            {
                if (bit(x[k], j) == 1)
                {
                    HSS_AddInput(TEMP, param.pk, TEMP, Idelta[i][k]);
                }
            }
            Ix[i][j] = TEMP;
        }
    }
}

void ProviderEnc(Mat<vec_ZZ> &Idelta, Mat<vec_ZZ> &Iy, Mat<ZZ> &Iv, const Para &param, const vec_ZZ &y, const vec_ZZ &v, const vector<vector<int>> &delta)
{
    int i, j;

    Iy.SetDims(param.m, param.t);
    for (i = 0; i < param.m; ++i)
    {
        for (j = 0; j < param.t; ++j)
        {
            HSS_Input(Iy[i][j], param.pk, ZZ(bit(y[i], j)));
        }
    }

    Iv.SetDims(param.k, param.pk.l + 1);
    for (j = 0; j < param.k; ++j)
    {
        HSS_Input(Iv[j], param.pk, v[j]);
    }

    Idelta.SetDims(param.m, param.n);
    for (int i = 0; i < param.m; ++i)
    {
        for (int j = 0; j < param.n; ++j)
        {
            HSS_Input(Idelta[i][j], param.pk, ZZ(delta[i][j]));
        }
    }
}

void HSSCMP(ZZ &c_b, int b, const Para &param, const EK &ekb, const Vec<vec_ZZ> &Ix, const Vec<vec_ZZ> &Iy, int &prf_key, const Vec<ZZ> &M1b)
{
    Vec<ZZ> Mxb, Mxyb, Mcb,Mcxb, Mcyb, Mcxyb;
    HSS_Mul(Mxb, param.pk, ekb, Ix[0], M1b, prf_key);
    HSS_Mul(Mxyb, param.pk, ekb, Iy[0], Mxb, prf_key);
    HSS_SubMemory(Mcb, param.pk, Mxb, Mxyb);

    for (int i = 1; i < param.t; ++i)
    {
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

void ClassificationGen(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, const Para &param, const EK &ekb, const Vec<ZZ> &cmp_res,
                       const Mat<ZZ> &Iv, const Vec<ZZ> &M1, int &prf_key)
{
    vec_ZZ vb;
    ZZ pc, v, r_0, r_1;
    int idx;

    for (int i = 0; i < param.k; ++i)
    {
        double tt;
        tt = GetTime();
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
        cout << tt *1000<< endl;
    }
}

void DTEvaluation(Vec<ZZ> &pc_b, Vec<ZZ> &vv_b, int b, const Para &param, const EK &ekb,
                  const Mat<vec_ZZ> &Ix, const Mat<vec_ZZ> &Iy, const Mat<ZZ> &Iv)
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

void ClDecryption(ZZ &res, const Para &param, const Vec<ZZ> &pc_0, const Vec<ZZ> &vv_0, const Vec<ZZ> &pc_1, const Vec<ZZ> &vv_1)
{
    for (int i = 0; i < param.k; ++i)
    {
        if (pc_1[i] - pc_0[i] == 0)
        {
            res = vv_1[i] - vv_0[i];
            // break;
        }
    }
}
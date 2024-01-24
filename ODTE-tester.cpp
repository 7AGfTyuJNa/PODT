#include "ODTE-tester.h"

void ODTE_TIME_TEST(int depth, int msgbit, int cyctimes, bool de bu g)
{
    double t1, t2;

    // Start: Setup Test
    t1 = GetTime();

    Para param;
    EK ek0, ek1;
    param.h = depth;
    param.k = 2 << (depth - 1);
    param.m = param.k - 1;
    param.t = msgbit;
    KeyGen(param, ek0, ek1);

    t2 = GetTime() - t1;
    cout << "Setup algo time: " << t2 * 1000 << " ms\n";
    // End: Setup Test

    // Start: Random Input
    Vec<ZZ> X, Y, V;
    X.SetLength(param.m); // Client
    Y.SetLength(param.m); // Provider
    V.SetLength(param.k); // Classification
    for (int i = 0; i < param.m; ++i)
    {
        RandomBits(X[i], param.t);
        RandomBits(Y[i], param.t);
        RandomBits(V[i], param.t);
    }
    RandomBits(V[param.m], param.t);

    int true_idx = 0;
    while (true_idx < param.m)
    {
        if (X[true_idx] > Y[true_idx])
        {
            true_idx = 2 * true_idx + 1;
        }
        else
        {
            true_idx = 2 * true_idx + 2;
        }
    }
    true_idx -= param.m; // V[true_idx] is final result.

    // End: Random Input

    // Start: Client Input Test
    t1 = GetTime();

    vector<Mat<ZZ>> Ix;
    if (debug)
    {
        for (int i = 0; i < X.length(); ++i)
        {
            Mat<ZZ> tmp;
            tmp.SetDims(param.t, param.pk.l + 1);
            for (int j = 0; j < param.t; ++j)
            {
                for (int k = 0; k < param.pk.l + 1; ++k)
                {
                    RandomBnd(tmp[j][k], param.pk.N2);
                }
            }
            Ix.push_back(tmp);
        }
    }
    else
    {
        Mat<ZZ> tmp;
        tmp.SetDims(param.t, param.pk.l + 1);
        for (int i = 0; i < param.m; ++i)
        {
            for (int j = 0; j < param.t; ++j)
            {
                HSS_Input(tmp[j], param.pk, ZZ(bit(X[i], j)));
            }
            Ix.push_back(tmp);
        }
    }

    t2 = GetTime() - t1;
    cout << "Client Input algo time: " << t2 * 1000 << " ms\n";
    // End: Client Input Test

    // Start: Provider Input Test
    t1 = GetTime();
    vector<Mat<ZZ>> Iy;
    Mat<ZZ> Iv;
    // ProviderEnc(Iy, Iv, param, Y, V);
    if (debug)
    {
        for (int i = 0; i < Y.length(); ++i)
        {
            Mat<ZZ> tmp;
            tmp.SetDims(param.t, param.pk.l + 1);
            for (int j = 0; j < param.t; ++j)
            {
                for (int k = 0; k < param.pk.l + 1; ++k)
                {
                    RandomBnd(tmp[j][k], param.pk.N);
                }
            }
            Iy.push_back(tmp);
        }

        Iv.SetDims(param.k, param.pk.l + 1);
        for (int j = 0; j < param.k; ++j)
        {
            for (int k = 0; k < param.pk.l + 1; ++k)
            {
                RandomBnd(Iv[j][k], param.pk.N);
            }
        }
    }
    else
    {
        Mat<ZZ> tmp;
        tmp.SetDims(param.t, param.pk.l + 1);
        for (int i = 0; i < param.m; ++i)
        {
            for (int j = 0; j < param.t; ++j)
            {
                HSS_Input(tmp[j], param.pk, ZZ(bit(Y[i], j)));
            }
            Iy.push_back(tmp);
        }

        Iv.SetDims(param.k, param.pk.l + 1);
        for (int j = 0; j < param.k; ++j)
        {
            HSS_Input(Iv[j], param.pk, V[j]);
        }
    }
    t2 = GetTime() - t1;
    cout << "Provider Input algo time: " << t2 * 1000 << " ms\n";
    // End: Provider Input Test

    // Start: Evaluation 0 Test
    t1 = GetTime();
    vec_ZZ pc0, vv0;
    int prf_key = 0;

    Vec<ZZ> M1_0, cmp_res0;
    HSS_M1Gen(M1_0, 0, param.pk, ek0, prf_key);

    cmp_res0.SetLength(param.m);
    double tx = GetTime() - t1;
    cout << "Previous Node eval algo time: " << tx * 1000 << " ms\n";

    for (int j = 0; j < param.m; ++j)
    {
        double tt = GetTime();
        Vec<ZZ> Mxb, Mxyb, Mcb;
        HSS_Mul(Mxb, param.pk, ek0, Ix[j][0], M1_0, prf_key);
        HSS_Mul(Mxyb, param.pk, ek0, Iy[j][0], Mxb, prf_key);
        HSS_SubMemory(Mcb, param.pk, Mxb, Mxyb);

        for (int i = 1; i < param.t; ++i)
        {
            Vec<ZZ> Mcxb, Mcyb, Mcxyb;

            HSS_Mul(Mcxb, param.pk, ek0, Ix[j][i], Mcb, prf_key);
            HSS_Mul(Mcyb, param.pk, ek0, Iy[j][i], Mcb, prf_key);
            HSS_Mul(Mcxyb, param.pk, ek0, Ix[j][i], Mcyb, prf_key);
            HSS_Mul(Mxb, param.pk, ek0, Ix[j][i], M1_0, prf_key);
            HSS_Mul(Mxyb, param.pk, ek0, Iy[j][i], Mxb, prf_key);

            HSS_SubMemory(Mcb, param.pk, Mcb, Mcxb);
            HSS_SubMemory(Mcb, param.pk, Mcb, Mcyb);
            HSS_AddMemory(Mcb, param.pk, Mcb, Mcxyb);
            HSS_AddMemory(Mcb, param.pk, Mcb, Mcxyb);
            HSS_AddMemory(Mcb, param.pk, Mcb, Mxb);
            HSS_SubMemory(Mcb, param.pk, Mcb, Mxyb);
        }
        cmp_res0[j] = Mcb[0];
        tt = GetTime() - tt;
        cout << "Each Node need time: " << tt * 1000 << " ms\n";
    }

    for (int j = 0; j < param.m; ++j)
    {
        RandomBnd(cmp_res0[j], param.pk.N);
    }
    double ttt = GetTime();
    ClassificationGen(pc0, vv0, 0, param, ek0, cmp_res0, Iv, M1_0, prf_key);
    ttt = GetTime() - ttt;
    cout << "Classification algo time: " << ttt * 1000 << " ms\n";

    t2 = GetTime() - t1;
    cout << "Evaluate0 algo time: " << t2 * 1000 << " ms\n";
    // End: Evaluation 0 Test

    // Start: Evaluation 1 Test
    vec_ZZ pc1, vv1;
    prf_key = 0;
    t1 = GetTime();
    Vec<ZZ> M1_1, cmp_res1;
    HSS_M1Gen(M1_1, 1, param.pk, ek1, prf_key);

    cmp_res1.SetLength(param.m);
    for (int j = 0; j < param.m; ++j)
    {
        Vec<ZZ> Mxb, Mxyb, Mcb;
        HSS_Mul(Mxb, param.pk, ek0, Ix[j][0], M1_1, prf_key);
        HSS_Mul(Mxyb, param.pk, ek0, Iy[j][0], Mxb, prf_key);
        HSS_SubMemory(Mcb, param.pk, Mxb, Mxyb);

        for (int i = 1; i < param.t; ++i)
        {
            Vec<ZZ> Mcxb, Mcyb, Mcxyb;

            HSS_Mul(Mcxb, param.pk, ek0, Ix[j][i], Mcb, prf_key);
            HSS_Mul(Mcyb, param.pk, ek0, Iy[j][i], Mcb, prf_key);
            HSS_Mul(Mcxyb, param.pk, ek0, Ix[j][i], Mcyb, prf_key);
            HSS_Mul(Mxb, param.pk, ek0, Ix[j][i], M1_1, prf_key);
            HSS_Mul(Mxyb, param.pk, ek0, Iy[j][i], Mxb, prf_key);

            HSS_SubMemory(Mcb, param.pk, Mcb, Mcxb);
            HSS_SubMemory(Mcb, param.pk, Mcb, Mcyb);
            HSS_AddMemory(Mcb, param.pk, Mcb, Mcxyb);
            HSS_AddMemory(Mcb, param.pk, Mcb, Mcxyb);
            HSS_AddMemory(Mcb, param.pk, Mcb, Mxb);
            HSS_SubMemory(Mcb, param.pk, Mcb, Mxyb);
        }
        cmp_res1[j] = Mcb[0];
    }

    ClassificationGen(pc1, vv1, 1, param, ek1, cmp_res1, Iv, M1_1, prf_key);

    t2 = GetTime() - t1;
    cout << "Evaluate1 algo time: " << t2 * 1000 << " ms\n";
    // End: Evaluation 1 Test

    // Start: Decryption Test
    ZZ Predict_v;
    t1 = GetTime();
    ClDecryption(Predict_v, param, pc0, vv0, pc1, vv1);
    t2 = GetTime() - t1;
    cout << "Decryption algo time: " << t2 * 1000 << " ms\n";
    // End: Decryption Test

    if (V[true_idx] == Predict_v)
    {
        cout << "YES\n";
        cout << V[true_idx] << endl;
        cout << Predict_v << endl;
    }
    else
    {
        cout << "NO\n";
        cout << V[true_idx] << endl;
        cout << Predict_v << endl;
    }
}
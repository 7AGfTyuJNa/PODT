#include "TDSC20.h"

void ParaGen(TDSCPara &pk, BeaT &BT0,
             BeaT &BT1,
             BeaMaT &BMT0,
             BeaMaT &BMT1)
{

    pk.k = 2 << (pk.h - 1);
    pk.m = pk.k - 1;

    pk.P = power2_ZZ(128);

    ZZ g1, g2, g3;
    RandomBnd(g1, pk.P);
    RandomBnd(g2, pk.P);
    MulMod(g3, g1, g2, pk.P);

    ShareGen(BT0.g1, BT1.g1, pk.P, g1);
    ShareGen(BT0.g2, BT1.g2, pk.P, g2);
    ShareGen(BT0.g3, BT1.g3, pk.P, g3);

    BMT0.G1.SetDims(pk.m, pk.n);
    BMT0.G2.SetLength(pk.n);
    BMT0.G3.SetLength(pk.m);

    BMT1.G1.SetDims(pk.m, pk.n);
    BMT1.G2.SetLength(pk.n);
    BMT1.G3.SetLength(pk.m);

    Mat<ZZ> G1;
    vec_ZZ G2, G3;
    G1.SetDims(pk.m, pk.n);
    G2.SetLength(pk.n);
    G3.SetLength(pk.m);
    for (int i = 0; i < pk.m; ++i)
    {
        for (int j = 0; j < pk.n; ++j)
        {
            RandomBnd(G1[i][j], pk.P);
            ShareGen(BMT0.G1[i][j], BMT1.G1[i][j], pk.P, G1[i][j]);
        }
    }

    for (int i = 0; i < pk.n; ++i)
    {
        RandomBnd(G2[i], pk.P);
        ShareGen(BMT0.G2[i], BMT1.G2[i], pk.P, G2[i]);
    }

    MatMul(G3, pk, G1, G2);
    for (int i = 0; i < pk.m; ++i)
    {
        ShareGen(BMT0.G3[i], BMT1.G3[i], pk.P, G3[i]);
    }
}

void ClientEnc(vec_ZZ &x0, vec_ZZ &x1, TDSCPara pk, vec_ZZ x)
{
    x0.SetLength(x.length());
    x1.SetLength(x.length());
    for (int i = 0; i < x.length(); ++i)
    {
        ShareGen(x0[i], x1[i], pk.P, x[i]);
    }
}

void ProviderEnc(vec_ZZ &y0, vec_ZZ &y1, vec_ZZ &v0, vec_ZZ &v1,
                 Mat<ZZ> &mat0, Mat<ZZ> &mat1, TDSCPara pk, vec_ZZ y, vec_ZZ v, Mat<ZZ> matt)
{
    y0.SetLength(y.length());
    y1.SetLength(y.length());
    v0.SetLength(v.length());
    v1.SetLength(v.length());
    mat0.SetDims(matt.NumRows(), matt.NumCols());
    mat1.SetDims(matt.NumRows(), matt.NumCols());

    for (int i = 0; i < y.length(); ++i)
    {
        ShareGen(y0[i], y1[i], pk.P, y[i]);
    }

    for (int i = 0; i < v.length(); ++i)
    {
        ShareGen(v0[i], v1[i], pk.P, v[i]);
    }

    for (int i = 0; i < matt.NumRows(); ++i)
    {
        for (int j = 0; j < matt.NumCols(); ++j)
        {
            ShareGen(mat0[i][j], mat1[i][j], pk.P, matt[i][j]);
        }
    }
}

void InputSelection(vec_ZZ &res0, vec_ZZ &res1, Mat<ZZ> M0, Mat<ZZ> M1, TDSCPara pk,
                    BeaMaT BMT0,
                    BeaMaT BMT1,
                    vec_ZZ x0, vec_ZZ x1)
{
    // Server 0;
    Mat<ZZ> e0;
    MatSub(e0, pk, M0, BMT0.G1);
    vec_ZZ f0;
    f0.SetLength(x0.length());
    for (int i = 0; i < x0.length(); ++i)
    {
        SubMod(f0[i], x0[i], BMT0.G2[i], pk.P);
    }

    // Server 1
    Mat<ZZ> e1;
    MatSub(e1, pk, M1, BMT1.G1);
    vec_ZZ f1;
    f1.SetLength(x1.length());
    for (int i = 0; i < x1.length(); ++i)
    {
        SubMod(f1[i], x1[i], BMT1.G2[i], pk.P);
    }

    // together
    Mat<ZZ> e;
    vec_ZZ f;
    MatAdd(e, pk, e0, e1);
    f.SetLength(f0.length());
    for (int i = 0; i < x1.length(); ++i)
    {
        AddMod(f[i], f1[i], f0[i], pk.P);
    }

    // Server 0
    Vec<ZZ> tmp0;
    MatMul(res0, pk, e, BMT0.G2);
    MatMul(tmp0, pk, BMT0.G1, f);
    for (int i = 0; i < res0.length(); ++i)
    {
        AddMod(res0[i], res0[i], tmp0[i], pk.P);
        AddMod(res0[i], res0[i], BMT0.G3[i], pk.P);
    }

    // Server 1
    Vec<ZZ> tmp1;
    MatMul(res1, pk, e, BMT1.G2);
    MatMul(tmp1, pk, BMT1.G1, f);
    for (int i = 0; i < res1.length(); ++i)
    {
        AddMod(res1[i], res1[i], tmp1[i], pk.P);
        AddMod(res1[i], res1[i], BMT1.G3[i], pk.P);
    }
    MatMul(tmp1, pk, e, f);
    for (int i = 0; i < res1.length(); ++i)
    {
        AddMod(res1[i], res1[i], tmp1[i], pk.P);
    }
}

void SSCMP(ZZ &cmp_res0, ZZ &cmp_res1, TDSCPara pk, BeaT BT0,
           BeaT BT1, ZZ x0, ZZ y0, ZZ x1, ZZ y1)
{
    // Server 0
    ZZ a0;
    SubMod(a0, y0, x0, pk.P);

    vec_ZZ p;
    p.SetLength(NumBits(a0));
    for (int i = 0; i < p.length(); ++i)
    {
        p[i] = ZZ(bit(a0, i));
    }

    // Server 1
    ZZ a1;
    SubMod(a1, y1, x1, pk.P);
    vec_ZZ q;
    q.SetLength(NumBits(a1));
    for (int i = 0; i < q.length(); ++i)
    {
        q[i] = ZZ(bit(a1, i));
    }

    // together
    ZZ c_0, c_1;
    SSZZMult(c_0, c_1, pk, BT0, BT1, p[0], ZZ(0), ZZ(0), q[0]);

    for (int i = 1; i < min(NumBits(a0), NumBits(a1)) - 1; ++i)
    {
        ZZ d_0, d_1;
        SSZZMult(d_0, d_1, pk, BT0, BT1, p[i], ZZ(0), ZZ(0), q[i]);
        AddMod(d_0, d_0, ZZ(1), pk.P);
        AddMod(d_1, d_1, ZZ(1), pk.P);

        ZZ e_0, e_1;
        SSZZMult(e_0, e_1, pk, BT0, BT1, p[i], c_0, q[i], c_1);
        AddMod(e_0, e_0, ZZ(1), pk.P);
        AddMod(e_1, e_1, ZZ(1), pk.P);

        SSZZMult(c_0, c_1, pk, BT0, BT1, e_0, d_0, e_1, d_1);
        AddMod(c_0, c_0, ZZ(1), pk.P);
        AddMod(c_1, c_1, ZZ(1), pk.P);
    }

    ZZ t_1, t_2;
    AddMod(t_1, p[NumBits(a0) - 1], c_0, pk.P);
    AddMod(t_2, q[NumBits(a1) - 1], c_1, pk.P);

    ZZ tmp0, tmp1;
    SSZZMult(tmp0, tmp1, pk, BT0, BT1, t_1, ZZ(0), ZZ(0), t_2);

    cmp_res0 = (t_1 - tmp0 - tmp0) % pk.P;

    SubMod(cmp_res1, t_2, tmp1, pk.P);
    SubMod(cmp_res1, cmp_res1, tmp1, pk.P);
}

void ClassificationGen(vec_ZZ &pc0, vec_ZZ &pc1, vec_ZZ &vv0, vec_ZZ &vv1, TDSCPara pk,
                       vec_ZZ cmp_res0, vec_ZZ cmp_res1, vec_ZZ v0, vec_ZZ v1)
{
    vec_ZZ vb;
    ZZ ppc0, ppc1, v, r_0, r_1;
    int idx;

    pc0.SetLength(pk.k);
    pc1.SetLength(pk.k);
    vv0.SetLength(pk.k);
    vv1.SetLength(pk.k);

    for (int i = 0; i < pk.k; ++i)
    {
        double tt = GetTime();
        // compute pc_k
        ppc0 = 0;
        ppc1 = 0;
        idx = pk.m + i;

        while (idx)
        {
            bool lr = (idx - 1) & 1; // & 1
            idx = (idx - 1) / 2;
            if (lr)
            {
                ppc0 += (1 - cmp_res0[idx]) % pk.P;
                ppc1 += cmp_res1[idx];
            }
            else
            {
                ppc0 += cmp_res0[idx];
                ppc1 += cmp_res1[idx];
            }
        }

        // random
        int prf_key = 0;
        r_0 = PRF_ZZ(prf_key++, pk.P);
        r_1 = PRF_ZZ(prf_key++, pk.P);

        MulMod(r_1, r_1, ppc0, pk.P);
        AddMod(vv0[i], r_1, v0[i], pk.P);
        MulMod(pc0[i], r_0, ppc0, pk.P);

        prf_key = 0;
        r_0 = PRF_ZZ(prf_key++, pk.P);
        r_1 = PRF_ZZ(prf_key++, pk.P);
        MulMod(r_1, r_1, ppc1, pk.P);
        AddMod(vv1[i], r_1, v1[i], pk.P);
        MulMod(pc1[i], r_0, ppc1, pk.P);
        tt = GetTime() - tt;
        cout << "Each leaf Node time: " << tt * 1000 << " ms\n";
    }
}

void DTEvaluation(vec_ZZ &pc0, vec_ZZ &pc1, vec_ZZ &vv0, vec_ZZ &vv1, TDSCPara pk,
                  BeaT BT0,
                  BeaT BT1,
                  BeaMaT BMT0,
                  BeaMaT BMT1,
                  vec_ZZ x0, vec_ZZ x1, vec_ZZ y0, vec_ZZ y1, vec_ZZ v0, vec_ZZ v1,
                  Mat<ZZ> mat0, Mat<ZZ> mat1)
{
    vec_ZZ newx0, newx1;
    newx0.SetLength(y0.length());
    newx1.SetLength(y1.length());

    double t1 = GetTime();
    InputSelection(newx0, newx1, mat0, mat1, pk, BMT0, BMT1, x0, x1);
    t1 = GetTime() - t1;
    cout << "Pre node takes :" << t1 * 1000 << "ms\n";
    return ;

    vec_ZZ cmp_res0, cmp_res1;
    cmp_res0.SetLength(y0.length());
    cmp_res1.SetLength(y1.length());

    for (int i = 0; i < pk.m; ++i)
    {
        double tt = GetTime();
        SSCMP(cmp_res0[0], cmp_res1[1], pk, BT0, BT1, newx0[i], y0[i], newx1[i], y1[i]);
        tt = GetTime() - tt;
        cout << "Each decision node takes :" << tt * 1000 << "ms\n";
    }
    ClassificationGen(pc0, pc1, vv0, vv1, pk, cmp_res0, cmp_res1, v0, v1);
}

void Decrypt(ZZ &res, vec_ZZ pc0, vec_ZZ pc1, vec_ZZ vv0, vec_ZZ vv1, TDSCPara pk)
{
    for (int i = 0; i < pk.k; ++i)
    {
        if (pc1[i] + pc0[i] == 0)
        {
            res = vv1[i] + vv0[i];
            break;
        }
    }
}

void ShareGen(ZZ &res0, ZZ &res1, ZZ modP, ZZ x)
{
    RandomBnd(res1, modP);
    SubMod(res0, x, res1, modP);
}

void SSZZMult(ZZ &res0, ZZ &res1, TDSCPara pk, BeaT BT0,
              BeaT BT1, ZZ a0, ZZ b0, ZZ a1, ZZ b1)
{
    ZZ e0, f0;
    // Server 0
    SubMod(e0, a0, BT0.g1, pk.P);
    SubMod(f0, b0, BT0.g2, pk.P);

    // Server 1
    ZZ e1, f1;
    SubMod(e1, a1, BT1.g2, pk.P);
    SubMod(f1, b1, BT1.g2, pk.P);

    // together
    ZZ e, f;
    AddMod(e, e0, e1, pk.P);
    AddMod(f, f0, f1, pk.P);

    // Server 0
    ZZ tmp0;
    MulMod(res0, e, BT0.g2, pk.P);
    MulMod(tmp0, f, BT0.g1, pk.P);
    AddMod(res0, res0, tmp0, pk.P);
    AddMod(res0, res0, BT0.g3, pk.P);

    // Server 1
    ZZ tmp1;
    MulMod(res1, e, BT1.g2, pk.P);
    MulMod(tmp1, f, BT1.g1, pk.P);
    AddMod(res1, res1, tmp1, pk.P);
    AddMod(res1, res1, BT1.g3, pk.P);
    MulMod(tmp1, f, e, pk.P);
    AddMod(res1, res1, tmp1, pk.P);
}

void MatMul(vec_ZZ &res, TDSCPara pk, Mat<ZZ> x, vec_ZZ y)
{
    if (!res.length())
    {
        res.SetLength(x.NumRows());
    }

    ZZ tmp;
    for (int i = 0; i < x.NumRows(); ++i)
    {
        res[i] = 0;
        for (int j = 0; j < x.NumCols(); ++j)
        {
            MulMod(tmp, x[i][j], y[j], pk.P);
            AddMod(res[i], res[i], tmp, pk.P);
        }
    }
}

void MatSub(Mat<ZZ> &res, TDSCPara pk, Mat<ZZ> x, Mat<ZZ> y)
{
    res.SetDims(x.NumRows(), x.NumCols());
    for (int i = 0; i < x.NumRows(); ++i)
    {
        for (int j = 0; j < x.NumCols(); ++j)
        {
            SubMod(res[i][j], x[i][j], y[i][j], pk.P);
        }
    }
}

void MatAdd(Mat<ZZ> &res, TDSCPara pk, Mat<ZZ> x, Mat<ZZ> y)
{
    res.SetDims(x.NumRows(), x.NumCols());
    for (int i = 0; i < x.NumRows(); ++i)
    {
        for (int j = 0; j < x.NumCols(); ++j)
        {
            AddMod(res[i][j], x[i][j], y[i][j], pk.P);
        }
    }
}

void TDSC20_TIME_TEST(int depth, int N_attribute, int msgbit, int cyctimes)
{
    TDSCPara pk;
    pk.h = depth;
    pk.n = N_attribute;
    pk.t = msgbit;
    BeaT BT0;
    BeaT BT1;
    BeaMaT BMT0;
    BeaMaT BMT1;
    auto *Time = new double[cyctimes];
    double time, mean, stdev;

    // Start: Setup Test
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        ParaGen(pk, BT0, BT1, BMT0, BMT1);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Setup algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Setup Test

    // Start: Random Input
    Vec<ZZ> X, Y, V;
    X.SetLength(pk.n); // Client
    Y.SetLength(pk.m); // Provider
    V.SetLength(pk.k); // Classification
    for (int i = 0; i < pk.n; ++i)
    {
        RandomBits(X[i], pk.t);
    }
    for (int i = 0; i < pk.m; ++i)
    {
        RandomBits(Y[i], pk.t);
        RandomBits(V[i], pk.t);
    }
    RandomBits(V[pk.m], pk.t);

    Mat<ZZ> matt;
    matt.SetDims(pk.m, pk.n);
    for (int i = 0; i < pk.m; ++i)
    {
        for (int j = 0; j < pk.n; ++j)
        {
            if (i == j)
            {
                matt[i][j] = ZZ(1);
            }
            else
            {
                matt[i][j] = ZZ(0);
            }
        }
    }

    int true_idx = 0;
    while (true_idx < pk.m)
    {
        if (true_idx > pk.n)
        {
            if (X[0] > Y[true_idx])
            {
                true_idx = 2 * true_idx + 1;
            }
            else
            {
                true_idx = 2 * true_idx + 2;
            }
        }
        else
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
    }
    true_idx -= pk.m; // V[true_idx] is final result.

    // End: Random Input

    // Start: Client Input Test
    vec_ZZ Ix0, Ix1;
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        ClientEnc(Ix0, Ix1, pk, X);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Client Input algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Client Input Test

    // Start: Provider Input Test
    vec_ZZ y0, y1, v0, v1;
    Mat<ZZ> mat0, mat1;
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        ProviderEnc(y0, y1, v0, v1, mat0, mat1, pk, Y, V, matt);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Provider Input algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Provider Input Test

    // Start: Evaluation  Test
    vec_ZZ pc0, pc1, vv0, vv1;
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        DTEvaluation(pc0, pc1, vv0, vv1, pk, BT0, BT1, BMT0, BMT1, Ix0, Ix1, y0, y1, v0, v1, mat0, mat1);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Evaluate algo time: " << mean * 1000 / 2 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Evaluation Test

    // Start: Decryption Test
    ZZ Predict_v;
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        Decrypt(Predict_v, pc0, pc1, vv0, vv1, pk);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Decryption algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
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
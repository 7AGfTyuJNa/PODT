#include "ODTE-tester.h"

void ODTE_SETUP_TEST(Para &param, EK &ek0, EK &ek1, int depth, int N_attribute, int msgbit, int cyctimes, bool debug)
{
    double *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; ++i)
    {
        time = GetTime();
        Para paramTEST;
        EK ek0TEST, ek1TEST;
        paramTEST.h = depth;
        paramTEST.k = 2 << (depth - 1);
        paramTEST.m = paramTEST.k - 1;
        paramTEST.t = msgbit;
        paramTEST.n = N_attribute;
        KeyGen(paramTEST, ek0TEST, ek1TEST);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Setup algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    param.h = depth;
    param.k = 2 << (depth - 1);
    param.m = param.k - 1;
    param.t = msgbit;
    param.n = N_attribute;
    KeyGen(param, ek0, ek1);
}

void ODTE_ProviderEnc_TEST(Mat<vec_ZZ> &Iy,
                           Mat<ZZ> &Iv,
                           Mat<vec_ZZ> &Idelta,
                           const vec_ZZ &Y,
                           const vec_ZZ &V,
                           const vector<vector<int>> &delta,
                           const Para &param, int cyctimes, bool debug)
{
    double *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; ++i)
    {
        Iy.kill();
        Iv.kill();
        Idelta.kill();
        if (!debug)
        {
            time = GetTime();
            ProviderEnc(Idelta, Iy, Iv, param, Y, V, delta);
            Time[i] = GetTime() - time;
        }
        else
        {
            time = GetTime();
            RandomProviderEnc(Iy, Iv, Idelta, param);
            Time[i] = GetTime() - time;
        }
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Provider Enc time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void RandomProviderEnc(Mat<vec_ZZ> &Iy,
                       Mat<ZZ> &Iv,
                       Mat<vec_ZZ> &Idelta,
                       const Para &param)
{

    int i, j, k;
    vec_ZZ tmp;
    tmp.SetLength(param.pk.l + 1);
    Iy.SetDims(param.m, param.t);
    for (i = 0; i < param.m; ++i)
    {
        for (j = 0; j < param.t; ++j)
        {
            for (int k = 0; k < param.pk.l + 1; ++k)
            {
                RandomBnd(tmp[k], param.pk.N);
            }
            Iy[i][j] = tmp;
        }
    }

    Iv.SetDims(param.k, param.pk.l + 1);
    for (int j = 0; j < param.k; ++j)
    {
        for (int k = 0; k < param.pk.l + 1; ++k)
        {
            RandomBnd(Iv[j][k], param.pk.N);
        }
    }

    Idelta.SetDims(param.m, param.n);
    for (int i = 0; i < param.m; ++i)
    {
        for (int j = 0; j < param.n; ++j)
        {
            for (int k = 0; k < param.pk.l + 1; ++k)
            {
                RandomBnd(tmp[k], param.pk.N);
            }
            Idelta[i][j] = tmp;
        }
    }
}

void ODTE_DATA_PREPARATION(Vec<ZZ> &X, Vec<ZZ> &Y, Vec<ZZ> &V, vector<vector<int>> &delta, const Para &param)
{
    X.SetLength(param.n); // Feature vector
    Y.SetLength(param.m); // Threshold value
    V.SetLength(param.k); // Classification
    delta.assign(param.m, std::vector<int>(param.n, 0));
    GenerateMatrix(param.m, param.n, delta);
    for (int i = 0; i < param.n; ++i)
    {
        RandomBits(X[i], param.t);
    }

    for (int i = 0; i < param.m; ++i)
    {
        RandomBits(Y[i], param.t);
    }
    for (int i = 0; i < param.k; ++i)
    {
        RandomBits(V[i], param.t);
    }
}

void ODTE_FeatureSelection2_TEST(Mat<vec_ZZ> &Ix,
                                 const Mat<vec_ZZ> &Idelta,
                                 const vec_ZZ &X,
                                 const Para &param, int cyctimes, bool debug)
{
    double *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; ++i)
    {
        Ix.kill();
        time = GetTime();
        FeatureSelection2(Ix, param, X, Idelta);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Feature Selection time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void ODTE_HSSCMP_TEST(const Para &param, int b, const EK &ek0, const EK &ek1, int cyctimes)
{
    int prf_key = 0;
    Vec<vec_ZZ> Ix, Iy;
    ZZ c_b, x, y;
    RandomBits(x, param.t);
    RandomBits(y, param.t);

    vec_ZZ tmp;
    for (int i = 0; i < param.t; ++i)
    {
        HSS_Input(tmp, param.pk, ZZ(bit(x, i)));
        Ix.append(tmp);
        HSS_Input(tmp, param.pk, ZZ(bit(y, i)));
        Iy.append(tmp);
    }

    EK ekb;
    ekb = b ? ek1 : ek0;
    Vec<ZZ> M1_b;
    HSS_M1Gen(M1_b, b, param.pk, ekb, prf_key);

    double *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; ++i)
    {
        time = GetTime();
        HSSCMP(c_b, b, param, ekb, Ix, Iy, prf_key, M1_b);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "HSSCMP time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void ODTE_ClassificationGen_TEST(const Para &param, int b, const EK &ek0, const EK &ek1, int cyctimes)
{
    Vec<ZZ> cmp_resb;
    cmp_resb.SetLength(param.m);

    for (int i = 0; i < param.m; ++i)
    {
        RandomBits(cmp_resb[i], 1);
    }

    int prf_key = 0;
    EK ekb;
    ekb = b ? ek1 : ek0;
    Vec<ZZ> M1_b;
    HSS_M1Gen(M1_b, b, param.pk, ekb, prf_key);
    vec_ZZ pcb, vvb;

    vec_ZZ V;
    V.SetLength(param.k); // Classification
    for (int i = 0; i < param.k; ++i)
    {
        RandomBits(V[i], param.t);
    }
    Mat<ZZ> Iv;
    Iv.SetDims(param.k, param.pk.l + 1);
    for (int j = 0; j < param.k; ++j)
    {
        for (int k = 0; k < param.pk.l + 1; ++k)
        {
            RandomBnd(Iv[j][k], param.pk.N);
        }
    }

    double *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; ++i)
    {
        pcb.kill();
        vvb.kill();
        time = GetTime();
        ClassificationGen(pcb, vvb, b, param, ekb, cmp_resb, Iv, M1_b, prf_key);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "ClassificationGen time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void ODTE_DTEvaluation_TEST(Vec<ZZ> &pc0, Vec<ZZ> &vv0, Vec<ZZ> &pc1, Vec<ZZ> &vv1,
                            const Para &param, const EK &ek0, const EK &ek1, const Mat<vec_ZZ> &Ix, const Mat<vec_ZZ> &Iy, const Mat<ZZ> &Iv,
                            int cyctimes)
{
    double *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; ++i)
    {
        pc0.kill();
        vv0.kill();
        time = GetTime();
        DTEvaluation(pc0, vv0, 0, param, ek0, Ix, Iy, Iv);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "DTevaluation 0 time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    for (int i = 0; i < cyctimes; ++i)
    {
        pc1.kill();
        vv1.kill();
        time = GetTime();
        DTEvaluation(pc1, vv1, 1, param, ek1, Ix, Iy, Iv);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "DTevaluation 1 time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void ODTE_Decryption_TEST(ZZ &res, const Para &param, const Vec<ZZ> &pc_0, const Vec<ZZ> &vv_0, const Vec<ZZ> &pc_1, const Vec<ZZ> &vv_1, int cyctimes)
{
    double *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; ++i)
    {
        time = GetTime();
        ClDecryption(res, param, pc_0, vv_0, pc_1, vv_1);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Decryption time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void ODTE_TIME_TEST(int depth, int N_attribute, int msgbit, int cyctimes, bool debug)
{
    // Start: Setup Test
    Para param;
    EK ek0, ek1;
    ODTE_SETUP_TEST(param, ek0, ek1, depth, N_attribute, msgbit, cyctimes, debug);
    // End: Setup Test

    // Start: Random Input
    Vec<ZZ> X, Y, V;
    vector<std::vector<int>> delta;
    ODTE_DATA_PREPARATION(X, Y, V, delta, param);
    // End: Random Input

    // Start: Provider Input Test
    Mat<vec_ZZ> Iy;
    Mat<ZZ> Iv;
    Mat<vec_ZZ> Idelta;
    ODTE_ProviderEnc_TEST(Iy, Iv, Idelta, Y, V, delta, param, cyctimes, debug);
    // End: Provider Input Test

    // Start: Feature Selection Test
    Mat<vec_ZZ> Ix;
    ODTE_FeatureSelection2_TEST(Ix, Idelta, X, param, cyctimes, debug);
    // End: Client Input Test
    // if (debug)
    // {
    //     for (int i = 0; i < Y.length(); ++i) // after feature selection, #X=#Y
    //     {
    //         Mat<ZZ> tmp;
    //         tmp.SetDims(param.t, param.pk.l + 1);
    //         for (int j = 0; j < param.t; ++j)
    //         {
    //             for (int k = 0; k < param.pk.l + 1; ++k)
    //             {
    //                 RandomBnd(tmp[j][k], param.pk.N2);
    //             }
    //         }
    //         Ix.push_back(tmp);
    //     }
    // }

    // Start: HSSCMP Test
    ODTE_HSSCMP_TEST(param, 0, ek0, ek1, cyctimes);
    // End: HSSCMP Test

    // Start: ClassificationGen Test
    ODTE_ClassificationGen_TEST(param, 0, ek0, ek1, cyctimes);
    // End: ClassificationGen Test

    // Start: Evaluation Test
    vec_ZZ pc0, vv0, pc1, vv1;
    ODTE_DTEvaluation_TEST(pc0, vv0, pc1, vv1, param, ek0, ek1, Ix, Iy, Iv, cyctimes);
    // End: Evaluation Test

    // Start: Decryption Test
    ZZ Predict_v;
    ODTE_Decryption_TEST(Predict_v, param, pc0, vv0, pc1, vv1, cyctimes);
    // End: Decryption Test
}
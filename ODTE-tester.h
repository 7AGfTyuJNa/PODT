#include "ODTE.h"

void ODTE_TIME_TEST(int depth, int N_attribute, int msgbit, int cyctimes, bool debug);

void ODTE_SETUP_TEST(Para &param, EK &ek0, EK &ek1, int depth, int N_attribute, int msgbit, int cyctimes, bool debug);
void ODTE_ProviderEnc_TEST(Mat<vec_ZZ> &Iy,
                           Mat<ZZ> &Iv,
                           Mat<vec_ZZ> &Idelta,
                           const vec_ZZ &Y,
                           const vec_ZZ &V,
                           const vector<vector<int>> &delta,
                           const Para &param, int cyctimes, bool debug);

void RandomProviderEnc(Mat<vec_ZZ> &Iy,
                       Mat<ZZ> &Iv,
                       Mat<vec_ZZ> &Idelta,
                       const Para &param);
void ODTE_DATA_PREPARATION(Vec<ZZ> &X, Vec<ZZ> &Y, Vec<ZZ> &V, vector<vector<int>> &delta, const Para &param);

void ODTE_HSSCMP_TEST(const Para &param, int b, const EK &ek0, const EK &ek1, int cyctimes);
void ODTE_ClassificationGen_TEST(const Para &param, int b, const EK &ek0, const EK &ek1, int cyctimes);
void ODTE_DTEvaluation_TEST(Vec<ZZ> &pc0, Vec<ZZ> &vv0, Vec<ZZ> &pc1, Vec<ZZ> &vv1,
                            const Para &param, const EK &ek0, const EK &ek1, const Mat<vec_ZZ> &Ix, const Mat<vec_ZZ> &Iy, const Mat<ZZ> &Iv,
                            int cyctimes);

void ODTE_Decryption_TEST(ZZ &res, const Para &param, const Vec<ZZ> &pc_0, const Vec<ZZ> &vv_0, const Vec<ZZ> &pc_1, const Vec<ZZ> &vv_1, int cyctimes);
#include "ODTE-tester.h"
#include "TDSC20.h"

int main(int, char **)
{
    int depth = 3;// 3，4，8，13，17
    int msgbit = 10;
    int cyctimes = 1;
    int N_attribute = 13; // 13，15，9，13，57
    bool debug = false;

    cout << "Number of non-leaf nodes: " << power2_ZZ(depth) - 1 << endl;

    ODTE_TIME_TEST(depth, N_attribute, msgbit, cyctimes, debug);

    // TDSC20_TIME_TEST(depth, N_attribute, msgbit, cyctimes);
    return 0;
}

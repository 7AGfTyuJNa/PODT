#include "ODTE-tester.h"
#include "TDSC20.h"

int main(int, char **)
{
    int depth = 4;       // 3，4，8，13，17
    int msgbit = 16;
    int cyctimes = 1;
    int N_attribute = 15; // 13，15，9，13，57
    bool debug = false;

    cout << power2_ZZ(depth) - 1 << endl;
    ODTE_TIME_TEST(depth, msgbit, cyctimes, debug);
    // TDSC20_TIME_TEST(depth, N_attribute, msgbit, cyctimes);


    return 0;
}

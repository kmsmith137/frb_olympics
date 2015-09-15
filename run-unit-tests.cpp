#include "frb_olympics.hpp"

using namespace std;
using namespace frb_olympics;


int main(int argc, char **argv)
{
    frb_rng::run_unit_tests();
    run_cf_unit_tests();

    // FIXME add bonsai unit tests
    return 0;
}

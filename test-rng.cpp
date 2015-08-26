#include "frb_olympics.hpp"

using namespace std;
using namespace frb_olympics;


static void test_rng_save(frb_rng &r, int n)
{
    frb_rng r2 = r;

    vector<double> buf(n);
    for (int i = 0; i < n; i++)
	buf[i] = r.gaussian();
    for (int i = 0; i < n; i++)
	assert(buf[i] == r2.gaussian());

    cerr << "  test_rng_save(" << n << "): pass\n";
}

int main(int argc, char **argv)
{
    static const int n = 10000000;  // 10^8
    static frb_rng r;

    test_rng_save(r, 3);
    test_rng_save(r, 4);
    test_rng_save(r, 6);
    test_rng_save(r, 5);

    double acc = 0.0;
    double acc2 = 0.0;
    double acc4 = 0.0;

    struct timeval tv1 = get_time();

    for (int i = 0; i < n; i++) {
	double x = r.gaussian();
	acc += x;
	acc2 += x*x;
	acc4 += x*x*x*x;
    }

    struct timeval tv2 = get_time();
    double dt = 1.0e9 * time_diff(tv1,tv2) / (double)n;

    cerr << "    generated " << n << " Gaussian values in " <<  dt << " nsec/value\n"
	 << "    should be close to zero: " << (acc/n) << "\n"
	 << "    should be close to 1: " << (acc2/n) << "\n"
	 << "    should be close to 3: " << (acc4/n) << "\n";
	
    return 0;
}

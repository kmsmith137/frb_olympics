// This toy program is used to make pulse plots, as visual sanity checks

#include "frb_olympics.hpp"
#include <iomanip>

using namespace std;
using namespace frb_olympics;


int main(int argc, char **argv)
{
    if (argc != 9) {
	fprintf(stderr, "usage: write-pulse <arrival_time> <intrinsic_width> <DM> <SM> <freq_lo> <freq_hi> <dt_sample> <nsamples_per_chunk>\n");
	exit(2);
    }

    double arrival_time = xlexical_cast<double> (argv[1], "arrival_time");
    double intrinsic_width = xlexical_cast<double> (argv[2], "intrinsic_width");
    double dispersion_measure = xlexical_cast<double> (argv[3], "DM");
    double scattering_measure = xlexical_cast<double> (argv[4], "SM");
    double freq_lo = xlexical_cast<double> (argv[5], "freq_lo");
    double freq_hi = xlexical_cast<double> (argv[6], "freq_hi");
    double dt_sample = xlexical_cast<double> (argv[7], "dt_sample");
    int nsamples_per_chunk = xlexical_cast<int> (argv[8], "nsamples_per_chunk");

    // fluence=1, spectral index=0
    frb_pulse p(1.0, arrival_time, intrinsic_width, dispersion_measure, scattering_measure, 0.0);

    double t0, t1;
    p.get_endpoints(t0, t1, freq_lo, freq_hi);

    assert(t0 <= t1);
    int i0 = int(t0 / dt_sample) - 2;
    int i1 = int(t1 / dt_sample) + 3;

    int nchunks = int(i1/nsamples_per_chunk) + 2;
    vector<float> v(nchunks * nsamples_per_chunk, 0.0);

    for (int ichunk = 0; ichunk < nchunks; ichunk++)
	p.add_to_timestream(freq_lo, freq_hi, &v[ichunk*nsamples_per_chunk], nsamples_per_chunk, dt_sample, ichunk);
    
    for (int i = 0; i < (int)v.size(); i++) {
	if (i >= i0 && i < i1) {
	    t0 = dt_sample * i;
	    t1 = dt_sample * (i+1);
	    cout << setprecision(10) << t0 << " " << t1 << " " << v[i] << "\n";
	}
	else
	    assert(v[i] == 0.0);
    }

    return 0;
}

#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


// helper for frb_search_params constructor
template<typename T> inline T kv_extract(const string &filename, boost::unordered_map<string,string> &kv_pairs, const char *key)
{
    boost::unordered_map<string,string>::iterator p = kv_pairs.find(key);
    T ret;

    if (p == kv_pairs.end()) {
	stringstream s;
	s << "Fatal: expected key '" << key << "' in " << filename << "\n";
	throw runtime_error(s.str().c_str());
    }

    ret = xlexical_cast<T> (p->second, p->first);
    kv_pairs.erase(p);
    return ret;
}


frb_search_params::frb_search_params(const string &filename)
{
    boost::unordered_map<string,string> kv_pairs;
    read_kv_pairs(filename, kv_pairs);

    this->dm_min = kv_extract<double>(filename, kv_pairs, "dm_min");
    this->dm_max = kv_extract<double>(filename, kv_pairs, "dm_max");
    xassert(dm_min <= dm_max);
    xassert(dm_min >= 0.0);

    this->sm_min = kv_extract<double>(filename, kv_pairs, "sm_min");
    this->sm_max = kv_extract<double>(filename, kv_pairs, "sm_max");
    xassert(sm_min <= sm_max);
    xassert(sm_min >= 0.0);

    this->beta_min = kv_extract<double>(filename, kv_pairs, "beta_min");
    this->beta_max = kv_extract<double>(filename, kv_pairs, "beta_max");
    xassert(beta_min <= beta_max);
    
    this->width_min = kv_extract<double>(filename, kv_pairs, "width_min");
    this->width_max = kv_extract<double>(filename, kv_pairs, "width_max");
    xassert(width_min <= width_max);
    xassert(width_min >= 0.0);

    this->nchan = kv_extract<int>(filename, kv_pairs, "nchan");
    xassert(nchan > 0);

    this->band_freq_lo_MHz = kv_extract<double>(filename, kv_pairs, "freq_lo");
    this->band_freq_hi_MHz = kv_extract<double>(filename, kv_pairs, "freq_hi");
    xassert(band_freq_lo_MHz < band_freq_hi_MHz);
    xassert(band_freq_lo_MHz > 100.0);

    this->dt_sample = kv_extract<double>(filename, kv_pairs, "dt_sample");
    this->nsamples_tot = kv_extract<int>(filename, kv_pairs, "nsamples_tot");
    this->nsamples_per_chunk = kv_extract<int>(filename, kv_pairs, "nsamples_per_chunk");

    xassert(dt_sample > 0.0);
    xassert(nsamples_tot > 0.0);
    xassert(nsamples_per_chunk > 0.0);
    xassert(nsamples_tot % nsamples_per_chunk == 0);

    // Include test to ensure that buffer is large enough to accommodate pulse at largest (DM, SM)
    double t0, t1;
    this->get_allowed_arrival_times(t0, t1, this->width_max, this->dm_max, this->sm_max);

    if (t0 >= t1) {
	// Buffer is too short!
	int nmin = nsamples_tot + int((t0-t1)/dt_sample) + 1;

	stringstream s;
	s << "Fatal: buffer is too short!\n"
	  << "    nsamples_tot = " << this->nsamples_tot << "\n"
	  << "    minimum nsamples_tot = " << nmin << "\n"
	  << "    suggest using a value larger than this by a factor ~2!";

	throw runtime_error(s.str().c_str());
    }

    // if this assert fails, the searched timestream is too short to accommodate the pulse
    xassert(t0 < t1);

    this->nchunks = nsamples_tot / nsamples_per_chunk;
}


void frb_search_params::write(ostream &os) const
{
    os << "dm_min = " << dm_min << endl
       << "dm_max = " << dm_max << endl
       << "sm_min = " << sm_min << endl
       << "sm_max = " << sm_max << endl
       << "beta_min = " << beta_min << endl
       << "beta_max = " << beta_max << endl
       << "width_min = " << width_min << endl
       << "width_max = " << width_max << endl
       << "nchan = " << nchan << endl
       << "freq_lo = " << band_freq_lo_MHz << endl
       << "freq_hi = " << band_freq_hi_MHz << endl
       << "dt_sample = " << dt_sample << endl
       << "nsamples_tot = " << nsamples_tot << endl
       << "nsamples_per_chunk = " << nsamples_per_chunk << endl
       << "# nchunks = " << nchunks << endl;
}


void frb_search_params::get_allowed_arrival_times(double &t0, double &t1, double intrinsic_width, double dm, double sm) const
{
    get_pulse_endpoints(t0, t1, intrinsic_width, dm, sm, this->band_freq_lo_MHz, this->band_freq_hi_MHz);
    t0 = dt_sample - t0;
    t1 = dt_sample*(nsamples_tot-1) - t1;
}


frb_pulse frb_search_params::make_random_pulse(frb_rng &r, double fluence) const
{
    double dm = r.uniform(dm_min, dm_max);
    double sm = r.uniform(sm_min, sm_max);
    double beta = r.uniform(beta_min, beta_max);
    double width = r.uniform(width_min, width_max);

    // determine range of allowed arrival times
    double t0, t1;
    this->get_allowed_arrival_times(t0, t1, width, dm, sm);

    // if this assert fails, the searched timestream is too short to accommodate the pulse
    xassert(t0 < t1);

    // now assign arrival time
    double arrival_time = r.uniform(t0, t1);
    return frb_pulse(fluence, arrival_time, width, dm, sm, beta);    
}


void frb_search_params::simulate_noise(frb_rng &r, float *timestream) const
{
    for (int i = 0; i < nchan * nsamples_per_chunk; i++)
	timestream[i] = r.gaussian();
}


void frb_search_params::add_pulse(const frb_pulse &p, float *timestream, int ichunk) const
{
    for (int i = 0; i < nchan; i++) {
	double freq_lo_MHz = this->freq_lo_of_channel(i);
	double freq_hi_MHz = this->freq_hi_of_channel(i);
	p.add_to_timestream(freq_lo_MHz, freq_hi_MHz, &timestream[i*this->nsamples_per_chunk], this->nsamples_per_chunk, this->dt_sample, ichunk);
    }
}

double frb_search_params::get_signal_to_noise_of_pulse(const frb_pulse &pulse) const
{
    double acc = 0.0;

    for (int ichan = 0; ichan < nchan; ichan++) {
	double freq_lo_MHz = this->freq_lo_of_channel(ichan);
	double freq_hi_MHz = this->freq_hi_of_channel(ichan);
	double sn = pulse.get_signal_to_noise_in_channel(freq_lo_MHz, freq_hi_MHz, this->dt_sample);
	acc += sn*sn;
    }

    return sqrt(acc);
}

// helper for frb_search_params::make_dm_table()
static inline double F(double DM, double DM0, double DM1, double epsilon)
{
    return 2.0/epsilon * (DM1/DM0) * asinh(sqrt(DM/DM1));
}

// helper for frb_search_params::make_dm_table()
static inline double Finv(double f, double DM0, double DM1, double epsilon)
{
    return DM1 * square(sinh(0.5 * epsilon * (DM0/DM1) * f));
}

void frb_search_params::make_dm_table(vector<double> &dm_table, double epsilon) const
{
    xassert(epsilon > 0.0);

    double nu = (band_freq_lo_MHz + band_freq_hi_MHz) / 2.0;      // central frequency    
    double dnu = (band_freq_hi_MHz - band_freq_lo_MHz) / nchan;   // channel width
    double DM0 = dt_sample / 4.148806e3 / (1.0/square(band_freq_lo_MHz) - 1.0/square(band_freq_hi_MHz));
    double DM1 = 0.5 / 4.148806e3 * dt_sample * (nu*nu*nu) / dnu;

    double F0 = F(dm_min, DM0, DM1, epsilon);
    double F1 = F(dm_max, DM0, DM1, epsilon);

    int ndm = max(round_up(F1-F0), 2);
    dm_table.resize(ndm);

    for (int i = 0; i < ndm; i++) {
	double F = ((ndm-1-i)*F0 + i*F1) / (double)(ndm-1);
	dm_table[i] = Finv(F, DM0, DM1, epsilon);
    }

    xassert(fabs(dm_table[0] - dm_min) < 1.0e-3);
    xassert(fabs(dm_table[ndm-1] - dm_max) < 1.0e-3);

    dm_table[0] = dm_min;
    dm_table[ndm-1] = dm_max;
}


// a heuristic guess for how to populate the spectral index table
void frb_search_params::make_spectral_index_table(vector<double> &beta_table, int nbeta) const
{
    if (beta_max < beta_min + 1.0e-10) {
	beta_table.resize(1);
	beta_table[0] = (beta_min + beta_max) / 2.;
	return;
    }

    nbeta = max(nbeta, 1);
    beta_table.resize(nbeta);

    for (int i = 0; i < nbeta; i++)
	beta_table[i] = ((nbeta-i-0.5)*beta_min + (i+0.5)*beta_max) / nbeta;
}


// helper for frb_search_params::make_sm_table()
static inline double G(double SM, double SM0)
{
    return log(SM+3*SM0);
}

// helper for frb_search_params::make_sm_table()
static inline double Ginv(double g, double SM0)
{
    return exp(g)-3*SM0;
}

// heuristic guess for how to populate the sm table
void frb_search_params::make_sm_table(vector<double> &sm_table, int nsm) const
{
    double nu_c = 0.5 * (band_freq_lo_MHz + band_freq_hi_MHz);
    double SM0 = dt_sample / scatter_broadening(1.0, nu_c);

    double Gmin = G(sm_min, SM0);
    double Gmax = G(sm_max, SM0);

    if (Gmax < Gmin + 1.0e-5) {
	sm_table.resize(1);
	sm_table[0] = (sm_min + sm_max) / 2.0;
	return;
    }

    nsm = max(nsm, 1);
    sm_table.resize(nsm);
    
    // We always put a trial SM at sm_min, but not necessarily one at sm_max.
    sm_table[0] = sm_min;

    for (int i = 1; i < nsm; i++) {
	double g = ((2*nsm-1-2*i)*Gmin + (2*i)*Gmax) / (2*nsm-1);
	sm_table[i] = Ginv(g, SM0);
    }
}



}  // namespace frb_olympics

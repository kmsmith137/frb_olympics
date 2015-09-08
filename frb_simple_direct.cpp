#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


struct frb_simple_direct_search_algorithm : public frb_search_algorithm_base
{
    double epsilon;    // determines spacing of DM table

    vector<double> dm_table;
    int ndm;           // always equal to dm_table.size()
 
    // arrival time search range for each DM (rounded to integer timesample)
    vector<int> it0_table;   // note: always negative
    vector<int> it1_table;   // can be positive or negative

    int min_it0;
    int max_it1;

    frb_simple_direct_search_algorithm(double epsilon);
    virtual ~frb_simple_direct_search_algorithm() { }

    virtual void  search_init(const frb_search_params &p);
    virtual void  search_start();
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();
};

    
frb_simple_direct_search_algorithm::frb_simple_direct_search_algorithm(double epsilon_)
    : epsilon(epsilon_), ndm(0), min_it0(0), max_it1(0)
{
    xassert(epsilon > 0.0);

    stringstream s;
    s << "simple-direct-" << epsilon; 
    this->name = s.str();
}


void frb_simple_direct_search_algorithm::search_init(const frb_search_params &p)
{
    search_params = p;
    xassert(p.nchunks == 1);   // incremental search not implemented

    p.make_dm_table(this->dm_table, epsilon);
    this->ndm = dm_table.size();

    this->it0_table.resize(ndm);
    this->it1_table.resize(ndm);

    for (int idm = 0; idm < ndm; idm++) {
	double t0, t1;
	p.get_allowed_arrival_times(t0, t1, 0.0, dm_table[idm], 0.0);

	it0_table[idm] = round_down(t0 / p.dt_sample);
	it1_table[idm] = round_up(t1 / p.dt_sample);
    }

    // some paranoid checking
    for (int idm = 0; idm< ndm; idm++) {
	xassert(it0_table[idm] >= it0_table[ndm-1]);
	xassert(it1_table[idm] <= it1_table[0]);
    }

    this->min_it0 = it0_table[ndm-1];
    this->max_it1 = it1_table[0];

    this->debug_buffer_ndm = ndm;
    this->debug_buffer_nt = max_it1 - min_it0;
    this->search_gb = 1.0e-9 * (max_it1-min_it0) * sizeof(float);
    this->search_result = 0.0;

    cout << name << " initialized " << this->ndm << " DM channels" << endl;
}


void frb_simple_direct_search_algorithm::search_start()
{
    this->search_result = -1.0e30;
}

void frb_simple_direct_search_algorithm::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    // incremental search not supported
    xassert(ichunk == 0);

    float w = 1.0 / sqrt(search_params.nchan);

    vector<float> buf_v(max_it1-min_it0, 0.0);
    float *buf = &buf_v[0];

    for (int idm = 0; idm < ndm; idm++) {
	double dm = dm_table[idm];
	int it0 = it0_table[idm];
	int it1 = it1_table[idm];

	int nt = it1-it0;
	memset(buf, 0, nt * sizeof(buf[0]));

	for (int ichan = 0; ichan < search_params.nchan; ichan++) {
	    double nu = (search_params.freq_lo_of_channel(ichan) + search_params.freq_hi_of_channel(ichan)) / 2.;
	    int d = it0 + (int)(dispersion_delay(dm,nu) / search_params.dt_sample);   // d = offset between buf and channel timestream
	    int s = ichan * search_params.nsamples_per_chunk + d;                     // s = offset between buf and chunk

	    xassert(d >= 0);
	    xassert(d+nt <= search_params.nsamples_per_chunk);
	    //int j0 = max(-d, 0);                                          // j0 = initial offset in buf
	    //int j1 = min(search_params.nsamples_per_chunk - d, nt);       // j1 = final offset in buf

	    for (int j = 0; j < nt; j++)
		buf[j] += chunk[s+j];
	}

	for (int j = 0; j < nt; j++)
	    this->search_result = max(this->search_result, w * buf[j]);

	if (debug_buffer) {
	    int d = it0 - min_it0;  // d = offset between buf and debug timestream

	    xassert(d >= 0);
	    xassert(d+nt <= debug_buffer_nt);

	    for (int j = 0; j < nt; j++)
		debug_buffer[idm*debug_buffer_nt + d + j] = w * buf[j];
	}
    }
}

void frb_simple_direct_search_algorithm::search_end()
{
    xassert(this->search_result > -1.0e30);
}

frb_search_algorithm_base *simple_direct(double epsilon)
{
    return new frb_simple_direct_search_algorithm(epsilon);
}


}  // namespace frb_olympics

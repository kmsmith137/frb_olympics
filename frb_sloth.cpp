#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


struct frb_sloth_search_algorithm : public frb_search_algorithm_base
{
    double epsilon;   // determines spacing of DM table
    int novs;         // oversampling parameter

    int ndm;
    vector<double> dm_table;  // length-ndm

    //
    // In an incremental search, we must save partially complete
    // values between calls to start_search().
    //
    //  save_i0[idm] = first arrival time (in "narrow" samples) saved,
    //     relative to start of chunk (always negative)
    //
    //  save_n[idm] = number of arrival times saved (in "narrow" samples)
    //
    vector<int> save_i0;   // length ndm
    vector<int> save_n;    // length ndm
    int max_save_n;
    
    //
    // work_buf is an array of shape (ndm, nwork)
    // 
    // At the beginning of every call to search_start(), 
    // work_buf[idm,0] corresponds to arrival time save_i0[idm],
    // in oversampled units.
    //
    int nwork;
    vector<float> work_buf;

    //
    // weight_arr is a shape-(ndm,novs) array, containing the optimal weight
    // for the given DM channel and subsampling phase (the subsampling phase
    // is defined mod novs)
    //
    vector<double> weight_arr;

    frb_sloth_search_algorithm(const frb_search_params &p, double epsilon, int noversample);

    virtual ~frb_sloth_search_algorithm() { }

    virtual void  search_start();
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();
};


frb_sloth_search_algorithm::frb_sloth_search_algorithm(const frb_search_params &p_, double epsilon_, int noversample)
    : frb_search_algorithm_base(p_), epsilon(epsilon_), novs(noversample)
{
    xassert(epsilon > 0.0);
    xassert(noversample >= 1);
    xassert(noversample <= 16);      // surely unintentional

    // FIXME current implementation numerically unstable as DM->0!
    if (p.dm_min < 0.1)
	throw runtime_error("Currently, can't specify dm_min < 0.1 with sloth!");

    double ts = p.dt_sample;

    stringstream s;
    s << "sloth-" << epsilon << "-" << noversample;
    this->name = s.str();

    p.make_dm_table(this->dm_table, epsilon);
    this->ndm = dm_table.size();

    // Compute save_i0, save_n, max_save_n
    this->save_i0.resize(ndm);
    this->save_n.resize(ndm);
    this->max_save_n = 0;

    for (int idm = 0; idm < ndm; idm++) {
	// Min, max dispersion delay across band, in "narrow samples"
	double dt0 = novs * dispersion_delay(dm_table[idm], p.band_freq_hi_MHz) / ts;
	double dt1 = novs * dispersion_delay(dm_table[idm], p.band_freq_lo_MHz) / ts;

	// Earliest pulse arrival time which overlaps "wide" sample 0
	int i0 = round_up(-dt1);
	xassert(i0 + dt1 > -1.0e-10);
	xassert(i0 + dt1 - 1 < 1.0e-10);

	// Earliest pulse arrival time which does not overlap negatively indexed "wide" samples
	int i1 = round_up(-dt0);
	xassert(i1 + dt0 > -1.0e-10);
	xassert(i1 + dt0 - 1 < 1.0e-10);

	save_i0[idm] = i0;
	save_n[idm] = i1 - i0;
	max_save_n = max(max_save_n, i1-i0);
    }

    //
    // Compute weight_arr
    //
    // FIXME there is some cut-and-paste with search_chunk() below
    // (should factor shared code to a function)
    //
    this->weight_arr.resize(ndm * novs);
    memset(&weight_arr[0], 0, ndm * novs * sizeof(weight_arr[0]));

    for (int ichan = 0; ichan < p.nchan; ichan++) {
	double nu_lo = p.freq_lo_of_channel(ichan);
	double nu_hi = p.freq_hi_of_channel(ichan);

	for (int idm = 0; idm < ndm; idm++) {
	    // Dispersion delay, in "narrow" samples
	    double dm = dm_table[idm];
	    double dt0 = novs * dispersion_delay(dm,nu_hi)/ts;
	    double dt1 = novs * dispersion_delay(dm,nu_lo)/ts;

	    int i0 = (int)dt0;
	    int i1 = (int)dt1;

	    for (int s = 0; s < novs; s++) {
		// Consider pulse with arrival time s, so that the pulse endpoints are (dt0+s,dt1+s)
		// These are the locations of the "wide" pulse boundary, but in narrow sample units
		int j0 = ((i0+s)/novs) * novs;
		int j1 = ((i1+s)/novs) * novs;

		if (j0 == j1)
		    weight_arr[idm*novs + s] += 1.0;
		else {
		    double w = 0.0;
		    w += square((j0+novs) - (dt0+s));   // contribution from first wide sample
		    w += (j1-j0-novs) * novs;           // contribution from middle samples (note no square on novs)
		    w += square((dt1+s) - j1);          // contribution from last wide sample

		    weight_arr[idm*novs + s] += w / square(dt1-dt0);
		}
	    }
	}
    }

    for (int i = 0; i < ndm*novs; i++)
	weight_arr[i] = 1.0 / sqrt(weight_arr[i]);

    //
    // How large does the work buffer need to be?
    //
    // In each call to start_search(), the earliest arrival time which
    // overlaps sample 0 is save_i0[idm].  The earliest arrival time which
    // does not overlap sample (N-1) is (save_i1[idm] + N*novs).  Thus the
    // work buffer size is (max_save_n + N*novs).
    //
    this->nwork = max_save_n + p.nsamples_per_chunk * novs;

    this->debug_buffer_ndm = ndm;
    this->debug_buffer_nt = novs * p.nsamples_tot;   // see logic in search_chunk() below
    this->search_gb = 1.0e-9 * ndm * nwork * sizeof(float);

    cout << name << " initialized " << this->ndm << " DM channels" << endl;
}


void frb_sloth_search_algorithm::search_start()
{
    this->search_result = -1.0e30;
    this->work_buf.resize((size_t)ndm * (size_t)nwork);
    memset(&this->work_buf[0], 0, (size_t)ndm * (size_t)nwork * sizeof(float));
}


void frb_sloth_search_algorithm::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    double ts = p.dt_sample;
    int ns = p.nsamples_per_chunk;

    //
    // The cumsum array will store the cumulative sum of the data (one
    // frequency channel at a time).  Indices correspond to "narrow" samples.
    //
    // We prepend some sentinel values (zeros) at the beginning, and append
    // some sentinels (non-zero) at the end, for convenience below.
    //
    int icsum = max_save_n + 1;               // the "+1" is to accommodate roundoff error
    int ncsum = ns*novs + 2*max_save_n + 2;   // the "+2" is to accommodate roundoff error

    vector<float> cumsum_v(ncsum, 0.0);
    float *cumsum = &cumsum_v[0];

    for (int ichan = 0; ichan < p.nchan; ichan++) {
	double nu_lo = p.freq_lo_of_channel(ichan);
	double nu_hi = p.freq_hi_of_channel(ichan);

	//
	// Fill cumsum, including sentinel values
	//
	for (int j = 0; j < icsum+1; j++)
	    cumsum[j] = 0.0;

	for (int i = 0; i < ns; i++)
	    for (int j = icsum + i*novs; j < icsum + (i+1)*novs; j++)
		cumsum[j+1] = cumsum[j] + chunk[ichan*ns + i];

	for (int j = icsum + ns*novs + 1; j < ncsum; j++)
	    cumsum[j] = cumsum[icsum + ns*novs];
	     
	// 
	// Loop over trial DM's
	//
	for (int idm = 0; idm < ndm; idm++) {
	    // Number of entries in the work buffer to be formally updated (actual number will be smaller)
	    int nt = save_n[idm] + ns*novs;
	    xassert(nt <= nwork);

	    // Min/max dispersion delay in channel, in oversampled units
	    double dm = dm_table[idm];
	    double dt0 = novs * dispersion_delay(dm,nu_hi)/ts;
	    double dt1 = novs * dispersion_delay(dm,nu_lo)/ts;

	    // Split into integer/remainder
	    int i0 = (int)dt0;
	    int i1 = (int)dt1;
	    double x0 = dt0 - (double)i0;
	    double x1 = dt1 - (double)i1;
	    
	    // Shift i0,i1 to offset between work_buf and cumsum
	    i0 = (save_i0[idm] + i0 + icsum);
	    i1 = (save_i0[idm] + i1 + icsum);

	    // Check that there are enough sentinels
	    xassert(i0 >= 0);
	    xassert(i0 <= i1);
	    xassert(i1+nt+1 <= ncsum);

	    // Interpolation weights
	    double w00 = (1.0-x0) / (dt1-dt0);
	    double w01 = (x0) / (dt1-dt0);
	    double w10 = (1.0-x1) / (dt1-dt0);
	    double w11 = (x1) / (dt1-dt0);

	    float *dst = &work_buf[idm*nwork];

	    //
	    // Because the cumsum array has sentinels, we can write the inner loop over
	    // arrival times in a straightforwarad way.
	    //
	    for (int it = 0; it < nt; it++) {
		dst[it] += ( w10 * cumsum[it+i1] + w11 * cumsum[it+i1+1]      // Leading edge
			    -w00 * cumsum[it+i0] - w01 * cumsum[it+i0+1]);   // Trailing edge
	    }
	}
    }

    //
    // Update search_result and shift buffer for next call to search_chunk()
    //
    for (int idm = 0; idm < ndm; idm++) {
	int n = novs * ns;      // number of finalized entries
	int m = save_n[idm];    // number of partially computed entries

	// Offset between work_buf and weight_arr (defined mod novs)
	int s = novs - ((-save_i0[idm]) % novs);

	// Which of the finalized entries correpond to allowed arrival times?
	double t0, t1;
	p.get_allowed_arrival_times(t0, t1, 0.0, dm_table[idm], 0.0);
	
	// Convert to oversampled units, relative to start of work buffer
	t0 = novs * (t0/ts - ichunk*ns) - save_i0[idm];
	t1 = novs * (t1/ts - ichunk*ns) - save_i0[idm];
	
	// Round and clamp to [0,n) to get the search range
	int it0 = max(int(t0), 0);
	int it1 = min(int(t1)+1, n);

	for (int it = it0; it < it1; it++) {
	    float w = weight_arr[idm*novs + ((it+s) % novs)];
	    this->search_result = max(this->search_result, w * work_buf[idm*nwork+it]);
	}

	if (debug_buffer != NULL) {
	    for (int it = it0; it < it1; it++) {
		float w = weight_arr[idm*novs + ((it+s) % novs)];
		debug_buffer[idm*debug_buffer_nt + ichunk*n + it] = w * work_buf[idm*nwork + it];
	    }
	}

	// Shift buffer, in preperation for next call to search_chunk()
	memmove(&work_buf[idm*nwork], &work_buf[idm*nwork + n], m * sizeof(float));
	memset(&work_buf[idm*nwork + m], 0, n * sizeof(float));
    }
}

void frb_sloth_search_algorithm::search_end()
{
    xassert(this->search_result > -1.0e30);
    deallocate(this->work_buf);
}


frb_search_algorithm_base *sloth(const frb_search_params &p, double epsilon, int noversample)
{
    return new frb_sloth_search_algorithm(p, epsilon, noversample);
}


}  // namespace frb_olympics

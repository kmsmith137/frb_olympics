#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


// -------------------------------------------------------------------------------------------------
//
// This class searches one SM


struct frb_sloth_sm_subsearch : public frb_search_algorithm_base
{
    // parameters specified at construction
    double epsilon_d;
    double sm;                      // fiducial SM for this search
    int nbeta0;                     // nominal number of trial spectral indices, specified at construction (but see "nbeta" below)
    int nups;                       // upsampling parameter
    bool flag_strict_incremental;   // flag determining whether incremental search is "strict", see below

    // scattering delay parameter "ascatt" in each frequency channel, in _upsampled_ units
    vector<double> ascatt_table;    // length-nchan

    // dm table to be searched
    int ndm;
    vector<double> dm_table;  // length-ndm

    // spectral index table to be searched
    int nbeta;                  // actual number of trial spectral indices, returned by make_spectral_index_table()
    vector<double> beta_table;  // length-nbeta

    // spectral index weighting
    vector<double> wbeta;       // shape (nchan, nbeta)
    
    // sbuf provides a little buffering between calls to search_chunk(), needed for scattering convolution
    int nt_sbuf;         // in non-upsampled units
    vector<float> sbuf;  // shape (nchan, nt_sbuf)

    //
    // Bookkeeping for dedispersion
    //
    // In an incremental search, we must save partially complete values between chunks.
    //
    //  save_i0[idm] = first arrival time (in "narrow" samples) saved,
    //     relative to start of dedispered block (always negative)
    //
    //  save_n[idm] = number of arrival times saved (in "narrow" samples)
    //
    vector<int> save_i0;   // length ndm
    vector<int> save_n;    // length ndm
    int max_save_n;
    
    //
    // Workspace for dedispersion
    //
    // work_buf is an array of shape (ndm, nbeta, nwork)
    // 
    // The indexing here is confusing!  In every call to search_chunk(),
    // work_buf[idm,:,0] corresponds to nominal arrival time
    //
    //    t0_chunk - nt_sbuf*dt_wide + save_i0[idm]*dt_narrow
    //
    // where t0_chunk is the start time of the chunk and save_i0[] is negative
    //
    // To match work_buf[idm,:,0] to the scatter-convolved timestream
    // (which includes the nt_sbuf offset), one should simply shift by
    // the following number of narrow samples:
    //
    //    save_i0[idm] + (dispersion_delay) / dt_narrow
    //
    int nwork;
    vector<float> work_buf;

    // Final weighting
    vector<float> weight_arr;    // shape (ndm, nbeta, nups)


    // A helper function used in several places
    // Note that the indexology here is consistent with the comment above
    inline void setup_dedispersion(int ichan, int idm, int &i0, int &i1, double &w00, double &w01, double &w10, double &w11)
    {
	double nu_lo = search_params.freq_lo_of_channel(ichan);
	double nu_hi = search_params.freq_hi_of_channel(ichan);
	double dt_narrow = search_params.dt_sample / nups;

	// Min/max dispersion delay in channel, in upsampled units
	double dm = dm_table[idm];
	double dt0 = dispersion_delay(dm,nu_hi) / dt_narrow;
	double dt1 = dispersion_delay(dm,nu_lo) / dt_narrow;
	
	// Split into integer/remainder
	i0 = (int)dt0;
	i1 = (int)dt1;
	double x0 = dt0 - (double)i0;
	double x1 = dt1 - (double)i1;

	// Shift i0,i1 to offset between work_buf and upsampled array
	i0 += save_i0[idm];
	i1 += save_i0[idm];

	// Interpolation weights
	w00 = (1.0-x0) / (dt1-dt0);
	w01 = (x0) / (dt1-dt0);
	w10 = (1.0-x1) / (dt1-dt0);
	w11 = (x1) / (dt1-dt0);
    }


    frb_sloth_sm_subsearch(double epsilon_d, double sm, int nbeta, int nups, bool flag_strict_incremental);
    virtual ~frb_sloth_sm_subsearch() { }

    virtual void  search_init(const frb_search_params &p);
    virtual void  search_start(int mpi_rank_within_node);
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();
};


frb_sloth_sm_subsearch::frb_sloth_sm_subsearch(double epsilon_d_, double sm_, int nbeta_, int nups_, bool flag_strict_incremental_)
    : epsilon_d(epsilon_d_), sm(sm_), nbeta0(nbeta_), nups(nups_), flag_strict_incremental(flag_strict_incremental_),
      ndm(0), nbeta(0), nt_sbuf(0), max_save_n(0), nwork(0)
{
    xassert(epsilon_d > 0.0);
    xassert(nups >= 1);
    xassert(nups <= 16);  // upsampling more than this is surely unintentional
    xassert(sm >= 0.0);

    stringstream s;
    s << "subsloth-sm" << sm;

    if (nbeta0 > 1)
	s << "-nbeta" << nbeta0;
    if (nups > 1)
	s << "-ups" << nups;
    if (flag_strict_incremental)
	s << "-strict";

    this->name = s.str();
}


void frb_sloth_sm_subsearch::search_init(const frb_search_params &p)
{
    this->search_params = p;

    // FIXME current implementation numerically unstable as DM->0!
    if (p.dm_min < 0.1)
	throw runtime_error("Currently, can't specify dm_min < 0.1 with sloth!");

    // Initialize ascatt_table
    this->ascatt_table.resize(p.nchan);
    for (int ichan = 0; ichan < p.nchan; ichan++) {
	double nu = p.freq_mid_of_channel(ichan);
	ascatt_table[ichan] = ascatt(this->sm, nu, p.dt_sample / nups);   // note factor of nups here
    }

    // Initialize nt_sbuf
    double t0, t1;
    get_pulse_endpoints(t0, t1, 0.0, 0.0, this->sm, p.band_freq_lo_MHz, p.band_freq_hi_MHz);
    nt_sbuf = (int)((t1-t0) / p.dt_sample);
    xassert(nt_sbuf >= 0);

    // Initialize DM table
    p.make_dm_table(this->dm_table, epsilon_d);
    this->ndm = dm_table.size();

#if 0
    cout << "   dm_table = [";
    for (int idm = 0; idm < ndm; idm++) {
	if (idm % 10 == 9)
	    cout << "\n          ";
	cout << " " << dm_table[idm];
    }
    cout << " ]" << endl;
#endif

    // Initialize spectral index table
    p.make_spectral_index_table(this->beta_table, nbeta0);
    this->nbeta = beta_table.size();

    // Initialize wbeta
    this->wbeta.resize(p.nchan * nbeta);
    for (int ichan = 0; ichan < p.nchan; ichan++) {
	for (int ibeta = 0; ibeta < nbeta; ibeta++) {
	    double nu = p.freq_mid_of_channel(ichan);	    
	    double nu0 = 0.5 * (p.band_freq_lo_MHz + p.band_freq_hi_MHz);
	    wbeta[ichan*nbeta + ibeta] = pow(nu/nu0, beta_table[ibeta]);
	}
    }

    // Compute save_i0, save_n, max_save_n
    this->save_i0.resize(ndm);
    this->save_n.resize(ndm);
    this->max_save_n = 0;

    double ts = p.dt_sample;

    for (int idm = 0; idm < ndm; idm++) {
	// Min, max dispersion delay across band, in "narrow samples"
	double dt0 = nups * dispersion_delay(dm_table[idm], p.band_freq_hi_MHz) / ts;
	double dt1 = nups * dispersion_delay(dm_table[idm], p.band_freq_lo_MHz) / ts;

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
    // How large does the work buffer need to be?
    //
    // In each call to start_search(), the earliest arrival time which
    // overlaps sample 0 is save_i0[idm].  The earliest arrival time which
    // does not overlap sample (N-1) is (save_i1[idm] + N*nups).  Thus the
    // work buffer size is (max_save_n + N*nups).
    //
    this->nwork = max_save_n + p.nsamples_per_chunk * nups;

    //
    // Compute weight_arr
    //
    this->weight_arr.resize(ndm * nbeta * nups);
    memset(&weight_arr[0], 0, weight_arr.size() * sizeof(weight_arr[0]));

    int lag_dedisp = 1;   // placeholder for later
    vector<double> cf_dedisp(nups * lag_dedisp, 0.0);

    for (int ichan = 0; ichan < p.nchan; ichan++) {
	for (int idm = 0; idm < ndm; idm++) {
	    int i0, i1;
	    double w00, w01, w10, w11;
	    
	    setup_dedispersion(ichan, idm, i0, i1, w00, w01, w10, w11);

	    // shift i0,i1 by multiple of nups, so that 0 <= i0 < nups
	    int shift = (i0 >= 0) ? (-(i0/nups)) : ((-i0+nups-1)/nups);
	    i0 += shift * nups;
	    i1 += shift * nups;

	    xassert(i0 >= 0);
	    xassert(i0 < nups);
	    xassert(i0 <= i1);

	    int nconv = i1+1;
	    vector<double> wconv(nconv, 0.0);

	    // a little ugly, interface could be improved here
	    wconv[i0] -= w01;
	    wconv[i1] -= w10;
	    for (int i = i0; i <= i1; i++)
		wconv[i] += w00 + w01;     // equal to 1/(dt1-dt0) in notation from setup_dedispersion()

	    int lag_data = lag_dedisp + nconv - 1;
	    vector<double> cf_data(nups * lag_data, 0.0);

	    init_cf(&cf_data[0], nups, lag_data, ascatt_table[ichan]);
	    convolve_cf(&cf_dedisp[0], nups, lag_dedisp, lag_data, nconv, &wconv[0], &cf_data[0]);

	    for (int ibeta = 0; ibeta < nbeta; ibeta++) {
		double ww = square(wbeta[ichan*nbeta+ibeta]);
		for (int s = 0; s < nups; s++) {
		    xassert(idm*nbeta*nups + ibeta*nups + s < (int)weight_arr.size());
		    weight_arr[idm*nbeta*nups + ibeta*nups + s] += ww * cf_dedisp[s*lag_dedisp];
		}
	    }
	}
    }

    for (unsigned int i = 0; i < weight_arr.size(); i++)
	weight_arr[i] = 1.0 / sqrt(weight_arr[i]);

#if 0
    for (int idm = 0; idm < ndm; idm += 50) {
	int ibeta = nbeta/2;
	cout << "  idm=" << idm << " ibeta=" << ibeta << ": weight_arr=[";
	for (int s = 0; s < nups; s++) 
	    cout << " " << weight_arr[idm*nbeta*nups + ibeta*nups + s];
	cout << " ]\n";
    }
#endif

    this->debug_buffer_ndm = ndm;
    this->debug_buffer_nt = nups * p.nsamples_tot;   // see logic in search_chunk() below

    this->search_gb = 0.0;
    this->search_gb += 1.0e-9 * p.nchan * nt_sbuf * sizeof(float);     // sbuf
    this->search_gb += 1.0e-9 * ndm * nbeta * nwork * sizeof(float);   // work_buf

    cout << "  " << this->name << ": initialized, ndm="<< ndm << ", beta = [";
    for (int ibeta = 0; ibeta < nbeta; ibeta++)
	cout << " " << beta_table[ibeta];
    cout << " ]" << endl;
}


void frb_sloth_sm_subsearch::search_start(int mpi_rank_within_node)
{
    int nchan = this->search_params.nchan;

    // Initialize sbuf
    this->sbuf.resize((size_t)nchan * (size_t)nt_sbuf);
    memset(&sbuf[0], 0, sbuf.size() * sizeof(sbuf[0]));

    // Initialize work_buf
    this->work_buf.resize((size_t)ndm * (size_t)nbeta * (size_t)nwork);
    memset(&work_buf[0], 0, work_buf.size() * sizeof(work_buf[0]));

    this->search_result = -1.0e30;
}


void frb_sloth_sm_subsearch::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    double dt_narrow = search_params.dt_sample / nups;

    int nt_chunk = search_params.nsamples_per_chunk;
    int nt_avail = nt_sbuf + nt_chunk;

    vector<float> rdata_v(nt_avail, 0.0);           // raw data (neither upsampled nor scatter-convolved) from one channel
    vector<float> cdata_v(nt_avail * nups, 0.0);    // upsampled, scatter-convolved data from one channel
    vector<float> tbuf_v(nwork, 0.0);               // temporary buffer for interpolation

    float *rdata = &rdata_v[0];
    float *cdata = &cdata_v[0];
    float *tbuf = &tbuf_v[0];

    //
    // The cumsum array will store the cumulative sum of "cdata".
    //
    // We prepend some sentinel values (zeros) at the beginning, and append
    // some sentinels (non-zero) at the end, for convenience below.
    //
    //  icsum = number of leading sentinels
    //  ncsum = total size of cumsum buf, including sentinels on both ends
    //
    int icsum = max_save_n + 1;                     // the "+1" is to accommodate roundoff error
    int ncsum = nt_chunk*nups + 2*max_save_n + 2;   // the "+2" is to accommodate roundoff error

    vector<float> cumsum_v(ncsum, 0.0);
    float *cumsum = &cumsum_v[0];

    for (int ichan = 0; ichan < search_params.nchan; ichan++) {
	double a = ascatt_table[ichan];

	// Fill rdata: raw data (neither upsampled nor scatter-convolved) from one channel
	memcpy(&rdata[0], &sbuf[ichan*nt_sbuf], nt_sbuf * sizeof(rdata[0]));
	memcpy(&rdata[nt_sbuf], &chunk[ichan*nt_chunk], nt_chunk * sizeof(rdata[0]));

	// It's convenient to update the sbuf here (note that sbuf for this channel won't be used again)
	memcpy(&sbuf[ichan*nt_sbuf], &rdata[nt_chunk], nt_sbuf * sizeof(sbuf[0]));
	
	//
	// Fill cdata = upsampled, scattering-convolved data from a single channel
	//
	for (int i = 0; i < nt_avail; i++)
	    for (int j = 0; j < nups; j++)
		cdata[i*nups+j] = rdata[i];
	
	for (int i = (nt_avail*nups)-2; i >= 0; i--)
	    cdata[i] = (1-a)*cdata[i] + a*cdata[i+1];

	if (flag_strict_incremental && (a > 0.0)) {
	    xassert(nt_sbuf > 0);

	    double an = pow(a, nt_sbuf*nups);
	    for (int i = 0; i < nt_chunk*nups; i++)
		cdata[i] = (cdata[i] - an*cdata[i+nt_sbuf*nups]) / (1-an);
	}

	//
	// Fill cumsum, including sentinel values
	//
	for (int j = 0; j < icsum+1; j++)
	    cumsum[j] = 0.0;

	for (int j = 0; j < nt_chunk*nups; j++)
	    cumsum[j+icsum+1] = cumsum[j+icsum] + cdata[j];

	int n0 = icsum + nt_chunk*nups;
	for (int j = n0+1; j < ncsum; j++)
	    cumsum[j] = cumsum[n0];

	// 
	// Main dedispersion loop over trial DM's
	//
	for (int idm = 0; idm < ndm; idm++) {
	    int i0, i1;
	    double w00, w01, w10, w11;

	    setup_dedispersion(ichan, idm, i0, i1, w00, w01, w10, w11);
	    i0 += icsum;
	    i1 += icsum;

	    // Number of entries in the work buffer to be formally updated
	    int nt = save_n[idm] + nt_chunk*nups;

	    xassert(i0 >= 0);
	    xassert(i0 <= i1);
	    xassert(i1+nt+1 <= ncsum);
	    xassert(nt <= nwork);

	    //
	    // Because the cumsum array has sentinels, we can write the inner loops over
	    // arrival times in a straightforwarad way.
	    //
	    for (int it = 0; it < nt; it++) {
		tbuf[it] = ( w10 * cumsum[it+i1] + w11 * cumsum[it+i1+1]     // Leading edge
			    -w00 * cumsum[it+i0] - w01 * cumsum[it+i0+1]);   // Trailing edge

	    }

	    for (int ibeta = 0; ibeta < nbeta; ibeta++) {
		float *dst = &work_buf[idm*nbeta*nwork + ibeta*nwork];
		double w = wbeta[ichan*nbeta + ibeta];

		for (int it = 0; it < nt; it++)
		    dst[it] += w * tbuf[it];
	    }
	}
    }

    //
    // Update search_result and shift buffer for next call to search_chunk()
    //
    for (int idm = 0; idm < ndm; idm++) {
	int n = nups * nt_chunk;      // number of finalized entries
	int m = save_n[idm];          // number of partially computed entries

	double t0, t1;
	search_params.get_allowed_arrival_times(t0, t1, 0.0, dm_table[idm], 0.0);
	
	// Convert allowed range to upsampled units, relative to start of work buffer
	double t0_chunk = (ichunk*nt_chunk - nt_sbuf) * nups;   // start of work buffer, note sbuf correction
	t0 = (t0/dt_narrow) - t0_chunk - save_i0[idm];
	t1 = (t1/dt_narrow) - t0_chunk - save_i0[idm];
	
	// Round and clamp to [0,n) to get the search range
	int it0 = max(int(t0), 0);
	int it1 = min(int(t1)+1, n);

	for (int ibeta = 0; ibeta < nbeta; ibeta++) {
	    for (int it = it0; it < it1; it++) {
		float w = weight_arr[idm*nbeta*nups + ibeta*nups + (it % nups)];
		float v = w * work_buf[idm*nbeta*nwork + ibeta*nwork + it];
		this->search_result = max(this->search_result, v);
	    }
	}

	if (debug_buffer != NULL) {
	    int ibeta = nbeta/2;
	    for (int it = it0; it < it1; it++) {
		float w = weight_arr[idm*nbeta*nups + ibeta*nups + (it % nups)];
		float v = w * work_buf[idm*nbeta*nwork + ibeta*nwork + it];
	        debug_buffer[idm*debug_buffer_nt + ichunk*n + it] = v;
	    }
	}

	// Shift buffer, in preperation for next call to search_chunk()
	for (int ibeta = 0; ibeta < nbeta; ibeta++) {
	    float *p = &work_buf[idm*nbeta*nwork + ibeta*nwork];
	    memmove(p, p+n, m * sizeof(float));
	    memset(p+m, 0, n * sizeof(float));
	}
    }
}

void frb_sloth_sm_subsearch::search_end()
{
    xassert(this->search_result > -1.0e30);
    deallocate(this->sbuf);
    deallocate(this->work_buf);
}


// -------------------------------------------------------------------------------------------------


struct frb_sloth_search_algorithm : public frb_search_algorithm_base
{
    // Parameters specified at construction
    double epsilon_d;
    int nsm0;      // nominal number of trial scattering measures, specified at construction (see "nsm" below)
    int nbeta0;    // nominal number of trial spectral indices, specified at construction 
    int nups;
    bool flag_strict_incremental;

    int nsm;       // actual number of trial scattering measures, returned by make_sm_table()
    vector<double> sm_table;                            // length-nsm
    vector<frb_sloth_sm_subsearch *> subsearch_table;   // length-nsm

    frb_sloth_search_algorithm(double epsilon_d, int nsm, int nbeta, int nups, bool flag_strict_incremental);
    virtual ~frb_sloth_search_algorithm();

    virtual void  search_init(const frb_search_params &p);
    virtual void  search_start(int mpi_rank_within_node);
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();
};


frb_sloth_search_algorithm::frb_sloth_search_algorithm(double epsilon_d_, int nsm_, int nbeta_, int nups_, bool flag_strict_incremental_)
    : epsilon_d(epsilon_d_), nsm0(nsm_), nbeta0(nbeta_), nups(nups_), flag_strict_incremental(flag_strict_incremental_), nsm(0)
{
    xassert(epsilon_d > 0.0);
    xassert(nups >= 1);
    xassert(nups <= 16);      // upsampling more than this is surely unintentional

    stringstream s;
    s << "sloth";

    if (nsm0 > 1)
	s << "-nsm" << nsm0;
    if (nbeta0 > 0)
	s << "-beta" << nbeta0;
    if (nups > 1)
	s << "-ups" << nups;
    if (flag_strict_incremental)
	s << "-strict";

    this->name = s.str();
}


void frb_sloth_search_algorithm::search_init(const frb_search_params &p)
{
    this->search_params = p;

    p.make_sm_table(this->sm_table, this->nsm0);
    this->nsm = sm_table.size();

    cout << name << ": initializing, nsm=" << nsm << endl;

    this->subsearch_table.resize(nsm);
    for (int ism = 0; ism < nsm; ism++) {
	this->subsearch_table[ism] = new frb_sloth_sm_subsearch(epsilon_d, sm_table[ism], nbeta0, nups, flag_strict_incremental);
	this->subsearch_table[ism]->search_init(p);
    }

    this->debug_buffer_ndm = subsearch_table[nsm/2]->debug_buffer_ndm;
    this->debug_buffer_nt = subsearch_table[nsm/2]->debug_buffer_nt;
    this->search_gb = 0.0;

    for (int ism = 0; ism < nsm; ism++)
	search_gb += subsearch_table[ism]->search_gb;

    cout << name << ": initialized" << endl;
}


void frb_sloth_search_algorithm::search_start(int mpi_rank_within_node)
{
    for (int ism = 0; ism < nsm; ism++)
	this->subsearch_table[ism]->search_start(mpi_rank_within_node);

    this->search_result = -1.0e30;
}


void frb_sloth_search_algorithm::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    for (int ism = 0; ism < nsm; ism++) {
	float *p = (ism == nsm/2) ? debug_buffer : NULL;
	this->subsearch_table[ism]->search_chunk(chunk, ichunk, p);
    }

    this->search_result = subsearch_table[0]->search_result;
    for (int ism = 1; ism < nsm; ism++)
	this->search_result = max(this->search_result, this->subsearch_table[ism]->search_result);
}


void frb_sloth_search_algorithm::search_end()
{
    for (int ism = 0; ism < nsm; ism++)
	this->subsearch_table[ism]->search_end();
}


frb_sloth_search_algorithm::~frb_sloth_search_algorithm()
{
    for (int ism = 0; ism < nsm; ism++) {
	delete subsearch_table[ism];
	subsearch_table[ism] = NULL;
    }
}


// -------------------------------------------------------------------------------------------------


frb_search_algorithm_base *sloth(double epsilon_d, int nsm, int nbeta, int nups, bool strict_incremental)
{
    return new frb_sloth_search_algorithm(epsilon_d, nsm, nbeta, nups, strict_incremental);
}

frb_search_algorithm_base *sloth_sm_subsearch(double epsilon_d, double sm, int nbeta, int nups, bool strict_incremental)
{
    return new frb_sloth_sm_subsearch(epsilon_d, sm, nbeta, nups, strict_incremental);
}


}  // namespace frb_olympics

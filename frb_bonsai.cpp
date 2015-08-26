#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


// -------------------------------------------------------------------------------------------------
//
// Low-level kernels for the tree transform


static inline int bit_reverse(int i, int depth)
{
    int ret = 0;

    while (i > 0) {
	ret += (1 << (depth-1)) / (i & ~(i-1));
	i = i & (i-1);
    }

    return ret;
}


//
// ts1[i] <- ts1[i] + ts2[i]          0 <= i < nt
// ts2[i] <- ts1[i+s] + ts2[i+s+1]    0 <= s < nt-s-1
//
static inline void tree_merge(float *ts1, float *ts2, int s, int nt, float *scratch)
{
    for (int i = 0; i < nt-s-1; i++)
	scratch[i] = ts1[i+s] + ts2[i+s+1];
    for (int i = 0; i < nt; i++)
	ts1[i] += ts2[i];
    for (int i = 0; i < nt-s-1; i++)
	ts2[i] = scratch[i];
}


static void do_tree_transform(float *buf, int depth, int nt)
{
    vector<float> scratch(nt, 0.0);

    for (int iter = 0; iter < depth; iter++) {
	int bs = 1 << iter;              // block size
	int np = 1 << (depth-iter-1);    // number of block pairs to be merged

	for (int ip = 0; ip < np; ip++)
	    for (int j = 0; j < bs; j++)
		tree_merge(&buf[((2*ip)*bs+j)*nt], &buf[((2*ip+1)*bs+j)*nt], ip, nt, &scratch[0]);
    }
}


// Fills the 'dst' array (length 2^depth).  The scratch array should have size 2^(depth-1).
static void get_tree_offsets(int *dst, int depth, int idm, int *scratch)
{
    assert(idm >= 0 && idm < (1 << depth));

    if (depth == 0) {
	dst[0] = 0;
	return;
    }

    int high_bit = 1 << (depth-1);
    get_tree_offsets(scratch, depth-1, idm & (high_bit-1), dst);
    
    if ((idm & high_bit) == 0) {
	for (int i = 0; i < high_bit; i++)
	    dst[2*i] = dst[2*i+1] = scratch[i];
    }
    else {
	for (int i = 0; i < high_bit; i++) {
	    dst[2*i] = i + scratch[i];
	    dst[2*i+1] = i + scratch[i] + 1;
	}
    }
}



//
// @offsets = output array from get_tree_offsets() above
// It definitely suffices for @scratch to have size ndm, although this is probably overkill!
//
// FIXME: a little slow (written for code clarity, not speed)
//
static inline double w2(double a, double b, int ia, int ib, int ndm, int s, int nups, const int *offsets, double *scratch)
{
    assert(ia >= 0);
    assert(ia <= ib);
    assert(ib < ndm);
    assert(s >= 0);
    assert(s < nups);

    if (ia == ib)
	return 1.0;

    int bin0 = (offsets[ia]+s) / nups;
    int bin1 = (offsets[ib]+s) / nups + 1;
    memset(scratch, 0, (bin1-bin0) * sizeof(scratch[0]));

    for (int i = ia; i <= ib; i++) {
	double ii = (double)i;
	int bin = (offsets[i]+s) / nups;
	scratch[bin-bin0] += (min(ii+1,b) - max(ii,a)) / (b-a);
    }

    double ret = 0.0;
    for (int i = 0; i < (bin1-bin0); i++)
	ret += square(scratch[i]);

    return ret;
}


// -------------------------------------------------------------------------------------------------


struct frb_bonsai_search_algorithm : public frb_search_algorithm_base
{
    int depth;
    int ntree;    // always equals 2^depth
    int nups;     // upsampling factor

    double dm1_index;   // converts DM -> transform index
    double dm1_delay;   // converts DM -> delay in lowest lambda^2-channel

    // arrays of length nchan which keep track of the channel -> lambda^2 mapping
    int nchan;
    vector<double> chan_lo_l2;
    vector<double> chan_hi_l2;
    vector<int> chan_lo_il2;
    vector<int> chan_hi_il2;

    // shape-(ntree,nups) array containing relative weight for each (DM, subsample_phase) pair
    vector<double> weight_arr;

    frb_bonsai_search_algorithm(const frb_search_params &p, int depth, int nups);

    virtual ~frb_bonsai_search_algorithm() { }

    virtual void  search_start();
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();

    void remap_channels(float *buf, const float *chunk) const;
    void run_unit_tests() const;

    static frb_search_algorithm_base::ptr_t create(const vector<string> &tokens, const frb_search_params &p);
};


frb_bonsai_search_algorithm::frb_bonsai_search_algorithm(const frb_search_params &p, int depth_, int nups_)
    : frb_search_algorithm_base(p), depth(depth_), ntree(1 << depth_), nups(nups_)
{
    assert(depth >= 1);
    assert(p.nchunks == 1);  // incremental search not implemented yet
    assert(nups >= 1);

    stringstream s;
    s << "bonsai-" << ntree;

    this->name = s.str();
    this->debug_buffer_ndm = ntree;
    this->debug_buffer_nt = p.nsamples_tot * nups;
    this->search_gb = 1.0e-9 * (double)(ntree) * (double)(p.nsamples_tot) * (double)nups * sizeof(float);

    double t0 = dispersion_delay(1.0, p.band_freq_hi_MHz);
    double t1 = dispersion_delay(1.0, p.band_freq_lo_MHz);

    this->dm1_index = (t1-t0)/(p.dt_sample/nups) * (ntree-1)/(double)ntree;
    this->dm1_delay = t0 + (t1-t0)/(2*ntree);

    double tree_dm_max = (ntree-1) / dm1_index;
    if (p.dm_max > tree_dm_max) {
	cerr << "bonsai: max DM in paramfile (=" << p.dm_max << ") exceeds max DM of tree (=" << tree_dm_max << ")\n";
	exit(1);
    }

    // 
    // Initialize channel->lambda2 mapping
    //
    this->nchan = p.nchan;
    this->chan_lo_l2.resize(nchan);
    this->chan_hi_l2.resize(nchan);
    this->chan_lo_il2.resize(nchan);
    this->chan_hi_il2.resize(nchan);

    for (int ichan = 0; ichan < nchan; ichan++) {
	double tt0 = dispersion_delay(1.0, p.freq_hi_of_channel(ichan));
	double tt1 = dispersion_delay(1.0, p.freq_lo_of_channel(ichan));

	double a = ntree * (tt0 - t0) / (t1 - t0);
	double b = ntree * (tt1 - t0) / (t1 - t0);

	assert(a > -1.0e-10);
	assert(a < b);
	assert(b < ntree + 1.0e-10);

	int ia = max((int)a, 0);
	int ib = min((int)b, ntree-1);
	assert(ia <= ib);

	chan_lo_l2[ichan] = a;
	chan_hi_l2[ichan] = b;
	chan_lo_il2[ichan] = ia;
	chan_hi_il2[ichan] = ib;
    }

    //
    // Initialize (DM, undersample_phase) weight array
    //
    this->weight_arr.resize(ntree * nups);
    memset(&weight_arr[0], 0, weight_arr.size() * sizeof(weight_arr[0]));

    vector<int> offsets(ntree);
    vector<int> scratch(ntree/2);
    vector<double> scratch2(ntree);

    for (int idm = 0; idm < ntree; idm++) {
	get_tree_offsets(&offsets[0], depth, idm, &scratch[0]);

	for (int s = 0; s < nups; s++) {
	    double w = 0.0;
	    for (int ichan = 0; ichan < nchan; ichan++) {
		w += w2(chan_lo_l2[ichan], chan_hi_l2[ichan], 
			chan_lo_il2[ichan], chan_hi_il2[ichan], 
			ntree, s, nups, &offsets[0], &scratch2[0]);
	    }

	    assert(w > 0.0);
	    weight_arr[idm*nups + s] = 1.0 / sqrt(w);
	}
    }
    
    //
    // Print "occupancy" of the tree (number of lambda^2 channels associated to 
    // each frequency channel).  I just like to look at this once in a while!
    //
    vector<int> occupancy;

    for (int ichan = 0; ichan < nchan; ichan++) {
	int j = chan_hi_il2[ichan] - chan_lo_il2[ichan];
	assert(j >= 0);

	while ((int)occupancy.size() <= j)
	    occupancy.push_back(0);

	occupancy[j]++;
    }

    cout << name << " initialized, tree_dm_max=" << tree_dm_max << ", occupancy = [";
    for (int i = 0; i < (int)occupancy.size(); i++)
	cout << " " << occupancy[i];
    cout << " ]" << endl;

    if (depth <= 5) {
	cout << name << " depth <= 5, running unit tests\n";
	this->run_unit_tests();
    }
}

void frb_bonsai_search_algorithm::remap_channels(float *buf, const float *chunk) const
{
    int nt_wide = p.nsamples_per_chunk;
    int nt_narrow = p.nsamples_per_chunk * nups;

    memset(buf, 0, ntree * nt_narrow * sizeof(float));

    for (int ichan = 0; ichan < nchan; ichan++) {
	double a = chan_lo_l2[ichan];
	double b = chan_hi_l2[ichan];
	int ia = chan_lo_il2[ichan];
	int ib = chan_hi_il2[ichan];

	for (int itree = ia; itree <= ib; itree++) {
	    double ii = (double)itree;  // std::max gets confused by types otherwise
	    double wt = (min(b,ii+1) - max(a,ii)) / (b-a);
	    
	    for (int it = 0; it < nt_wide; it++)
		for (int it2 = it*nups; it2 < (it+1)*nups; it2++)
		    buf[itree*nt_narrow + it2] += wt * chunk[ichan*nt_wide + it];
	}
    }
}

void frb_bonsai_search_algorithm::search_start()
{
    this->search_result = -1.0e30;
}

void frb_bonsai_search_algorithm::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    // incremental search not supported yet
    assert(ichunk == 0);

    int nt_narrow = p.nsamples_per_chunk * nups;
    double dt_narrow = p.dt_sample / nups;

    vector<float> buf(ntree*nt_narrow, 0.0);

    remap_channels(&buf[0], chunk);
    do_tree_transform(&buf[0], depth, nt_narrow);

    // range of DM indices to be searched
    int idm0 = (int)(p.dm_min * dm1_index);      // round down
    int idm1 = (int)(p.dm_max * dm1_index) + 1;  // round up

    // loop over trial DM's, keeping in mind that the tree uses bit-reversed indexing
    for (int itree = 0; itree < ntree; itree++) {
	int idm = bit_reverse(itree, depth);
	if ((idm < idm0) || (idm >= idm1))
	    continue;

	double dm = idm / dm1_index;

	double t0, t1;
	p.get_allowed_arrival_times(t0, t1, 0.0, dm, 0.0);

	// convert to buffer indices (still floating-point though)
	t0 = (t0 + dm * dm1_delay) / dt_narrow;
	t1 = (t1 + dm * dm1_delay) / dt_narrow;

	// clamp and round
	int it0 = max((int)(t0), 0);
	int it1 = min((int)(t1+2), nt_narrow-idm);

	for (int it = it0; it < it1; it++) {
	    float w = weight_arr[idm*nups + (it%nups)];
	    this->search_result = max(this->search_result, w * buf[itree*nt_narrow+it]);
	}

	if (debug_buffer != NULL) {
	    for (int it = it0; it < it1; it++) {
		float w = weight_arr[idm*nups + (it%nups)];
		debug_buffer[idm*debug_buffer_nt + it] = w * buf[itree*nt_narrow + it];
	    }
	}
    }
}

void frb_bonsai_search_algorithm::search_end()
{
    assert(this->search_result > -1.0e30);
}

void frb_bonsai_search_algorithm::run_unit_tests() const
{
    const double a = 31.3849;
    const double b = 17.3218;

    int nt = ntree + 10;

    vector<int> scratch(ntree);
    vector<int> all_offsets(ntree*ntree);
    vector<float> buf(ntree*nt);
    vector<float> scratch2(nt);

    for (int idm = 0; idm < ntree; idm++)
	get_tree_offsets(&all_offsets[idm*ntree], depth, idm, &scratch[0]);

    for (int ichan = 0; ichan < ntree; ichan++) {
	memset(&buf[0], 0, buf.size() * sizeof(buf[0]));
	for (int it = 0; it < nt; it++)
	    buf[ichan*nt + it] = a*ichan + b*it;

	do_tree_transform(&buf[0], depth, nt);

	for (int itree = 0; itree < ntree; itree++) {
	    int idm = bit_reverse(itree, depth);
	    int delay = all_offsets[idm*ntree + ichan];

	    for (int it = 0; it < nt-idm; it++) {
		float expected = a*ichan + b*(it+delay);
		if (fabs(buf[itree*nt+it] - expected) < 1.0e-6)
		    continue;

		cout << "test_offsets FAILED: depth=" << depth << " idm=" << idm << " ichan=" << ichan
		     << " it=" << it << " delay=" << delay << " expected=" << expected 
		     << " actual=" << buf[itree*nt+it] << endl;

		exit(1);
	    }
	}
    }

    cout << "    test_offsets: pass (depth=" << depth << ")" << endl;
    
    // And now a ridiculously slow way to compute the weight_arr

    int nt_wide = p.nsamples_per_chunk;
    int nt_narrow = nt_wide * nups;

    vector<float> chunk(nchan * nt_wide);
    vector<float> buf2(ntree * nt_narrow);
    vector<double> acc(ntree * nups, 0.0);

    for (int ichan = 0; ichan < nchan; ichan++) {
	int max_samp = (chan_hi_il2[ichan] + nups - 1) / nups;	
	assert(max_samp < nt_wide);

	for (int isamp = 0; isamp <= max_samp; isamp++) {
	    memset(&chunk[0], 0, chunk.size() * sizeof(chunk[0]));
	    chunk[ichan*nt_wide + isamp] = 1.0;

	    this->remap_channels(&buf2[0], &chunk[0]);
	    do_tree_transform(&buf2[0], depth, nt_narrow);

	    for (int itree = 0; itree < ntree; itree++) {
		int idm = bit_reverse(itree, depth);

		for (int s = 0; s < nups; s++) 
		    acc[idm*nups + s] += square(buf2[itree*nt_narrow + s]);
	    }
	}
    }

    for (int idm = 0; idm < ntree; idm++) {
	for (int s = 0; s < nups; s++) {
	    assert(acc[idm*nups+s] > 0.0);
	    
	    double w1 = weight_arr[idm*nups+s];
	    double w2 = 1.0 / sqrt(acc[idm*nups+s]);

	    if (fabs(w1-w2) < 1.0e-3)
		continue;

	    cout << "test_weights FAILED: idm=" << idm << ", s=" << s << ": w1=" << w1 << ", w2=" << w2 << ", diff=" << fabs(w1-w2) << endl;
	    exit(1);
	}
    }

    cout << "    test_weights: pass (depth=" << depth << ", nchan=" << nchan << ", nups=" << nups << ")" << endl;
}

// Registry boilerplate follows

static const char *usage = "    bonsai <ntree> [upsample NN]\n";

// static member function
frb_search_algorithm_base::ptr_t frb_bonsai_search_algorithm::create(const vector<string> &tokens_, const frb_search_params &p)
{
    vector<string> tokens = tokens_;
    int nups = 1;

    for (int i = 1; i < (int)tokens.size() - 1; i++) {
	if (tokens[i] != string("upsample"))
	    continue;

	nups = xlexical_cast<int> (tokens[i+1], "bonsai nupsample");
	assert(nups >= 1);

	tokens.erase(tokens.begin()+i, tokens.begin()+i+2);
	break;
    }

    if (tokens.size() != 1)
	return frb_search_algorithm_base::ptr_t();

    int ntree = xlexical_cast<int> (tokens[0], "bonsai ntree");

    if (!is_power_of_two(ntree)) {
	cerr << "bonsai: expected ntree to be a power of two (got ntree=" << ntree << ")\n";
	exit(1);
    }

    int depth = integer_log2(ntree);
    return boost::make_shared<frb_bonsai_search_algorithm>(p, depth, nups);
}

// register tree algorithm at library load time
namespace {
    struct _initializer {
	_initializer()
	{
	    frb_search_algorithm_base::register_algorithm("bonsai", frb_bonsai_search_algorithm::create, usage);
	}
    } _init;
}


}  // namespace frb_olympics

#include "frb_olympics.hpp"
#include <jstree.h>

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


struct frb_simple_tree_search_algorithm : public frb_search_algorithm_base
{
    int depth;
    int ndm;       // always equal to 2^depth
    int nsquish;   // number of time squishings

    float *channel_freqs;       // length p.nchan
    float *dm_freqs;            // length ndm
    double output_normalization;

    double tree_dm_min;
    double tree_dm_max;
    double dt_squished;
    int ndata_squished;       // always equal to nsamples_per_chunk / 2^nsquish

    float **ring_buffer;
    int ring_t0;  // keeps track of position in ring buffer
    int nring;

    frb_simple_tree_search_algorithm(const frb_search_params &p, int depth, int nsquish);
    virtual ~frb_simple_tree_search_algorithm();

    // Helper functions for search_chunk
    float **dedisperse_chunk_incremental(float **buf);
    float **dedisperse_chunk_non_incremental(float **buf);
    void postprocess_chunk_incremental(int ichunk, int idm, int it0, int it1, float *debug_buffer);
    void postprocess_chunk_non_incremental(float **buf, int idm, int it0, int it1, float *debug_buffer);

    virtual void  search_start();
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();

    static frb_search_algorithm_base::ptr_t create(const vector<string> &tokens, const frb_search_params &p);
};


frb_simple_tree_search_algorithm::frb_simple_tree_search_algorithm(const frb_search_params &p_, int depth_, int nsquish_)
    : frb_search_algorithm_base(p_), depth(depth_), ndm(1 << depth_), nsquish(nsquish_)
{ 
    xassert(depth >= 1);
    xassert(nsquish >= 0);
    xassert(nsquish <= 6);
    xassert(p.nsamples_per_chunk % (1 << nsquish) == 0);

    stringstream s;
    s << "simple_tree-" << ndm << "-" << nsquish;
    this->name = s.str();

    this->channel_freqs = jstree::get_chime_channels(p.nchan, p.band_freq_lo_MHz, p.band_freq_hi_MHz);
    this->dm_freqs = jstree::get_dm_channels(this->depth, channel_freqs[0], channel_freqs[p.nchan-1]);
    this->output_normalization = 1.0 / sqrt(p.nchan * (1<<nsquish));

    this->dt_squished = p.dt_sample * (1 << nsquish);
    this->ndata_squished = p.nsamples_per_chunk / (1 << nsquish);

    // No "lagging" in this version!
    this->tree_dm_min = 0.0;
    this->tree_dm_max = jstree::get_diagonal_dm(dt_squished, dm_freqs);

    this->ring_buffer = NULL;
    this->ring_t0 = 0;
    this->nring = 0;

    this->debug_buffer_ndm = ndm;
    this->debug_buffer_nt = p.nsamples_tot / (1 << nsquish);

    this->search_gb = 0.0;
    this->search_gb += 1.0e-9 * p.nchan * ndata_squished * sizeof(float);
    this->search_gb += 2.0e-9 * ndm * ndata_squished * sizeof(float);

    if (p.nchunks > 1)
	this->search_gb += 2.0e-9 * ndm * ndm * sizeof(float);
	
    cout << this->name << " initialized, tree_dm_min=" << tree_dm_min << ", tree_dm_max=" << tree_dm_max << endl;

    if ((p.dm_min < tree_dm_min-0.1) || (p.dm_max > tree_dm_max+0.1))
	throw runtime_error("Fatal: DM range in search_params exceeds range searched by tree\n");
}


frb_simple_tree_search_algorithm::~frb_simple_tree_search_algorithm()
{
    free(channel_freqs);
    free(dm_freqs);

    channel_freqs = NULL;
    dm_freqs = NULL;
}


void frb_simple_tree_search_algorithm::search_start()
{
    this->search_result = -1.0e30;

    if (p.nchunks > 0) {
	this->ring_t0 = 0;
	this->nring = ndata_squished + ndm;
	this->ring_buffer = jstree::matrix(ndm, nring);
    }
}


float **frb_simple_tree_search_algorithm::dedisperse_chunk_incremental(float **buf)
{
    int ndata2 = ndata_squished + ndm;

    // Copy to overallocated buffer
    float **buf2 = jstree::matrix(ndm, ndata2);
    for (int i = 0; i < ndm; i++)
	memcpy(buf2[i], buf[i], ndata_squished * sizeof(float));

    // Reallocate original buffer at larger size
    free(buf[0]);
    free(buf);
    buf = jstree::matrix(ndm, ndata2);

    jstree::dedisperse_lagged(buf2, buf, ndm, ndata_squished);
    jstree::update_ring_buffer(buf, ring_buffer, ndm, ndata_squished, nring, &ring_t0);

    free(buf2[0]);
    free(buf2);

    return buf;
}


float **frb_simple_tree_search_algorithm::dedisperse_chunk_non_incremental(float **buf)
{
    xassert(p.nchunks == 1);

    float **buf2 = jstree::matrix(ndm, ndata_squished);
    jstree::dedisperse(buf, buf2, ndm, ndata_squished);
    free(buf[0]);
    free(buf);

    return buf2;
}


void frb_simple_tree_search_algorithm::postprocess_chunk_incremental(int ichunk, int idm, int it0, int it1, float *debug_buffer)
{
    xassert(idm >= 0 && idm < ndm);
    xassert(nring >= ndata_squished + ndm);

    // Formal start time of ring buffer
    int it_ring = (ichunk+1)*ndata_squished - nring;

    // Clamp to "valid" ring buffer range
    it0 = max(it0, 0);
    it0 = max(it0, it_ring + nring - idm - ndata_squished);
    it1 = min(it1, it_ring + nring - idm);
    
    //
    // Slow way to implement the update (better to avoid %) but
    // this isn't really supposed to be an optimized implementation
    //
    for (int it = it0; it < it1; it++) {
	int j = (ring_t0 + it - it_ring) % nring;   // first arugment is always positive
	this->search_result = max(this->search_result, ring_buffer[idm][j]);
    }

    if (debug_buffer != NULL) {
	for (int it = it0; it < it1; it++) {
	    int j = (ring_t0 + it - it_ring) % nring;
	    debug_buffer[idm*debug_buffer_nt + it] = ring_buffer[idm][j];
	}
    }
	
}


void frb_simple_tree_search_algorithm::postprocess_chunk_non_incremental(float **buf, int idm, int it0, int it1, float *debug_buffer)
{
    xassert(p.nchunks == 1);
    xassert(idm >= 0 && idm < ndm);

    // clamp
    it0 = max(it0, 0);
    it1 = min(it1, ndata_squished - idm);

    for (int it = it0; it < it1; it++)
	this->search_result = max(this->search_result, buf[idm][it]);
	
    if (debug_buffer != NULL) {
	for (int it = it0; it < it1; it++)
	    debug_buffer[idm*debug_buffer_nt + it] = buf[idm][it];
    }	    
}


void frb_simple_tree_search_algorithm::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    float **data = jstree::matrix(p.nchan, p.nsamples_per_chunk);
    for (int i = 0; i < p.nchan; i++)
	for (int j = 0; j < p.nsamples_per_chunk; j++)
	    data[i][j] = output_normalization * chunk[i*p.nsamples_per_chunk+j];   // Note output_normalization here

    for (int isquish = 0; isquish < nsquish; isquish++) {
	int ndata = p.nsamples_per_chunk / (1 << isquish);
	float tsamp = p.dt_sample * (1 << isquish);

	float **tmp = data;
	data = jstree::time_squish_data(data, channel_freqs, p.nchan, ndata, tsamp, 0.0);   // no lagging in this version!
	free(tmp[0]); 
	free(tmp);
    }

    float unused = 0.0;
    float **buf = jstree::remap_data(ndata_squished, p.nchan, channel_freqs, depth, dt_squished, &unused, data);
    free(data[0]); 
    free(data);

    // At this point, buf is an array with shape (ndm, ndata_squished)

    if (p.nchunks > 1)
	buf = this->dedisperse_chunk_incremental(buf);
    else
	buf = this->dedisperse_chunk_non_incremental(buf);

    // This is the DM search range, shifted and normalized to output array indices
    float dm0 = (p.dm_min - tree_dm_min) / (tree_dm_max - tree_dm_min) * (ndm-1);
    float dm1 = (p.dm_max - tree_dm_min) / (tree_dm_max - tree_dm_min) * (ndm-1);
    
    // This is the integer range of DM indices to be searched
    int idm0 = max(0, int(dm0));         // round down
    int idm1 = min(ndm, int(dm1+1.0));   // round up

    for (int idm = idm0; idm < idm1; idm++) {
	float dm = ((ndm-1-idm)*tree_dm_min + idm*tree_dm_max) / (float)(ndm-1);
	
	//
	// Using the same logic as the pulse simulator itself, determine the allowed
	// range of arrival times for this DM.  (The pulse simulator restricts the
	// range so that the entire pulse lies within the simulated timestream.)
	//
	double t0, t1;
	p.get_allowed_arrival_times(t0, t1, 0.0, dm, 0.0);
	
	//
	// Detail: the times t0,t1 are undispersed arrival times, whereas the
	// tree dedispersion code is indexing based on arrival time in the lowest
	// DM channel, so calculate the difference between the two.
	//	
	double dt = dispersion_delay(dm, dm_freqs[0]);
	
	// Index range for time search (relative to start of timestream, unclamped)
 	int it0 = int((t0+dt)/dt_squished);
	int it1 = int((t1+dt)/dt_squished+2.0);

	if (p.nchunks > 1)
	    this->postprocess_chunk_incremental(ichunk, idm, it0, it1, debug_buffer);
	else
	    this->postprocess_chunk_non_incremental(buf, idm, it0, it1, debug_buffer);
    }

    free(buf[0]); 
    free(buf);
}

void frb_simple_tree_search_algorithm::search_end()
{
    xassert(this->search_result > -1.0e30);

    this->ring_t0 = 0;
    this->nring = 0;

    if (this->ring_buffer) {
	free(ring_buffer[0]);
	free(ring_buffer);
	ring_buffer = NULL;
    }
}


// Registry boilerplate follows

// static member function
frb_search_algorithm_base::ptr_t frb_simple_tree_search_algorithm::create(const vector<string> &tokens, const frb_search_params &p)
{
    if (tokens.size() != 2)
	return frb_search_algorithm_base::ptr_t();

    int ndm = xlexical_cast<int> (tokens[0], "simple_tree ndm");
    int nsquish = xlexical_cast<int> (tokens[1], "simple_tree nsquish");

    xassert(is_power_of_two(ndm));
    int depth = integer_log2(ndm);
    return boost::make_shared<frb_simple_tree_search_algorithm>(p, depth, nsquish);
}

static const char *usage = (
"    simple_tree <ndm> <nsquish>\n"
"       note that the downsampling factor is 2^nsquish\n"
);

// register tree algorithm at library load time
namespace {
    struct _initializer {
	_initializer()
	{
	    frb_search_algorithm_base::register_algorithm("simple_tree", frb_simple_tree_search_algorithm::create, usage);
	}
    } _init;
}

}  // namespace frb_olympics

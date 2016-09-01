#include <memory>
#include <cstring>
#include "frb_olympics.hpp"

#include "bonsai.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


//
// Each frb_olympics search algorithm being compared has a 'name' string for plot labels etc.
//
// Since I often include the bonsai algorithm several times with different config files,
// I wrote this function to produce a name string from the filename, with some cosmetic
// edits like replacing underscores with hyphens (the underscores confuse matplotlib if
// texstrings=true).
//
static string name_from_filename(const string &filename)
{
    char *s = strdup(filename.c_str());

    int n = strlen(s);

    if ((n > 5) && !strcmp(s+n-5, ".hdf5")) {
	s[n-5] = 0;
	n -= 5;
    }
    else if ((n > 3) && !strcmp(s+n-3, ".h5")) {
	s[n-3] = 0;
	n -= 3;
    }

    for (int i = 0; i < n; i++) {
	if (s[i] == '_')
	    s[i] = '-';
    }

    string ret = s;
    free(s);

    return ret;
}


//
// A simple subclass of bonsai::dedispserser
//
//   - trigger processing just keeps track of the global maximum.
//
//   - if frb_olympics is being run in 'debug' mode, meaning that the triggers are being
//     written to an output array for visual inspection or later postprocessing, we simply
//     copy coarse-grained triggers to the debug array.
//
struct frb_olympics_dedisperser : public bonsai::dedisperser {    
    float global_max_trigger;

    // If the frb_olympics run uses a debug_buffer, then we stash the bare pointer
    // here during the call to dedisperser::run(), so that it can be accessed during
    // dedisperser::process_triggers().
    float *debug_buffer;
    
    // These fields are taken from the frb_search_params.
    int debug_buffer_ndm;
    int debug_buffer_nt;

    // Position in debug buffer (starts at zero and increases)
    int debug_buffer_tpos;

    frb_olympics_dedisperser(const bonsai::config_params &cp, int ibeam)
	: bonsai::dedisperser(cp, ibeam, true),     // (params, ibeam, init_weights)
	  global_max_trigger(-1.0e30), 
	  debug_buffer(nullptr), debug_buffer_ndm(0), 
	  debug_buffer_nt(0), debug_buffer_tpos(0)
    { }

    virtual void preprocess_data(float *data, float *weights, int data_stride) { return; }
    virtual void finalize() { return; }

    virtual void process_triggers(int itree, const bonsai::trigger_set &ts)
    {
	for (int i = 0; i < ts.ntm_tot; i++)
	    global_max_trigger = max(global_max_trigger, ts.tm[i]);

	if (!this->debug_buffer || (itree != 0))
	    return;

	int ts_ndm = ts.ntree_tot / ts.ndm_per_trigger;
	int ts_nt = ts.ntime_triggers_per_transform;

	xassert(debug_buffer_ndm == ts_ndm);
	xassert(debug_buffer_tpos >= 0);
	xassert(debug_buffer_tpos + ts_nt <= debug_buffer_nt);
	
	// This (ism, ibeta) pair will be used to fill the debug_buffer
	int ism = 0;
	int ibeta = ts.nbeta/2;

	// Overall offset and DM stride in trigger_set array
	int s0 = (ism*ts.nbeta + ibeta) * ts_nt;
	int dstride = ts.nsm * ts.nbeta * ts_nt;

	for (int idm = 0; idm < ts_ndm; idm++)
	    for (int it = 0; it < ts_nt; it++)
		this->debug_buffer[idm*debug_buffer_nt + debug_buffer_tpos + it] = ts.tm[idm*dstride + s0 + it];
	
	this->debug_buffer_tpos += ts_nt;
    }
};


// -------------------------------------------------------------------------------------------------


struct frb_bonsai_search_algorithm : frb_search_algorithm_base
{
    // Initialized at construction
    shared_ptr<bonsai::config_params> config;

    //
    // We allocate a new dedisperser object in search_start() whenever a new simulation is
    // analyzed, and destroy it in search_end().  This is so that all memory buffers are freed
    // before going on to the next algorithm.  Eventually I'll implement dynamic buffer freeing/
    // reallocation in libbonsai.  Then this logic can be improved by creating the dedisperser
    // object once, and allocating/freeing its buffers for each simulation.
    //
    shared_ptr<frb_olympics_dedisperser> dp;

    frb_bonsai_search_algorithm(const string &hdf5_filename);
    virtual ~frb_bonsai_search_algorithm() { }

    // Devirtualize frb_search_algorithm_base
    virtual void  search_init(const frb_search_params &p);
    virtual void  search_start(int mpi_rank_within_node);
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();
};


// Constructor just reads the config file
frb_bonsai_search_algorithm::frb_bonsai_search_algorithm(const string &hdf5_filename)
{
    bonsai::hdf5_file f(hdf5_filename);
    bonsai::hdf5_group g(f, ".");

    this->name = name_from_filename(hdf5_filename);
    this->config = make_shared<bonsai::config_params> (g, true);
}


//
// Called when the frb_olympics search_params are determined.
//
// We just check consistency with the bonsai config file, and set bonsai::config_params::nt_data
// (the number of time samples which bonsai expects to receive in each input chunk) to the chunk
// size in the frb_olympics run.
//
void frb_bonsai_search_algorithm::search_init(const frb_search_params &p)
{
    xassert(config->nchan == p.nchan);
    xassert(fabs(config->freq_lo_MHz - p.band_freq_lo_MHz) < 0.01);
    xassert(fabs(config->freq_hi_MHz - p.band_freq_hi_MHz) < 0.01);
    xassert(fabs(config->dt_sample - p.dt_sample) < 1.0e-5);

    for (int itree = 0; itree < config->ntrees; itree++) {
	int nt_tree = (config->nt_tree[itree] * config->nds[itree]) / config->nups[itree];
	xassert(p.nsamples_tot % nt_tree == 0);
    }

    config->set_nt_data(p.nsamples_per_chunk);

    // Tree 0 is used to set debug_buffer parameters
    int ndm = config->tree_size[0];
    int nt = (p.nsamples_tot / config->nds[0]) * config->nups[0];
    int ndm_per_trigger = config->ndm_per_trigger[0];
    int nt_per_trigger = config->nt_per_trigger[0];

    xassert(ndm % ndm_per_trigger == 0);
    xassert(nt % nt_per_trigger == 0);

    this->search_params = p;
    this->debug_buffer_ndm = ndm / ndm_per_trigger;
    this->debug_buffer_nt = nt / nt_per_trigger;

    // FIXME loose end here: enable memory profiling
    this->search_gb = 0.1;
}


//
// Called whenever a new simulation is processed.
// As described above, we allocate the bonsai::dedisperser object here.
//
void frb_bonsai_search_algorithm::search_start(int mpi_rank_within_node)
{
    if (this->dp)
	throw runtime_error("frb_bonsai_search_algorithm::search_start(): dedisperser pointer is set?!");

    this->dp = make_shared<frb_olympics_dedisperser> (*config, mpi_rank_within_node);
    this->search_result = -1.0e30;

    dp->debug_buffer_ndm = this->debug_buffer_ndm;
    dp->debug_buffer_nt = this->debug_buffer_nt;

    // Must always be called after calling dedisperser constructor
    dp->spawn_slave_threads();
}


//
// Called for each chunk of data in the frb_olympics simulation.
//
void frb_bonsai_search_algorithm::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    // libbonsai channel ordering is reversed relative to frb_olympics, so use negative stride
    const float *data = &chunk[(search_params.nchan-1) * search_params.nsamples_per_chunk];
    int data_stride = -search_params.nsamples_per_chunk;

    // FIXME this is a kludge.
    vector<float> weight_vec(search_params.nchan * search_params.nsamples_per_chunk, 1.0);
    float *weights = &weight_vec[(search_params.nchan-1) * search_params.nsamples_per_chunk];
    
    // The frb_olympics_dedisperser doesn't actually modify the chunk in its preprocess_data() virtual, so const_cast is OK here.
    dp->debug_buffer = debug_buffer;
    dp->run(const_cast<float *> (data), weights, data_stride);
    dp->debug_buffer = nullptr;

    // FIXME it would make more sense to put this in search_end(), but this currently doesn't work.
    this->search_result = dp->global_max_trigger;
}


//
// Called at the end of each simulation.
// As described above, we deallocate the bonsai::dedisperser object here.
//
void frb_bonsai_search_algorithm::search_end()
{
    if (!this->dp)
	throw runtime_error("frb_bonsai_search_algorithm::search_end(): dedisperser pointer is not set?!");

    // It's important to call dedisperser::terminate() before destroying the dedisperser object
    this->dp->terminate();

    this->dp = shared_ptr<frb_olympics_dedisperser> ();
}


frb_search_algorithm_base *bonsai(const string &hdf5_filename)
{
    return new frb_bonsai_search_algorithm(hdf5_filename);
}


}  // namespace frb_olympics

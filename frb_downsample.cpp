#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


struct frb_downsample_search : public frb_search_algorithm_base
{
    frb_search_algorithm_base::ptr_t child;
    int nds;

    frb_downsample_search(const frb_search_params &p, frb_search_algorithm_base::ptr_t child, int nds);

    virtual ~frb_downsample_search() { }

    virtual void  search_start();
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();

    static frb_search_algorithm_base::ptr_t create(const vector<string> &tokens, const frb_search_params &p);
};

    
frb_downsample_search::frb_downsample_search(const frb_search_params &p_, frb_search_algorithm_base::ptr_t child_, int nds_)
    : frb_search_algorithm_base(p_), child(child_), nds(nds_)
{ 
    // more asserts could be added...
    xassert(nds >= 1);
    xassert(p.nsamples_per_chunk == child->p.nsamples_per_chunk * nds);
    xassert(p.nsamples_tot == child->p.nsamples_tot * nds);
    xassert(p.nchunks == child->p.nchunks);
    xassert(fabs(p.dt_sample - child->p.dt_sample/nds) < 1.0e-10*p.dt_sample);

    stringstream s;
    s << "downsample-" << nds << "-" << child->name;
    this->name = s.str();

    this->debug_buffer_ndm = child->debug_buffer_ndm;
    this->debug_buffer_nt = child->debug_buffer_nt;
    this->search_gb = child->search_gb;
    this->search_gb += 1.0e-9 * p.nchan * child->p.nsamples_per_chunk * sizeof(float);
}


void frb_downsample_search::search_start()
{
    child->search_start();
    this->search_result = child->search_result;
}


void frb_downsample_search::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    int nt = child->p.nsamples_per_chunk;  // number of downsampled samples in chunk
    vector<float> buf(p.nchan * nt, 0.0);

    double w = 1.0 / sqrt(nds);

    for (int ichan = 0; ichan < p.nchan; ichan++)
	for (int it = 0; it < nt; it++)
	    for (int j = it*nds; j < (it+1)*nds; j++)
		buf[ichan*nt + it] += w * chunk[ichan*nt*nds + j];

    child->search_chunk(&buf[0], ichunk, debug_buffer);
    this->search_result = child->search_result;
}

void frb_downsample_search::search_end()
{
    child->search_end();
}


// Registry boilerplate follows

// static member function
frb_search_algorithm_base::ptr_t frb_downsample_search::create(const vector<string> &tokens, const frb_search_params &p)
{
    if (tokens.size() < 2)
	return frb_search_algorithm_base::ptr_t();

    int nds = xlexical_cast<int> (tokens[0], "downsample nds");
    
    xassert(nds >= 1);
    xassert(p.nsamples_per_chunk % nds == 0);

    frb_search_params child_search_params = p;
    child_search_params.nsamples_per_chunk = p.nsamples_per_chunk / nds;
    child_search_params.nsamples_tot = p.nsamples_tot / nds;
    child_search_params.dt_sample = p.dt_sample * nds;

    vector<string> tokens2 = tokens;
    tokens2.erase(tokens2.begin(), tokens2.begin()+1);
    frb_search_algorithm_base::ptr_t child = frb_search_algorithm_base::parse_line(tokens2, child_search_params);

    return boost::make_shared<frb_downsample_search>(p, child, nds);
}
 
static const char *usage = (
"    downsample <downsample_factor> ...\n"
"        Here, '...' stands for any other algorithm (including all command-line parameters)\n"
);

namespace {
    struct _initializer {
	_initializer()
	{
	    frb_search_algorithm_base::register_algorithm("downsample", frb_downsample_search::create, usage);
	}
    } _init;
}

}  // namespace frb_olympics

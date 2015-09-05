#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


struct frb_rechunked_search : public frb_search_algorithm_base
{
    frb_search_algorithm_base::ptr_t child;
    vector<float> buf;

    frb_rechunked_search(const frb_search_params &p, frb_search_algorithm_base::ptr_t child);

    virtual ~frb_rechunked_search() { }

    virtual void  search_start();
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer);
    virtual void  search_end();

    static frb_search_algorithm_base::ptr_t create(const vector<string> &tokens, const frb_search_params &p);
};


    
frb_rechunked_search::frb_rechunked_search(const frb_search_params &p_, frb_search_algorithm_base::ptr_t child_)
    : frb_search_algorithm_base(p_), child(child_)
{ 
    // more asserts could be added...
    xassert(child->p.nchan == p.nchan);
    xassert(child->p.nsamples_tot == p.nsamples_tot);
    xassert(child->p.nsamples_per_chunk % p.nsamples_per_chunk == 0);
    xassert(child->p.nsamples_tot == child->p.nchunks * child->p.nsamples_per_chunk);

    stringstream s;
    s << "rechunked-" << child->p.nsamples_per_chunk << "-" << child->name;
    this->name = s.str();

    this->debug_buffer_ndm = child->debug_buffer_ndm;
    this->debug_buffer_nt = child->debug_buffer_nt;
    this->search_gb = child->search_gb;
    this->search_gb += 1.0e-9 * p.nchan * child->p.nsamples_per_chunk * sizeof(float);
}


void frb_rechunked_search::search_start()
{
    child->search_start();
    buf.resize(p.nchan * child->p.nsamples_per_chunk);
    this->search_result = child->search_result;
}


void frb_rechunked_search::search_chunk(const float *chunk, int ichunk, float *debug_buffer)
{
    int n = p.nsamples_per_chunk;
    int m = child->p.nsamples_per_chunk / p.nsamples_per_chunk;

    for (int ichan = 0; ichan < p.nchan; ichan++)
	memcpy(&buf[ichan*m*n + (ichunk%m)*n], &chunk[ichan*n], n * sizeof(buf[0]));

    if ((ichunk % m) == (m-1)) {
	child->search_chunk(&buf[0], ichunk/m, debug_buffer);
	this->search_result = child->search_result;
    }
}

void frb_rechunked_search::search_end()
{
    deallocate(buf);
    child->search_end();
}


// Registry boilerplate follows

// static member function
frb_search_algorithm_base::ptr_t frb_rechunked_search::create(const vector<string> &tokens, const frb_search_params &p)
{
    if (tokens.size() < 2)
	return frb_search_algorithm_base::ptr_t();

    int child_nsamples_per_chunk = xlexical_cast<int> (tokens[0], "frb rechunk factor");
    xassert(child_nsamples_per_chunk % p.nsamples_per_chunk == 0);
    xassert(p.nsamples_tot % child_nsamples_per_chunk == 0);

    frb_search_params child_search_params = p;
    child_search_params.nsamples_per_chunk = child_nsamples_per_chunk;
    child_search_params.nchunks = p.nsamples_tot / child_nsamples_per_chunk;

    vector<string> tokens2 = tokens;
    tokens2.erase(tokens2.begin(), tokens2.begin()+1);
    frb_search_algorithm_base::ptr_t child = frb_search_algorithm_base::parse_line(tokens2, child_search_params);

    return boost::make_shared<frb_rechunked_search>(p, child);
}
 
static const char *usage = (
"    rechunk <rechunk_factor> ...\n"
"        Here, '...' stands for any other algorithm (including all command-line parameters)\n"
);

namespace {
    struct _initializer {
	_initializer()
	{
	    frb_search_algorithm_base::register_algorithm("rechunk", frb_rechunked_search::create, usage);
	}
    } _init;
}

}  // namespace frb_olympics

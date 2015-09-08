#ifndef _FRB_OLYMPICS_HPP
#define _FRB_OLYMPICS_HPP

#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

#include <vector>
#include <iostream>
#include <exception>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/noncopyable.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>

// osx
#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
#endif

#ifndef xassert
#define xassert(cond) xassert2(cond, __LINE__)
#define xassert2(cond,line) do { if (!(cond)) throw std::runtime_error("Assertion '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")"); } while (0)
#endif

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


//
// Note: throughout the code, all times are in seconds, and all frequencies 
// are in MHz, with the exception of the scattering measure (SM).  We define
// the SM to be the scatter-broadening in MILLISECONDS at 1 GHz.
//

inline double dispersion_delay(double dm, double freq_MHz)
{
    return 4.148806e3 * dm / (freq_MHz * freq_MHz);
}

inline double scatter_broadening(double sm, double freq_MHz)
{
    return 1.0e9 * sm / (freq_MHz * freq_MHz * freq_MHz * freq_MHz);
}


// boilerplate for error message generation
template<typename T> inline const char *tname();
template<> inline const char *tname<int>() { return "integer"; }
template<> inline const char *tname<double>() { return "double"; }

// C++ doesn't provide a one-liner to deallocate a vector!
template<typename T> inline void deallocate(std::vector<T> &v)
{
    std::vector<T> w;
    v.swap(w);
}

// version of boost::lexical_cast which generates a more useful error message
template<typename T> inline T xlexical_cast(const std::string &s, const std::string &prefix)
{
    try {
	return boost::lexical_cast<T>(s);
    } catch (...) {
	std::stringstream s;
	s << prefix << ": couldn't parse token " << s << " (expected type: " << tname<T>() << ")\n";
	throw std::runtime_error(s.str().c_str());
    }
}

template<typename T> static inline T square(const T &x) 
{
    return (x*x);
}

// works for positive or negative x
static inline int round_down(double x)
{
    int i = (int)x;
    return (i <= x) ? i : (i-1);
}

static inline int round_up(double x)
{
    int i = (int)x;
    return (i >= x) ? i : (i+1);
}

// frb_misc.cpp
extern void get_pulse_endpoints(double &t0, double &t1, double intrinsic_width, double dm, double sm, double freq_lo_MHz, double freq_hi_MHz);
extern void tokenize_file(const std::string &filename, std::vector<std::vector<std::string> > &tokens);
extern void read_kv_pairs(const std::string &filename, boost::unordered_map<std::string,std::string> &kv_pairs);
extern bool is_power_of_two(int n);
extern int  round_up_to_power_of_two(int n);
extern int  integer_log2(int n);     // aborts if n is not a power of two
extern double time_diff(const struct timeval &tv1, const struct timeval &tv2);
extern struct timeval get_time();


// -------------------------------------------------------------------------------------------------


class frb_rng {
protected:
    boost::random::mt19937 generator;
    boost::random::uniform_01<> uniform_dist;
    boost::random::normal_distribution<> gaussian_dist;

public:
    frb_rng();

    double uniform(double lo, double hi);
    double gaussian();   // zero mean, unit variance
};


// Lightweight structure containing parameters for a single FRB event
struct frb_pulse {
    //
    // Note: the fluence and spectral index use whatever normalization
    // the timesterams use.  (E.g. if the timestreams are in units of
    // sky temperature, this includes a 1/nu^2 weighting compared to
    // the usual radio astronomy normalization in Janskies, so the spectral
    // index is shifted by 2.)
    //
    double   fluence;   // at 1 GHz
    double   spectral_index;

    //
    //
    // Note: scattering is implemented by applying a Gaussian convolution
    // with width proportional to nu^(-4), with the central arrival time
    // of the pulse unchanged (i.e. there is no frequency-dependent
    // scattering delay).  Are these the correct physical assumptions?
    //
    // We define the "scattering measure" SM to be the scatter-broadening
    // in MILLISECONDS at 1 GHz (this is the only exception to the statement
    // that all times are in seconds and all frequencies in MHz).
    //
    double   arrival_time;          // seconds
    double   intrinsic_width;       // seconds
    double   dispersion_measure;    // pc cm^{-3}
    double   scattering_measure;    // scatter broadening in msec at 1 GHz

    frb_pulse();  // default constructor catering to cython

    frb_pulse(double fluence, double arrival_time, double intrinsic_width,
	      double dispersion_measure, double scattering_measure, double spectral_index);

    void    get_endpoints(double &t0, double &t1, double freq_lo_MHz, double freq_hi_MHz) const;
    void    add_to_timestream(double freq_lo_MHz, double freq_hi_MHz, float *timestream, int nsamples_per_chunk, double dt_sample, int ichunk) const;
    double  get_signal_to_noise_in_channel(double freq_lo_MHz, double freq_hi_MHz, double dt_sample) const;
};


struct frb_search_params {
    // Search space
    double  dm_min, dm_max;              // dispersion measure
    double  sm_min, sm_max;              // scattering measure
    double  beta_min, beta_max;          // spectral index
    double  width_min, width_max;        // intrinsic width

    // Frequency channels
    // For now, only equal spacing in frequency is supported
    int     nchan;
    double  band_freq_lo_MHz;
    double  band_freq_hi_MHz;

    // Time sampling
    double  dt_sample;     // seconds
    int     nsamples_tot;  // always equal to (nsamples_per_chunk * nchunks)
    int     nsamples_per_chunk;
    int     nchunks;

    // Default constructor needed temporarily in frb_search_algorithm_base
    frb_search_params() { memset(this, 0, sizeof(*this)); }

    // Construct from file (see frb_search_params.cpp for a description of file format)
    explicit frb_search_params(const std::string &filename);
    
    // Arrival time will be randomly chosen so that entire pulse lies in timestream
    frb_pulse make_random_pulse(frb_rng &r, double fluence) const;

    // @timestream should point to a buffer of size nchan * nsamples_per_chunk
    void    get_allowed_arrival_times(double &t0, double &t1, double intrinsic_width, double dm, double sm) const;
    void    simulate_noise(frb_rng &r, float *timestream) const;
    void    add_pulse(const frb_pulse &pulse, float *timestream, int ichunk) const;
    double  get_signal_to_noise_of_pulse(const frb_pulse &pulse) const;

    void write(std::ostream &os) const;
    void write() const { this->write(std::cout); }

    inline double freq_lo_of_channel(int i) const
    { 
	return ((nchan-i)*band_freq_lo_MHz + (i)*band_freq_hi_MHz) / (double)nchan;
    }

    inline double freq_hi_of_channel(int i) const
    { 
	return ((nchan-i-1)*band_freq_lo_MHz + (i+1)*band_freq_hi_MHz) / (double)nchan;
    }

    // A useful helper function for direct searches
    void make_dm_table(std::vector<double> &dm_table, double epsilon) const;
};


struct frb_search_algorithm_base : boost::noncopyable
{
    // initialized in constructor
    std::string name;

    // initialized in search_init()
    frb_search_params search_params;
    int debug_buffer_ndm;
    int debug_buffer_nt;
    double search_gb;

    // initialized in search_start()
    float search_result;

    frb_search_algorithm_base()
	: name(), search_params(), debug_buffer_ndm(0), debug_buffer_nt(0), search_gb(0.0), search_result(0.0)
    { }

    virtual ~frb_search_algorithm_base() { }

    //
    // Each FRB search algorithm should subclass frb_search_algorithm_base, and define the following virtuals.
    //
    // The basic call sequence for an FRB search is
    //    - one call to search_start(), this will probably just allocate buffers
    //    - N calls to search_chunk(), where N is the number of "chunks" that the timestream is divided into
    //    - one call to search_end(), which should deallocate all buffers
    //
    // The virtual function search_estimate_gb() should return an estimate of the total buffer size in GB.
    //
    virtual void  search_init(const frb_search_params &p) = 0;
    virtual void  search_start() = 0;
    virtual void  search_chunk(const float *chunk, int ichunk, float *debug_buffer) = 0;
    virtual void  search_end() = 0;
};


// Algorithms
extern frb_search_algorithm_base *simple_direct(double epsilon);
extern frb_search_algorithm_base *simple_tree(int depth, int nsquish);
extern frb_search_algorithm_base *sloth(double epsilon, int nupsample);
extern frb_search_algorithm_base *bonsai(int depth, int nupsample);


}  // namespace frb_olympics


#endif  // _FRB_OLYMPICS_HPP

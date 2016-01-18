from libcpp.string cimport string
from libcpp cimport bool


cdef extern from "frb_olympics.hpp" namespace "frb_olympics":
    cdef cppclass frb_rng:
        frb_rng() except +
        frb_rng(frb_rng &) except +

        double uniform(double lo, double hi) except +
        double gaussian() except +


    cdef cppclass frb_pulse:
        frb_pulse()

        frb_pulse(double fluence, double arrival_time, double intrinsic_width,
                  double dispersion_measure, double scattering_measure, double spectral_index) except +
        
        double fluence
        double spectral_index
        double arrival_time
        double intrinsic_width
        double dispersion_measure
        double scattering_measure

        double get_signal_to_noise_in_channel(double freq_lo_MHz, double freq_hi_MHz, double dt_sample) except +

        void get_endpoints(double &t0, double &t1, double freq_lo_MHz, double freq_hi_MHz) except +
        void add_to_timestream(double freq_lo_MHz, double freq_hi_MHz, float *timestream, int nsamples_per_chunk, double dt_sample, int ichunk) except +


    cdef cppclass frb_search_params:
        frb_search_params(const frb_search_params &) except +
        frb_search_params(string &filename) except +

        double  dm_min, dm_max
        double  sm_min, sm_max
        double  beta_min, beta_max
        double  width_min, width_max

        int     nchan
        double  band_freq_lo_MHz
        double  band_freq_hi_MHz
        
        double  dt_sample
        int     nsamples_tot
        int     nsamples_per_chunk
        int     nchunks

        frb_pulse  make_random_pulse(frb_rng &r, double fluence) except +
        void       get_allowed_arrival_times(double &t0, double &t1, double intrinsic_width, double dm, double sm) except +
        void       simulate_noise(frb_rng &r, float *timestream) except +
        void       add_pulse(frb_pulse &pulse, float *timestream, int ichunk) except +
        double     get_signal_to_noise_of_pulse(frb_pulse &pulse) except +
        void       write() except +


    cdef cppclass frb_search_algorithm_base:
        frb_search_params  search_params
        string             name
        int                debug_buffer_ndm
        int                debug_buffer_nt
        float              search_gb
        float              search_result

        void   search_init(const frb_search_params &p) except +
        void   search_start(int mpi_rank_within_node) except +
        void   search_chunk(const float *chunk, int ichunk, float *debug_buffer) except +
        void   search_end() except +

    cdef frb_search_algorithm_base *simple_direct(double epsilon) except +
    cdef frb_search_algorithm_base *sloth(double epsilon_d, int nsm, int nbeta, int nupsample, bool strict_incremental) except +
    cdef frb_search_algorithm_base *bonsai(const string &hdf5_filename) except +

    cdef frb_search_algorithm_base *sloth_sm_subsearch(double epsilon_d, double sm, int nbeta, int nupsample, bool strict_incremental) except +

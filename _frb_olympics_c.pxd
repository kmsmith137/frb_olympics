cdef extern from "frb_olympics.hpp" namespace "frb_olympics":

    cdef cppclass frb_pulse:
        frb_pulse(double fluence, double arrival_time, double intrinsic_width,
                  double dispersion_measure, double scattering_measure, double spectral_index) except +
        
        double fluence
        double spectral_index
        double arrival_time
        double intrinsic_width
        double dispersion_measure
        double scattering_measure

        double get_signal_to_noise_in_channel(double freq_lo_MHz, double freq_hi_MHz, double dt_sample) except +
        
        void add_to_timestream(double freq_lo_MHz, double freq_hi_MHz, float *timestream, int nsamples_per_chunk, double dt_sample, int ichunk) except +

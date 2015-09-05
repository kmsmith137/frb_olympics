import numpy as np
cimport numpy as np
cimport _frb_olympics_c


cdef class frb_pulse:
    cdef _frb_olympics_c.frb_pulse *_pulse

    def __cinit__(self):
        self._pulse = NULL

    def __init__(self, fluence, arrival_time, intrinsic_width, dispersion_measure, scattering_measure, spectral_index):
        self._pulse = new _frb_olympics_c.frb_pulse(fluence, arrival_time, intrinsic_width, dispersion_measure, scattering_measure, spectral_index)
     	 
    def __dealloc__(self):
        if self._pulse != NULL:
            del self._pulse
            self._pulse = NULL


    def get_signal_to_noise_in_channel(self, freq_lo_MHz, freq_hi_MHz, dt_sample):
        assert self._pulse != NULL
        return self._pulse.get_signal_to_noise_in_channel(freq_lo_MHz, freq_hi_MHz, dt_sample)
     
    def add_to_timestream(self, freq_lo_MHz, freq_hi_MHz, np.ndarray[float,ndim=1,mode='c'] timestream not None, dt_sample, ichunk):
        assert self._pulse != NULL
        self._pulse.add_to_timestream(freq_lo_MHz, freq_hi_MHz, &timestream[0], timestream.shape[0], dt_sample, ichunk)

    def __repr__(self):
        return ('frb_pulse(fluence=%s, arrival_time=%s, width=%s, dm=%s, sm=%s, beta=%s)' 
                % (self.fluence, self.arrival_time, self.intrinsic_width, self.dispersion_measure,
                   self.scattering_measure, self.spectral_index))


    property fluence:
        def __get__(self):
            assert self._pulse != NULL
            return self._pulse.fluence
        def __set__(self, x):
            assert self._pulse != NULL
            self._pulse.fluence = <double> x

    property spectral_index:
        def __get__(self):
            assert self._pulse != NULL
            return self._pulse.spectral_index
        def __set__(self, x):
            assert self._pulse != NULL
            self._pulse.spectral_index = <double> x

    property arrival_time:
        def __get__(self):
            assert self._pulse != NULL
            return self._pulse.arrival_time
        def __set__(self, x):
            assert self._pulse != NULL
            self._pulse.arrival_time = <double> x

    property intrinsic_width:
        def __get__(self):
            assert self._pulse != NULL
            return self._pulse.intrinsic_width
        def __set__(self, x):
            assert self._pulse != NULL
            self._pulse.intrinsic_width = <double> x

    property dispersion_measure:
        def __get__(self):
            assert self._pulse != NULL
            return self._pulse.dispersion_measure
        def __set__(self, x):
            assert self._pulse != NULL
            self._pulse.dispersion_measure = <double> x

    property scattering_measure:
        def __get__(self):
            assert self._pulse != NULL
            return self._pulse.scattering_measure
        def __set__(self, x):
            assert self._pulse != NULL
            self._pulse.scattering_measure = <double> x

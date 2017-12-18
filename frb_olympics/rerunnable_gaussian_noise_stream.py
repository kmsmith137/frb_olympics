import copy
import numpy as np

import rf_pipelines


class rerunnable_gaussian_noise_stream(rf_pipelines.wi_stream):
    """
    Similar to rf_pipelines.gaussian_noise_stream, but allows the stream to be rerun by saving its state:
    
        s = rerunnable_gaussian_noise_stream(...)
        saved_state = s.get_state()
           # ... run stream ...
        s.set_state(saved_state)
           # ... rerunning stream will give same output ...

    If 'no_noise_flag' is True, then the stream will output zeroes instead of Gaussian random numbers.
    """

    def __init__(self, nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, simulate_noise=True, state=None, nt_chunk=None):
        rf_pipelines.wi_stream.__init__(self, 'rereunnable_gaussian_noise_stream')

        if nt_tot <= 0:
            raise RuntimeError('rerunnable_gaussian_noise_stream constructor: nt_tot must be > 0')
        if nt_chunk is None:
            nt_chunk = min(1024, nt_tot)

        # These are members of the rf_pipelines.wi_stream base class.
        self.nfreq = nfreq
        self.nt_chunk = nt_chunk

        self.freq_lo_MHz = freq_lo_MHz
        self.freq_hi_MHz = freq_hi_MHz
        self.dt_sample = dt_sample
        self.simulate_noise = simulate_noise
        self.nt_tot = nt_tot
        self.set_state(state)


    def _bind_stream(self, json_attrs):
        json_attrs['freq_lo_MHz'] = self.freq_lo_MHz
        json_attrs['freq_hi_MHz'] = self.freq_hi_MHz
        json_attrs['dt_sample'] = self.dt_sample


    def _fill_chunk(self, intensity, weights, pos):
        assert intensity.shape == weights.shape == (self.nfreq, self.nt_chunk)

        if self.simulate_noise:
            intensity[:,:] = self.state.standard_normal((self.nfreq, self.nt_chunk))
        else:
            intensity[:,:] = np.zeros((self.nfreq, self.nt_chunk), dtype=np.float32)

        weights[:,:] = np.ones((self.nfreq, self.nt_chunk), dtype=np.float32)

        return (pos + self.nt_chunk) < self.nt_tot
        

    def get_state(self):
        """Returns the current RNG state."""
        return copy.copy(self.state)


    def set_state(self, state):
        """Restores the RNG state to a previous value, obtained by calling get_state()."""

        if state is None:
            self.state = np.random.RandomState()
        else:
            assert isinstance(state, np.random.RandomState)
            self.state = copy.copy(state)

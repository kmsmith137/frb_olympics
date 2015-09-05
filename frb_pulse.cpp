#include "frb_olympics.hpp"

#include <fftw3.h>
#include <boost/shared_array.hpp>

using namespace std;
using namespace boost;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


// -------------------------------------------------------------------------------------------------


struct frb_sampled_pulse {
    int     nt;   // number of samples
    double  dt;   // sample duration
    double  t0;   // start time; first sample covers time range (t0,t0+dt)
    shared_array<double> buf;

    // Note: make sure to allocate some zero padding, to avoid artifacts from periodic boundary conditions!
    frb_sampled_pulse(int nt, double t0, double dt);

    void downsample(int new_nt);
    void truncate(int new_nt);
};

frb_sampled_pulse::frb_sampled_pulse(int nt_, double dt_, double t0_)
    : nt(nt_), dt(dt_), t0(t0_)
{
    xassert(nt >= 2);
    xassert(dt > 0.0);
    buf = shared_array<double> (new double[nt]);
    memset(buf.get(), 0, nt * sizeof(double));
}

void frb_sampled_pulse::downsample(int new_nt)
{
    xassert(new_nt > 0);
    xassert(nt % new_nt == 0);

    shared_array<double> new_buf = shared_array<double> (new double[new_nt]);
    memset(new_buf.get(), 0, new_nt * sizeof(double));

    int m = nt / new_nt;
    for (int i = 0; i < new_nt; i++)
	for (int j = 0; j < m; j++)
	    new_buf[i] += buf[i*m+j] / (double)m;

    dt *= m;
    buf = new_buf;
}

void frb_sampled_pulse::truncate(int new_nt)
{
    xassert(new_nt > 0);
    xassert(new_nt <= nt);

    shared_array<double> new_buf = shared_array<double> (new double[new_nt]);
    memcpy(new_buf.get(), buf.get(), new_nt * sizeof(double));
    
    nt = new_nt;
    buf = new_buf;
}


// -------------------------------------------------------------------------------------------------


struct frb_pulse_fft {
    int     nt;
    double  dt;
    double  t0;
    shared_array<complex<double> > buf;   // length (nt/2+1)

    frb_pulse_fft(int nt, double dt, double t0);
    
    void  fill_with_gaussian(double arrival_time, double fluence, double intrinsic_width);
    void  apply_dispersion(double dt_dispersion);
    void  apply_scattering(double dt_scattering);

    frb_sampled_pulse fft();
};


frb_pulse_fft::frb_pulse_fft(int nt_, double dt_, double t0_)
    : nt(nt_), dt(dt_), t0(t0_)
{
    xassert(nt >= 2);
    xassert(dt > 0.0);
    buf = shared_array<complex<double> > (new complex<double> [nt/2+1]);
    memset(buf.get(), 0, (nt/2+1) * sizeof(complex<double>));
}


// Returns Fourier transformed dispersion kernel: (e^{-ikt} - 1) / (-ikt)
static inline complex<double> dispersion_fft(double kt)
{
    if (kt*kt < 1.0e-20)
	return complex<double> (1.0, -0.5*kt);   // first two terms in Taylor series

    complex<double> x(0, -kt);   // -ikt
    return (exp(x)-1.0)/x;
}


void frb_pulse_fft::fill_with_gaussian(double arrival_time, double fluence, double intrinsic_width)
{
    xassert(arrival_time >= t0);
    xassert(arrival_time <= t0 + nt*dt);
    xassert(fluence > 0.0);
    xassert(intrinsic_width >= 0.0);

    double w = intrinsic_width;
    double t = arrival_time - t0;

    // boxcar-averaging included as dispersion kernel
    for (int i = 0; i < nt/2+1; i++) {
	double k = 2*M_PI / (nt*dt) * i;
	buf[i] = fluence * exp(complex<double>(-0.5*w*w*k*k, -t*k)) * dispersion_fft(k*dt);
    }

    // Nyquist
    if (nt % 2 == 0)
	buf[nt/2] = complex<double>(0);
}

void frb_pulse_fft::apply_dispersion(double dt_dispersion)
{
    for (int i = 0; i < (nt+1)/2; i++) {
	double k = 2*M_PI / (nt*dt) * i;
	buf[i] *= dispersion_fft(k * dt_dispersion);
    }

    // Nyquist
    if (nt % 2 == 0)
	buf[nt/2] = complex<double>(0);
}

void frb_pulse_fft::apply_scattering(double dt_scattering)
{
    // Note: Fourier transformed scattering kernel is 1/(1+ikt)
    for (int i = 0; i < (nt+1)/2; i++) {
	double k = 2*M_PI / (nt*dt) * i;
	buf[i] /= complex<double>(1.0, k * dt_scattering);
    }

    // Nyquist
    if (nt % 2 == 0)
	buf[nt/2] = complex<double>(0);
}


frb_sampled_pulse frb_pulse_fft::fft()
{
    static const int max_rank = 18;   // 262144-element FFT is largest allowed
    static fftw_plan *p = NULL;

    if (!is_power_of_two(nt))
	throw runtime_error("only power-of-two FFT is supported for now");

    int rank = integer_log2(nt);
    if (rank > max_rank) {
	stringstream s;
	s << "FFT requested is too large (n=" << nt << ")";
	throw runtime_error(s.str());
    }

    frb_sampled_pulse ret(nt, dt, t0);
    fftw_complex *src = reinterpret_cast<fftw_complex *> (this->buf.get());
    double *dst = reinterpret_cast<double *> (ret.buf.get());

    if (!p) {
	p = reinterpret_cast<fftw_plan *> (fftw_malloc((max_rank+1) * sizeof(fftw_plan)));
	xassert(p != NULL);
	memset(p, 0, (max_rank+1) * sizeof(fftw_plan));
    }

    if (!p[rank]) {
	p[rank] = fftw_plan_dft_c2r_1d(nt, src, dst, FFTW_ESTIMATE);
	xassert(p[rank] != NULL);
    }

    fftw_execute_dft_c2r(p[rank], src, dst);
    memset(src, 0, (nt/2+1) * sizeof(complex<double>));
    
    for (int i = 0; i < nt; i++)
	dst[i] /= (dt*nt);
    
    return ret;
}


// -------------------------------------------------------------------------------------------------


//
// Helper function for frb_pulse::get_endpoints() and frb_pulse::get_signal_to_noise_in_channel()
// Simulates pulse in one frequency channel
//
//  @nt = length of returned timestream
//  @dt = sample size of returned timestream
//  @t0 = start time of returned timestream
//
static frb_sampled_pulse make_pulse(const frb_pulse &p, double freq_lo_MHz, double freq_hi_MHz, int nt, double dt, double t0)
{
    xassert(freq_lo_MHz <= freq_hi_MHz); 
    xassert(nt > 0);
    xassert(dt > 1.0e-4);
    
    // zero-pad by factor of two
    int nt_pad = round_up_to_power_of_two(2*nt);

    // oversampling factor
    // FIXME currently hardcoded to a large value; should depend on tscatt and tdisp
    int novs = 64;

    double freq_mid_MHz = (freq_lo_MHz + freq_hi_MHz) / 2.0;
    double tdisp0 = dispersion_delay(p.dispersion_measure, freq_hi_MHz);
    double tdisp1 = dispersion_delay(p.dispersion_measure, freq_lo_MHz);
    double tscatt = scatter_broadening(p.scattering_measure, freq_mid_MHz);
    double feff = p.fluence * pow(freq_mid_MHz/1000., p.spectral_index);

    frb_pulse_fft pf(nt_pad * novs, dt/novs, t0);
    pf.fill_with_gaussian(p.arrival_time + tdisp0, feff, p.intrinsic_width);
    pf.apply_dispersion(tdisp1 - tdisp0);
    pf.apply_scattering(tscatt);

    frb_sampled_pulse sp = pf.fft();
    sp.downsample(nt_pad);
    sp.truncate(nt);

    return sp;
}



// -------------------------------------------------------------------------------------------------


frb_pulse::frb_pulse(double fluence_, double arrival_time_, double intrinsic_width_,
		     double dispersion_measure_, double scattering_measure_, double spectral_index_)
    : fluence(fluence_), 
      spectral_index(spectral_index_), 
      arrival_time(arrival_time_),
      intrinsic_width(intrinsic_width_), 
      dispersion_measure(dispersion_measure_), 
      scattering_measure(scattering_measure_)
{
    xassert(fluence > 0.0);
    xassert(intrinsic_width >= 0.0);
    xassert(dispersion_measure >= 0.0);
    xassert(scattering_measure >= 0.0);
}


void frb_pulse::get_endpoints(double &t0, double &t1, double freq_lo_MHz, double freq_hi_MHz) const
{
    get_pulse_endpoints(t0, t1, this->intrinsic_width, this->dispersion_measure, this->scattering_measure, freq_lo_MHz, freq_hi_MHz);
    t0 += this->arrival_time;
    t1 += this->arrival_time;
}


void frb_pulse::add_to_timestream(double freq_lo_MHz, double freq_hi_MHz, float *timestream, int nsamples_per_chunk, double dt_sample, int ichunk) const
{
    xassert(freq_lo_MHz <= freq_hi_MHz); 
    xassert(nsamples_per_chunk > 0);
    xassert(ichunk >= 0);
    
    double t0, t1;
    this->get_endpoints(t0, t1, freq_lo_MHz, freq_hi_MHz);

    xassert(t0 > 0.0);
    xassert(t0 < t1);

    // pulse endpoints in integer samples
    int i0 = (int)(t0 / dt_sample);
    int i1 = (int)(t1 / dt_sample) + 1;

    // pulse endpoints clamped to current chunk
    int j0 = max(i0, ichunk*nsamples_per_chunk);
    int j1 = min(i1, (ichunk+1)*nsamples_per_chunk);

    if (j0 >= j1)
	return;   // pulse does not overlap current chunk

    frb_sampled_pulse sp = make_pulse(*this, freq_lo_MHz, freq_hi_MHz, i1-i0, dt_sample, i0*dt_sample);

    for (int j = j0; j < j1; j++)
	timestream[j - ichunk*nsamples_per_chunk] += sp.buf[j-i0];
}


// FIXME still a little cut-and-paste with frb_pulse::add_to_timestream() to clean up
double frb_pulse::get_signal_to_noise_in_channel(double freq_lo_MHz, double freq_hi_MHz, double dt_sample) const
{
    xassert(freq_lo_MHz <= freq_hi_MHz);
    xassert(dt_sample > 1.0e-4);
    
    double t0, t1;
    this->get_endpoints(t0, t1, freq_lo_MHz, freq_hi_MHz);

    xassert(t0 > 0.0);
    xassert(t0 < t1);

    // pulse endpoints in integer samples
    int i0 = (int)(t0 / dt_sample);
    int i1 = (int)(t1 / dt_sample) + 1;

    frb_sampled_pulse sp = make_pulse(*this, freq_lo_MHz, freq_hi_MHz, i1-i0, dt_sample, i0*dt_sample);
    
    double acc = 0.0;
    for (int i = 0; i < i1-i0; i++)
	acc += square(sp.buf[i]);

    return sqrt(acc);
}


}  // namespace frb_olympics

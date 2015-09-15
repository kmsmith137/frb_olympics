#include "frb_olympics.hpp"
#include <fstream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


void get_pulse_endpoints(double &t0, double &t1, double intrinsic_width, double dm, double sm, double freq_lo_MHz, double freq_hi_MHz)
{
    xassert(freq_lo_MHz <= freq_hi_MHz);

    t0 = dispersion_delay(dm,freq_hi_MHz) - 5*intrinsic_width;
    t1 = dispersion_delay(dm,freq_lo_MHz) + 5*intrinsic_width + 12.5*scatter_broadening(sm,freq_lo_MHz);
} 


void tokenize_file(const string &filename, vector<vector<string> > &tokens)
{
    string line;
    vector<string> line_tokens;

    tokens.clear();

    ifstream f(filename.c_str());

    if (f.fail()) {
	stringstream s;
	s << ": couldn't open file " << filename << "\n";
	throw runtime_error(s.str().c_str());
    }

    for (;;) {
	if (!std::getline(f, line))
	    return;

	// remove comment starting with '#'
	string::size_type i = line.find_first_of("#");
	if (i != string::npos)
	    line.resize(i);

	// remove leading and trailing spaces
	boost::trim(line);

	if (line.size() == 0)
	    continue;

	boost::split(line_tokens, line, boost::is_space(), boost::token_compress_on);

	if (line_tokens.size() > 0)
	    tokens.push_back(line_tokens);
    }
}


void read_kv_pairs(const string &filename, boost::unordered_map<string,string> &kv_pairs)
{
    kv_pairs.clear();

    vector<vector<string> > tokens;
    tokenize_file(filename, tokens);

    for (unsigned int i = 0; i < tokens.size(); i++) {
	if ((tokens[i].size() != 3) || (tokens[i][1] != "=")) {
	    stringstream s;

	    s << filename << ": parse error in line:";

	    for (unsigned int j = 0; j < tokens[i].size(); j++)
		s << " " << tokens[i][j];

	    throw runtime_error(s.str().c_str());
	}

	if (kv_pairs.find(tokens[i][0]) != kv_pairs.end()) {
	    stringstream s;
	    s << filename << ": key '" << tokens[i][0] << "' occurs twice";
	    throw runtime_error(s.str().c_str());
	}

	kv_pairs[tokens[i][0]] = tokens[i][2];
    }
}

bool is_power_of_two(int n)
{
    xassert(n >= 1);
    return (n & (n-1)) == 0;
}

int round_up_to_power_of_two(int n)
{
    xassert(n >= 1);

    int ret = 1;
    while (ret < n) {
	ret <<= 1;
	xassert(ret > 0);  // detect signed overflow
    }

    return ret;
}

int integer_log2(int n)
{
    if (!is_power_of_two(n))
	throw runtime_error("integer_log2 called with non-power-of-two argument");

    int ret = 0;
    while ((1 << ret) < n)
	ret++;

    xassert(n == (1 << ret));
    return ret;
}

static inline bool time_before(const struct timeval &tv1, const struct timeval &tv2)
{
    if (tv1.tv_sec < tv2.tv_sec)
	return true;
    if (tv1.tv_sec > tv2.tv_sec)
	return false;
    return (tv1.tv_usec <= tv2.tv_usec);   // returns true if times are equal
}

double time_diff(const struct timeval &tv1, const struct timeval &tv2)
{
    xassert(time_before(tv1, tv2));
    return (tv2.tv_sec - tv1.tv_sec) + 1.0e-6 * (tv2.tv_usec - tv1.tv_usec);
}

struct timeval get_time()
{
    struct timeval ret;
    if (gettimeofday(&ret, NULL) < 0)
	throw runtime_error("gettimeofday() failed");
    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// correlation function utils


//
// This gives the correlation function of a scattered timestream.
//
// @a is a number between 0 and 1 parametrizing the strength of the scattering
// as follows.  The scattered timestream is given by 
//    S(t) = (1-a) * (U(t) + a U(t+1) + a^2 U(t+2) + ...) 
//
// where U(t) is the unscattered timestream, and times are in "upsampled" units.
// This is the same quantity returned by ascatt(), provided that ascatt() is called
// with the _upsampled_ dt_sample.
//
// The output is a shape (nups,lag) array containing the correlation funciton
//     zeta[m,n] = < S(m) S(m+n) >
//
void init_cf(double *out, int nups, int lag, double a)
{
    xassert(nups >= 1);
    xassert(lag >= 1);
    xassert(a >= 0.0 && a < 1.0);

    int maxpow = 2*nups + lag;
    vector<double> apow(maxpow + 1);

    apow[0] = 1.0;
    for (int i = 0; i < maxpow; i++)
	apow[i+1] = apow[i] * a;

    for (int s = 0; s < nups; s++) {
	for (int t = 0; t < lag; t++) {
	    // start and end of frame containing (s+t)
	    int i = ((s+t) / nups) * nups;
	    int j = i + nups;

	    double w0 = 1 - pow(a,nups);
	    double w1 = (i > s) ? (w0 * pow(a,i-s)) : (1 - pow(a,j-s));
	    double w2 = 1 - pow(a,j-s-t);

	    xassert(2*j-2*s-t <= maxpow);
	    out[s*lag+t] = w1*w2 + w0*w0 * pow(a,2*j-2*s-t) / (1 - pow(a,2*nups));
	}
    }
}


//
// Applies a convolution kernel to a correlation function.
//
// The input @in is a shape (nups, lag_in) array containing the correlation function
//   zeta_in[m,n] = < U(m) U(m+n) >
//
// where U is an unconvolved timestream.  The input @wconv is a shape (nconv,) array
// containing a convolution kernel K, such that the convolved timestream C is related
// to the unconvolved timestream U by
//   C(m) = sum_{t=0}^{nconv-1} W(t) U(m+t)
//
// The output @out is a shape (nups, lag_out) array containing the correlation function
//   zeta_out[m,n] = < C(m) C(m+n) >
//
// This routine is written in a brute-force way which scales as O(nconv^2) (with other input
// parameters fixed).  This could be improved to O(nconv), but I usually call it with small
// nconv, so I haven't done this yet.
//
void convolve_cf(double *out, int nups, int lag_out, int lag_in, int nconv, const double *wconv, const double *in)
{
    xassert(lag_in >= lag_out + nconv - 1);

    memset(out, 0, nups * lag_out * sizeof(*out));

    for (int s = 0; s < nups; s++) {
	for (int t = 0; t < lag_out; t++) {
	    for (int i = 0; i < nconv; i++) {
		for (int j = 0; j < nconv; j++) {
		    int ii = min(s+i, s+t+j);
		    int jj = max(s+i, s+t+j);
		    int ss = ii % nups;
		    int tt = jj - ii;

		    out[s*lag_out + t] += wconv[i] * wconv[j] * in[ss*lag_in + tt];
		}
	    }
	}
    }
}

// -------------------------------------------------------------------------------------------------
//
// unit tests for correlation function utils


static vector<double> randvec(frb_rng &r, int n)
{
    vector<double> ret(n);
    for (int i = 0; i < n; i++)
	ret[i] = r.gaussian();
    return ret;
}

static double compare(int n, const double *v1, const double *v2)
{
    double num = 0.0;
    double den = 0.0;

    for (ssize_t i = 0; i < n; i++) {
	num += square(v1[i] - v2[i]);
	den += square(v1[i]) + square(v2[i]);
    }

    return (den > 0.0) ? sqrt(num/den) : 0.0;
}

static double compare(const vector<double> &v1, const vector<double> &v2)
{
    xassert(v1.size() == v2.size());
    return compare(v1.size(), &v1[0], &v2[0]);
}


static void init_cf_slow(double *out, int nups, int lag, double a)
{
    xassert(nups >= 1);
    xassert(lag >= 1);
    xassert(a >= 0 && a < 1.0);

    int nt = (int)(log(1.0e-20)/log(a)) + 2;

    int nrows = nups + lag - 1;
    int ncols = (nrows + nt) / nups + 1;
    vector<double> w(nrows * ncols, 0.0);

    for (int i = 0; i < nrows; i++) {
	for (int j = 0; j < nt; j++) {
	    int s = (i+j) / nups;
	    xassert(s >= 0 && s < ncols);
	    w[i*ncols + s] += (1-a) * pow(a,j);
	}
    }

    memset(out, 0, nups * lag * sizeof(*out));

    for (int s = 0; s < nups; s++)
	for (int t = 0; t < lag; t++)
	    for (int i = 0; i < ncols; i++)
		out[s*lag+t] += w[s*ncols+i] * w[(s+t)*ncols+i];
}


static void test_init_cf(int nups, int lag, double a)
{
    vector<double> cf1(nups*lag, 0.0);
    vector<double> cf2(nups*lag, 0.0);

    init_cf(&cf1[0], nups, lag, a);
    init_cf_slow(&cf2[0], nups, lag, a);

    cout << "test_init_cf(nups=" << nups << ",lag=" << lag 
	 << ",a=" << a << "): " << compare(cf1,cf2) << endl;
}


static void test_convolve_associativity(int nups, int lag, int nconv1, int nconv2)
{
    frb_rng rng;

    int lag1 = lag + nconv1 - 1;
    int lag12 = lag + nconv1 + nconv2 - 2;

    vector<double> cf = randvec(rng, nups * lag12);
    vector<double> wconv1 = randvec(rng, nconv1);
    vector<double> wconv2 = randvec(rng, nconv2);

    vector<double> t1(nups * lag1);
    vector<double> t2(nups * lag);
    convolve_cf(&t1[0], nups, lag1, lag12, nconv2, &wconv2[0], &cf[0]);
    convolve_cf(&t2[0], nups, lag, lag1, nconv1, &wconv1[0], &t1[0]);

    int nconv12 = nconv1 + nconv2 - 1;
    vector<double> wconv12(nconv12, 0.0);
    for (int i = 0; i < nconv1; i++)
	for (int j = 0; j < nconv2; j++)
	    wconv12[i+j] += wconv1[i] * wconv2[j];

    convolve_cf(&t1[0], nups, lag, lag12, nconv12, &wconv12[0], &cf[0]);

    cout << "test_convolve_associativity(" << nups << "," << lag << "," << nconv1 
	 << "," << nconv2 << "): " << compare(lag, &t1[0], &t2[0]) << endl;
}

static void test_init_convolve_consistency(int nups, int lag, double a)
{
    int nconv = (int)(log(1.0e-20)/log(a)) + 2;
    int lag0 = lag + nconv - 1;

    vector<double> wconv(nconv);
    for (int i = 0; i < nconv; i++)
	wconv[i] = (1-a) * pow(a,i);
    
    vector<double> cf0(nups * lag0);
    for (int s = 0; s < nups; s++)
	for (int t = 0; t < lag0; t++)
	    cf0[s*lag0+t] = ((s+t) < nups) ? 1.0 : 0.0;

    vector<double> cf1(nups * lag);
    convolve_cf(&cf1[0], nups, lag, lag0, nconv, &wconv[0], &cf0[0]);
    
    vector<double> cf2(nups * lag);
    init_cf(&cf2[0], nups, lag, a);

    cout << "test_init_convolve_consistency(" << nups << "," << lag 
	 << "," << a << "): " << compare(cf1, cf2) << endl;
}


void run_cf_unit_tests()
{
    test_init_cf(5, 12, 0.72321);
    test_convolve_associativity(7, 14, 17, 21);
    test_init_convolve_consistency(6, 11, 0.7813672);
}


}  // namespace frb_olympics

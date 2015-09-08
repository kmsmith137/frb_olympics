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


}  // namespace frb_olympics

#include "frb_olympics.hpp"

using namespace std;

namespace frb_olympics {
#if 0
} // emacs pacifier
#endif


frb_search_algorithm_base::frb_search_algorithm_base(const frb_search_params &p_)
    : p(p_), memhack(1), name(), debug_buffer_ndm(0), debug_buffer_nt(0), search_gb(0.0)
{ }


// --------------------  Registry stuff  --------------------


struct reg_entry {
    frb_search_algorithm_base::creator_t creator;
    string usage;

    reg_entry(frb_search_algorithm_base::creator_t creator_, const string &usage_)
	: creator(creator_), usage(usage_) { }
};

typedef std::map<string,reg_entry> registry_t;

static registry_t &get_registry()
{
    static registry_t *registry = NULL;

    if (!registry)
	registry = new registry_t;
    return (*registry);
}

// static member function
void frb_search_algorithm_base::register_algorithm(const string &algo_name, creator_t f, const string &usage)
{
    registry_t &registry = get_registry();

    if (registry.find(algo_name) != registry.end()) {
	cerr << "fatal: search algorithm name '" << algo_name << "' was registered twice\n";
	exit(1);
    }
    
    registry.insert(registry_t::value_type(algo_name, reg_entry(f,usage)));
}

// static member function
void frb_search_algorithm_base::show_algorithm_usage(ostream &os)
{
    registry_t &registry = get_registry();
    os << "Syntax of registered algorithms follows:\n";

    for (registry_t::iterator pp = registry.begin(); pp != registry.end(); pp++)
	os << "\n" << pp->second.usage;

    os << "\n"
       << "    Putting [memhack N] anywhere on the line will increase available memory per MPI task by a factor N,\n"
       << "    by running in a mode in which each node leaves (N-1) out of N cores idle" << endl;
}

// static member function
frb_search_algorithm_base::ptr_t frb_search_algorithm_base::parse_line(const vector<string> &tokens_, const frb_search_params &p)
{
    vector<string> tokens = tokens_;

    int memhack = 1;
    for (unsigned int i = 0; i < tokens.size()-1; i++) {
	if (tokens[i] != string("memhack"))
	    continue;

	memhack = xlexical_cast<int> (tokens[i+1], "memhack");
	assert(memhack >= 1);

	tokens.erase(tokens.begin()+i, tokens.begin()+i+2);
	break;
    }

    assert(tokens.size() > 0);
    string algo_name = tokens[0];
    tokens.erase(tokens.begin());

    registry_t &registry = get_registry();
    registry_t::iterator pp = registry.find(algo_name);

    if (pp == registry.end()) {
	cerr << "fatal: search algorithm '" << algo_name << "' requested but not found in registry\n";
	exit(1);
    }

    frb_search_algorithm_base::ptr_t ret = pp->second.creator(tokens, p);
    
    if (!ret) {
	cerr << "parse error on line:";
	for (unsigned int i = 0; i < tokens_.size(); i++)
	    cerr << " " << tokens_[i];
	cerr << "\n" << algo_name << " command-line syntax is:\n" << pp->second.usage;
	exit(1);
    }

    return ret;
}


// static member function
void frb_search_algorithm_base::read_file(const string &filename, const frb_search_params &p, vector<ptr_t> &algo_list)
{
    algo_list.clear();

    vector<vector<string> > tokens;
    tokenize_file(filename, tokens);

    for (unsigned int i = 0; i < tokens.size(); i++)
	algo_list.push_back(parse_line(tokens[i],p));
}


}  // namespace frb_olympics

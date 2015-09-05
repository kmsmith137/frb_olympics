#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <fstream>
#include "frb_olympics.hpp"

using namespace std;
using namespace frb_olympics;


// -------------------------------------------------------------------------------------------------
//
// Node info (used in memhack)


static int mpi_tasks_per_node = 0;
static int mpi_rank_within_node = 0;


static void init_node_info()
{
    if (mpi_tasks_per_node > 0)
	return;   // already initialized

    int mpi_rank = MPI::COMM_WORLD.Get_rank();
    int mpi_size = MPI::COMM_WORLD.Get_size();

    vector<char> hname_buf1(HOST_NAME_MAX+2, 0);
    vector<char> hname_buf2(mpi_size * (HOST_NAME_MAX+2), 0);

    if (gethostname(&hname_buf1[0], HOST_NAME_MAX+1) < 0)
	throw runtime_error("rng_initializer: gethostname() failed");

    MPI::COMM_WORLD.Allgather(&hname_buf1[0], HOST_NAME_MAX+2, MPI::CHAR,
			      &hname_buf2[0], HOST_NAME_MAX+2, MPI::CHAR);

    string my_hname = string(&hname_buf1[0]);
    boost::unordered_map<string,int> hname_counts;

    for (int i = 0; i < mpi_size; i++) {
	string s = string(&hname_buf2[i * (HOST_NAME_MAX+2)]);

	if (s == my_hname) {
	    mpi_tasks_per_node++;
	    if (i < mpi_rank)
		mpi_rank_within_node++;
	}

	boost::unordered_map<string,int>::iterator p = hname_counts.find(s);
	int c = (p == hname_counts.end()) ? 0 : p->second;
	hname_counts[s] = c+1;
    }

    // check
    for (boost::unordered_map<string,int>::iterator p = hname_counts.begin(); p != hname_counts.end(); p++) {
	if (p->second != mpi_tasks_per_node)
	    throw runtime_error("init_node_info: expected all hostnames to have equal numbers of MPI tasks");
    }

    cout << "rank_within_node=" << mpi_rank_within_node << ", tasks_per_node=" << mpi_tasks_per_node << endl;
}


static void memhack_start(int memhack)
{
    if (memhack == 1)
	return;

    init_node_info();   // no-ops if already initialized
    xassert(memhack >= 1);
    xassert(mpi_tasks_per_node % memhack == 0);

    int nbarriers = (mpi_rank_within_node % memhack) + 1;
    for (int i = 0; i < nbarriers; i++)
	MPI::COMM_WORLD.Barrier();
}


static void memhack_end(int memhack)
{
    if (memhack == 1)
	return;

    init_node_info();   // no-ops if already initialized
    xassert(memhack >= 1);
    xassert(mpi_tasks_per_node % memhack == 0);

    int nbarriers = memhack - (mpi_rank_within_node % memhack);
    for (int i = 0; i < nbarriers; i++)
	MPI::COMM_WORLD.Barrier();
}


// -------------------------------------------------------------------------------------------------


static void usage()
{
    cerr << "usage: frb-compare <search_params.txt> <algo_list.txt> <nmc_noise> <nmc_pulse> <output_stem> [target_sn]\n"
	 << "       target_sn defaults to 30.0\n"
	 << "\n";

    frb_search_algorithm_base::show_algorithm_usage(cerr);
    exit(2);
}


int main(int argc, char **argv)
{
    MPI::Init(argc, argv);

    if ((argc < 6) || (argc > 7))
	usage();

    int mpi_rank = MPI::COMM_WORLD.Get_rank();
    int mpi_size = MPI::COMM_WORLD.Get_size();
    
    const char *search_params_filename = argv[1];
    const char *algo_list_filename = argv[2];
    int nmc_noise_tot = xlexical_cast<int> (argv[3], "frb-compare command line");
    int nmc_pulse_tot = xlexical_cast<int> (argv[4], "frb-compare command line");
    const char *output_stem = argv[5];
    const double target_sn = (argc == 7) ? xlexical_cast<double>(argv[6],"target_sn") : 30.0;

    // Redirect stdout to log file which is unique to this MPI rank
    if (mpi_size > 1) {
	stringstream slog;
	slog << output_stem << "_log" << mpi_rank;

	string log_filename = slog.str();
	cout << "redirecting stdout -> " << log_filename << endl;

	int log_fd = open(log_filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
	
	if (log_fd < 0) {
	    cerr << "couldn't create log file " << log_filename << ": " << strerror(errno) << endl;
	    exit(1);
	}

	// stdout is file descriptor 1
	if (dup2(log_fd, 1) < 0)
	    cerr << "dup2() failed: " << strerror(errno) << endl;
    }

    cout << "command line:";
    for (int i = 0; i < argc; i++) 
	cout << " " << argv[i];
    cout << endl;

    xassert(nmc_noise_tot >= 0);
    xassert(nmc_pulse_tot >= 0);
    xassert(nmc_noise_tot % mpi_size == 0);
    xassert(nmc_pulse_tot % mpi_size == 0);

    int nmc_noise_loc = nmc_noise_tot / mpi_size;
    int nmc_pulse_loc = nmc_pulse_tot / mpi_size;

    frb_rng rng;
    frb_search_params search_params(search_params_filename);

    vector<frb_search_algorithm_base::ptr_t> algo_list;
    frb_search_algorithm_base::read_file(algo_list_filename, search_params, algo_list);
    int nalgo = algo_list.size();

    ofstream f_noise;
    ofstream f_pulse;

    if ((mpi_rank == 0) && (nmc_noise_tot + nmc_pulse_tot > 0)) {
	//
	// Open output files
	// Note: filename convention must be kept in sync with frb-compare-postprocess.py
	//
	// FIXME should move files out of the way if they already exist (currently silently clobbers)
	//
	stringstream s1;
	s1 << output_stem << "_noise.txt";

	string filename1 = s1.str();
	f_noise.open(filename1.c_str());

	if (f_noise.fail()) {
	    cerr << "Fatal: couldn't open file " << filename1 << endl;
	    exit(1);
	}

	stringstream s2;
	s2 << output_stem << "_pulse.txt";

	string filename2 = s2.str();
	f_pulse.open(filename2.c_str());

	if (f_pulse.fail()) {
	    cerr << "Fatal: couldn't open file " << filename2 << endl;
	    exit(1);
	}
    }

    cout << "--------------------  Search params  --------------------\n";
    search_params.write(cout);
    cout << endl;

    cout << "--------------------  Memory analysis  --------------------\n";

    double timestream_gb = 1.0e-9 * search_params.nchan * search_params.nsamples_per_chunk * sizeof(float);
    cout << "timestream chunk size = " << timestream_gb << " GB\n";

    double search_gb = 0.0;
    for (int ialgo = 0; ialgo < nalgo; ialgo++) {
	frb_search_algorithm_base::ptr_t algo = algo_list[ialgo];
	cout << "    " << algo->name << " " << algo->search_gb << " GB\n";
	search_gb = max(search_gb, algo->search_gb);
    }

    cout << "total size = " << (timestream_gb + search_gb) << " GB\n";

    vector<double> score_noise(nmc_noise_loc * nalgo);
    vector<double> score_pulse(nmc_pulse_loc * nalgo);
    vector<double> pulse_tab(nmc_pulse_loc * 5);

    // --------------------  noise sims  --------------------

    for (int inoise = 0; inoise < nmc_noise_loc; inoise++) {
	cout << "starting noise sim " << inoise << "/" << nmc_noise_loc << endl;

	vector<float> chunk(search_params.nchan * search_params.nsamples_per_chunk);
	frb_rng rsave(rng);   // save

	for (int ialgo = 0; ialgo < nalgo; ialgo++) {
	    frb_search_algorithm_base::ptr_t algo = algo_list[ialgo];
	    struct timeval tv0 = get_time();

	    memhack_start(algo->memhack);
	    algo->search_start();

	    rng = rsave;     // restore
	    for (int ichunk = 0; ichunk < search_params.nchunks; ichunk++) {
		cout << "     " << algo->name 
		     << ": starting chunk " << ichunk << "/" << search_params.nchunks 
		     << "  [time=" << time_diff(tv0,get_time()) << " secs]" << endl;

		search_params.simulate_noise(rng, &chunk[0]);
		algo->search_chunk(&chunk[0], ichunk, NULL);
	    }

	    cout << "     " << algo->name 
		 << ": search_result = " << algo->search_result 
		 << "  [time=" << time_diff(tv0,get_time()) << " secs]" << endl;

	    score_noise[inoise*nalgo + ialgo] = algo->search_result;

	    algo->search_end();
	    memhack_end(algo->memhack);
	}
    }

    // --------------------  pulse sims  --------------------

    for (int ipulse = 0; ipulse < nmc_pulse_loc; ipulse++) {
	cout << "starting pulse sim " << ipulse << "/" << nmc_pulse_loc << endl;

	frb_pulse p = search_params.make_random_pulse(rng, 1.0);
	pulse_tab[5*ipulse] = p.arrival_time;
	pulse_tab[5*ipulse+1] = p.intrinsic_width;
	pulse_tab[5*ipulse+2] = p.dispersion_measure;
	pulse_tab[5*ipulse+3] = p.scattering_measure;
	pulse_tab[5*ipulse+4] = p.spectral_index;

	double sn = search_params.get_signal_to_noise_of_pulse(p);
	p.fluence *= (target_sn / sn);
	
	vector<float> chunk(search_params.nchan * search_params.nsamples_per_chunk);	

	for (int ialgo = 0; ialgo < nalgo; ialgo++) {
	    frb_search_algorithm_base::ptr_t algo = algo_list[ialgo];	    
	    struct timeval tv0 = get_time();

	    memhack_start(algo->memhack);
	    algo->search_start();

	    for (int ichunk = 0; ichunk < search_params.nchunks; ichunk++) {
		cout << "     " << algo->name 
		     << ": starting chunk " << ichunk << "/" << search_params.nchunks 
		     << "  [time=" << time_diff(tv0,get_time()) << " secs]" << endl;

		memset(&chunk[0], 0, chunk.size() * sizeof(chunk[0]));
		search_params.add_pulse(p, &chunk[0], ichunk);
		algo->search_chunk(&chunk[0], ichunk, NULL);
	    }

	    cout << "     " << algo->name 
		 << ": search_result = " << algo->search_result 
		 << "  [time=" << time_diff(tv0,get_time()) << " secs]" << endl;
	    
	    score_pulse[ipulse*nalgo + ialgo] = algo->search_result;

	    algo->search_end();
	    memhack_end(algo->memhack);
	}
    }

    vector<double> score_noise_all(nmc_noise_tot * nalgo, 0.0);
    vector<double> score_pulse_all(nmc_pulse_tot * nalgo);
    vector<double> pulse_tab_all(nmc_pulse_tot * 5, 0.0);

    MPI::COMM_WORLD.Gather(&score_noise[0], nmc_noise_loc * nalgo, MPI::DOUBLE,
			   &score_noise_all[0], nmc_noise_loc * nalgo, MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Gather(&score_pulse[0], nmc_pulse_loc * nalgo, MPI::DOUBLE,
			   &score_pulse_all[0], nmc_pulse_loc * nalgo, MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Gather(&pulse_tab[0], 5 * nmc_pulse_loc, MPI::DOUBLE,
			   &pulse_tab_all[0], 5 * nmc_pulse_loc, MPI::DOUBLE, 0);

    if (mpi_rank == 0) {
	for (int i = 0; i < nmc_noise_tot; i++) {
	    for (int j = 0; j < nalgo; j++)
		f_noise << " " << score_noise_all[i*nalgo+j];
	    f_noise << "\n";
	}

	for (int i = 0; i < nmc_pulse_tot; i++) {
	    for (int j = 0; j < 5; j++)
		f_pulse << " " << pulse_tab_all[5*i+j];
	    for (int j = 0; j < nalgo; j++)
		f_pulse << " " << score_pulse_all[i*nalgo+j];
	    f_pulse << "\n";
	}
    }

    MPI::Finalize();
    return 0;
}

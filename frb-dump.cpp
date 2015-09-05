#include "frb_olympics.hpp"
#include <hdf5.h>

// backwards compatibility hacks for HDF5 1.6                                                                                                                                                                      
#if H5_VERS_MINOR == 6
#  define H5Aiterate1 H5Aiterate
#  define H5Acreate1 H5Acreate
#  define H5Dopen1 H5Dopen
#  define H5Dcreate1 H5Dcreate
#  define H5Gcreate1 H5Gcreate
#  define H5Eset_auto1 H5Eset_auto
#  define H5Ewalk1 H5Ewalk
#  define H5E_error1_t H5E_error_t
#  define H5Gopen1 H5Gopen
#endif

using namespace std;
using namespace frb_olympics;


static void usage()
{
    cerr << "usage: frb-dump <search_params.txt> <outfile.hdf5> <algo...>\n";
    cerr << "\n";
    frb_search_algorithm_base::show_algorithm_usage(cerr);
    exit(2);
}


int main(int argc, char **argv)
{
    if (argc < 4)
	usage();

    frb_search_params search_params(argv[1]);
    const char *hdf5_filename = argv[2];

    vector<string> tokens;
    for (int i = 3; i < argc; i++)
	tokens.push_back(argv[i]);

    frb_search_algorithm_base::ptr_t algo = frb_search_algorithm_base::parse_line(tokens, search_params);
    cout << "Debug buffer shape = (" << algo->debug_buffer_ndm << "," << algo->debug_buffer_nt << ")" << endl;

    double gb = algo->search_gb;
    gb += 1.0e-9 * search_params.nchan * search_params.nsamples_per_chunk * sizeof(float);  // timestream
    gb += 1.0e-9 * algo->debug_buffer_ndm * algo->debug_buffer_nt * sizeof(float);
    cout << "Estimated memory usage: " << gb << " GB" << endl;

    vector<double> dm_list;
    vector<frb_pulse> pulse_list;
    
    double dm0 = search_params.dm_min;
    double dm1 = search_params.dm_max;
    dm_list.push_back(0.97*dm0 + 0.03*dm1);
    dm_list.push_back(0.50*dm0 + 0.50*dm1);
    dm_list.push_back(0.03*dm0 + 0.97*dm1);

    for (unsigned int i = 0; i < dm_list.size(); i++) {
	double t0, t1;
	search_params.get_allowed_arrival_times(t0, t1, 1.0e-5, dm_list[i], 0.0);
	pulse_list.push_back(frb_pulse(1.0, 0.97*t0 + 0.03*t1, 1.0e-5, dm_list[i], 0.0, 0.0));
	pulse_list.push_back(frb_pulse(1.0, 0.50*t0 + 0.50*t1, 1.0e-5, dm_list[i], 0.0, 0.0));
	pulse_list.push_back(frb_pulse(1.0, 0.03*t0 + 0.97*t1, 1.0e-5, dm_list[i], 0.0, 0.0));
    }

    vector<float> debug_buffer((size_t)algo->debug_buffer_ndm * (size_t)algo->debug_buffer_nt, -1.0e30);

    algo->search_start();

    for (int ichunk = 0; ichunk < search_params.nchunks; ichunk++) {
	cout << algo->name << " processing chunk " << ichunk << "/" << search_params.nchunks << endl;

	vector<float> timestream(search_params.nchan * search_params.nsamples_per_chunk, 0.0);
	for (unsigned int i = 0; i < pulse_list.size(); i++)
	    search_params.add_pulse(pulse_list[i], &timestream[0], ichunk);

	algo->search_chunk(&timestream[0], ichunk, &debug_buffer[0]);
    }

    algo->search_end();
    
    cout << "writing " << hdf5_filename << endl;
    
    vector<hsize_t> shape(2);
    shape[0] = algo->debug_buffer_ndm;
    shape[1] = algo->debug_buffer_nt;

    hid_t file_id = H5Fcreate(hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
	cerr << "Fatal: couldn't create HDF5 file " << hdf5_filename << endl;
	exit(1);
    }

    hid_t group_id = H5Gopen1(file_id, ".");
    xassert(group_id >= 0);

    hid_t dataspace_id = H5Screate(H5S_SIMPLE);
    xassert(dataspace_id >= 0);

    int ret = H5Sset_extent_simple(dataspace_id, 2, &shape[0], &shape[0]);
    xassert(ret >= 0);

    hid_t dataset_id = H5Dcreate1(group_id, "DATA", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
    xassert(dataset_id >= 0);

    ret = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &debug_buffer[0]);
    xassert(ret >= 0);

    H5Fclose(file_id);
    return 0;
}

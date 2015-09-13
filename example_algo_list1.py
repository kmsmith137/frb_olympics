# A minimal algo_list.py file, which will compare the bonsai search algorithm
# with and without downsampling.  (More discussion in the README file.)

# Could also import sloth, simple_direct, simple_tree, fdmt
# (See example_algo_list2.py, which includes an example of each algorithm.)
from frb_olympics import add_algo, bonsai, downsample

#
# Some search algorithms are written in C++ and you'll need to refer to frb_olympics_c.pyx
# to see their constructor syntax.  For example, the bonsai source code can be found in
# frb_bonsai.cpp and its constructor syntax (from frb_olympics_c.pyx) is
#    bonsai(ntree, nupsample=1)
#
add_algo(bonsai(1024))

#
# Some algorithms are written in python, e.g. the downsample "algorithm", which is
# just a wrapper which takes an existing search algorithm and passes it downsampled data.
# See frb_downsample.py for source code, including constructor syntax.
#
add_algo(downsample(bonsai(512), 2))

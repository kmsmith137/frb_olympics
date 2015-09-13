# A reference algo_list.py file, which includes one example of each search algorithm.
#
# This will be slow to run, due to use of the direct search algorithms (sloth,simple_tree)!
# I suggest either commenting out the direct search algorithms, or running on a large machine 
# with MPI (either a many-core machine or a cluster),

from frb_olympics import add_algo, simple_direct, simple_tree, sloth, bonsai, downsample, fdmt

# constructor syntax: simple_direct(epsilon)
# epsilon is a parameter controlling the spacing of the DM table
add_algo(simple_direct(0.1))

# constructor syntax: sloth(epsilon, nupsample=1)
add_algo(sloth(0.1, nupsample=2))

# constructor syntax: simple_tree(ntree, ndownsample=1)
# The 'downsample' optional arg isn't really necessary since one can use the downsample() wrapper, see below.
add_algo(simple_tree(1024))

# constructor syntax: bonsai(ntree, nupsample=1)
add_algo(bonsai(1024))

# constructor syntax: fdmt(maxDT, recdepth, doWeighting=True)
add_algo(fdmt(256,2))

# example of using the downsample() wrapper, in this case to run the bonsai algorithm with factor-of-two downsampling.
add_algo(downsample(bonsai(512), 2))

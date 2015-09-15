# Note: same algo list used for all 3 runs in this memo

from frb_olympics import add_algo, sloth, sloth_sm_subsearch

# Let's try to make something pretty optimal
add_algo(sloth(0.1, nupsample=2))

# What happens when we assume SM=beta=0?
add_algo(sloth_sm_subsearch(0.0, 0.1, 1.0e10, nupsample=2))

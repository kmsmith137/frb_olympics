from frb_olympics import add_algo, sloth

# This is very optimal
add_algo(sloth(0.1,nupsample=2))

# How suboptimal are we with no upsampling?
add_algo(sloth(0.1,nupsample=1))

# this seems to get ~90% !
add_algo(sloth(1.0,nupsample=1))

# just for fun
add_algo(sloth(4.0,nupsample=1))

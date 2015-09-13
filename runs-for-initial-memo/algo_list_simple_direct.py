from frb_olympics import add_algo, sloth, simple_direct, downsample

add_algo(sloth(0.1,nupsample=2))

# isolate effect of downsampling, without decreasing epsilon
#  simple_direct 0.1
#  downsample 2 simple_direct 0.1
#  downsample 4 simple_direct 0.1

add_algo(simple_direct(1.0))
add_algo(downsample(simple_direct(1.0),2))
add_algo(downsample(simple_direct(1.0),4))


from frb_olympics import add_algo, sloth, simple_tree, downsample

add_algo(sloth(0.1,nupsample=2))
add_algo(simple_tree(1024))
add_algo(downsample(simple_tree(512),2))
add_algo(downsample(simple_tree(256),4))


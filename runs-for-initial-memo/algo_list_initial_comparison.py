from frb_olympics import add_algo, simple_direct, sloth, simple_tree, bonsai

add_algo(simple_direct(0.1))
add_algo(sloth(0.1,nupsample=2))
add_algo(simple_tree(1024))
add_algo(bonsai(1024))

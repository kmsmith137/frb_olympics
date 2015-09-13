from frb_olympics import add_algo, sloth, bonsai, downsample

add_algo(sloth(0.1,nupsample=2))
add_algo(bonsai(4096,nupsample=4))
add_algo(bonsai(2048,nupsample=2))
add_algo(bonsai(1024))
add_algo(downsample(bonsai(512),2))
add_algo(downsample(bonsai(512),4))

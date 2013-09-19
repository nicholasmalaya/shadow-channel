#!/bin/py
#
# run the channel and create plots comparable to keefe (92)
#
import sys
sys.path.append("../channel/")
import channel 

# try subprocess to 
# redirect all the screen outputs in C  to a buffer
#  which we can access in Python through channel.stdout
#
# http://docs.python.org/2/library/subprocess.html

# will segfault, at present
# channel.init(16,32,16,1.6,1.6,32,1,0.01,1,5,1069,0)
channel.init(16,32,16,1.6,1.6,32,1,0.01,1,5,100,0)

U = channel.getsoln(0)

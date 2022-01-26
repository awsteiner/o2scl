import numpy

def fun(x):
    ret=0
    for i in range(0,len(x)):
        ret=ret+x[i]*numpy.pi
    return ret

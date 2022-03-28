import numpy

def fun(x):
    y=[]
    for i in range(0,len(x)):
        y.append(x[i]*numpy.pi)
    return y

def fun2(x):
    y=[]
    for i in range(0,len(x)):
        y.append(x[i]*numpy.pi)
    return 1,y

from random import uniform
from math import atan

# monte carlo
# n - (optional) number of samples
def approxPi(N=4000000):
    insideCircle = 0
    R = 1 # const
    RSQ = 1**2
    for _ in range(N):
        xSq = uniform(0,R)**2
        ySq = uniform(0,R)**2
        if xSq+ySq < RSQ:
            insideCircle +=1
    return 4*(insideCircle/N)

# monte carlo integration
# f(a) must be a local min or max and f(b) must be a local min or max
# in range a -> b
# n - (optional) number of samples
def approxIntegral(a,b,fn, N=10000000):
    sign = 1
    lessThan = 0
    ymin = fn(a)
    ymax = fn(b)
    if ymin > ymax:
        ymin, ymax = ymax, ymin
    for _ in range(N):
        y = uniform(ymin,ymax)
        x = uniform(a,b)
        fx = fn(x)
        if fx > y:
            lessThan +=1
    constArea = ymin*(b-a)
    unkArea = (ymax-ymin)*(b-a)*(lessThan/N)
    return (constArea+unkArea)*sign
    
    

# Leibniz formula (4 * arctan(1) = pi)
# N - (optional) number of iterations
def approxPi2(N=100000000):
    quarterPi = 0
    add = True
    for i in range(1,N+1,2):
        if add:
            quarterPi += 1/i
        else:
            quarterPi -= 1/i
        add ^= True
    return 4*quarterPi

# simmilar to approxPi2 uses related Leibniz-Madhava series
# x - angle in radians
# N - (optional) number of iterations
def approxAtan(x, N=10000000):
    if x <= -1 or x >= 1:
        raise ValueError("only applicable for x: -1 to 1")
    add = True
    y = 0
    for i in range(1,N+1,2):
        if add:
            y += x**i/i
        else:
            y -= x**i/i
        add ^= True
    return y

def addAtan(x0,x1):
    if x0 == x0 == 1:
        raise ValueError("inputs cannot both equal 1")
    newX = (x0 + x1)/ (1-x0*x1)
    return atan(newX)

def demo():
    print(f"Monte Carlo approximation of π: {approxPi()}")
    print(f"Monte Carlo integration (4√(1-x²) approximation of π:", approxIntegral(0,1,lambda x: sqrt(1-x**2))*4)
    print(f"Leibniz approximation of π: {approxPi2()}")
    print(f"Leibniz-Madhava approximation of atan(.5): {approxAtan(.5)}")
    print(f"atan(.99) + atan(.99) = {addAtan(.99,.99)}")

demo()
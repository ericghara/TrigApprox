from random import uniform
from math import atan, sqrt

# monte carlo approximation of pi via euclidean geometry
#
# N - (optional) number of samples
def approxPi(N=4000000):
    insideCircle = 0
    R = 1 # const
    RSQ = R**2
    for _ in range(N):
        xSq = uniform(0,R)**2
        ySq = uniform(0,R)**2
        if xSq+ySq < RSQ:
            insideCircle +=1
    return 4*(insideCircle/N)

# monte carlo integration
# f(a) must be a local min or max and f(b) must be a local min or max
# in range from a to b
#
# a - start of integration range
# b - end of integration range
# fn - the function to integrate. It must be a single variable function
# N - (optional) number of samples
def approxIntegral(a,b,fn, N=10000000):
    lessThanfx = 0
    ymin = fn(a)
    ymax = fn(b)
    if ymin > ymax:
        ymin, ymax = ymax, ymin
    for _ in range(N):
        y = uniform(ymin,ymax)
        x = uniform(a,b)
        if fn(x) > y:
            lessThan +=1
    p = lessThanfx / N # probability distribution
    constArea = ymin*(b-a) # region less than f(x) over a, b with 100% probability
    dynamicArea = (ymax-ymin)*(b-a)*p # area of unknown region calculated from p
    return constArea+dynamicArea
    
# Leibniz formula approximates 4 * arctan(1),
# which ≡ pi
#
# N - (optional) number of iterations
def approxPi2(N=100000000):
    quarterPi = 0
    add = True
    for i in range(1,2*N+1,2):
        if add:
            quarterPi += 1/i
        else:
            quarterPi -= 1/i
        add ^= True
    return 4*quarterPi

# simmilar to approxPi2 uses related Leibniz-Madhava series. Only
# applicable for calculating arctan of angles -1 to 1 (radians).
# 
# x - angle in radians
# N - (optional) number of iterations
def approxAtan(x, N=10000000):
    if x <= -1 or x >= 1:
        raise ValueError("only applicable for x: -1 to 1")
    add = True
    y = 0
    for i in range(1,2*N+1,2):
        if add:
            y += x**i/i
        else:
            y -= x**i/i
        add ^= True
    return y

# Computes the sum of the arctan of two angles.  The angles must not both
# equal 1.  Function should provide a more accurate result than calculating
# both arctans separately and adding.
#
# x0 - angle in radians
# x1 - angle in radians
# result - atan(x0) + atan(x1)
def addAtan(x0,x1):
    if x0 == x0 == 1:
        raise ValueError("inputs cannot both equal 1")
    newX = (x0 + x1)/ (1-x0*x1)
    return atan(newX)

def demo():
    print("***-->>> Function Demo <<<---***")
    print(f"Monte Carlo approximation of π: {approxPi()}")
    print(f"Monte Carlo integration (4√(1-x²) approximation of π:", approxIntegral(0,1,lambda x: sqrt(1-x**2))*4)
    print(f"Leibniz approximation of π: {approxPi2()}")
    print(f"Leibniz-Madhava approximation of atan(.5): {approxAtan(.5)}")
    print(f"atan(.99) + atan(.99) = {addAtan(.99,.99)}")

demo()

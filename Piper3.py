import matplotlib.pyplot as plt
import math
import time


# Parametre
A = 976/2980 # R1 / R2
B = 0.000005 # Reverse current
V_0 = 1.1*0.02585 # Ideality factor, n, multiplied by thermal voltage, V_t, n is typically between 1 and 2
R_0 = 98.6
L = 0.243
C = 0.00000001

n = 10000
start = 0
end = 0.1
h = (end-start)/n
res = 10000

# Prøver å havne i grensesyklus så fort som mulig for å kunne telle frekvens mer nøyaktig
x0 = 0.4
y0 = 0

t = [start]
x = [x0]
y = [y0]


def gInverse(y):
    return A*y + B*math.sinh(y/(V_0))


def makeTable(f, start, end):
    variable = []
    function = []
    h = (end-start)/res

    i = start
    while i < end:
        variable.append(i)
        function.append(f(i))
        i += h
    
    return variable, function

"""
# Old g, bad runtime
def g(val, var, func):
    minDiff = float("inf")
    index = len(func)-1

    for i in range(len(func)):
        if abs(val - func[i]) < minDiff:
            minDiff = abs(val-func[i])
            index = i
        
    return var[index]
"""

def g(val, var, func):
    minDiff = float("inf")
    index = int(val*res)
    run = True
    while run:
        if abs(val - func[index+1]) < minDiff:
            minDiff = abs(val - func[index+1])
            index += 1
        elif abs(val - func[index-1]) < minDiff:
            minDiff = abs(val - func[index-1])
            index -= 1
        else:
            run = False
            break
    return var[index]


def gDer(y, var, func):
    return 1/(A + B*math.cosh(g(y, var, func)/(V_0))/(V_0))


def nextX(x, y):
    return x + h*y


def nextY(x, y, var, func):
    return y - h*((R_0/L)*(1-gDer(x, var, func))*y+(1/(L*C))*x)


var, gInv = makeTable(gInverse, -2, 2)


def countFreq(t, x):
    prevX = x[0]
    foundStart = False
    for i in range(1, len(t)):
        if prevX < 0 and x[i] > 0 and not foundStart:
            start = t[i]
            prevX = x[i]
            foundStart = True
        elif prevX > 0 and x[i] < 0 and foundStart:
            end = t[i]
            break
    return 1/(2*(end-start))

begin = time.time()
for i in range(n):
    t.append(start + h*i)
    x0 = nextX(x0, y0)
    y0 = nextY(x0, y0, var, gInv)
    x.append(x0)
    y.append(y0)
end = time.time()
print(end-begin)


#print(countFreq(t[-1000:], x[-1000:]))
#plt.plot(x, y)
plt.plot(t[-1000:], x[-1000:])
plt.show()
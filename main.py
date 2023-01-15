import numpy
import matplotlib.pyplot as plt
import scipy

n = int(input())
length = 2
h = length/n
# i - indeks obliczanego przez nas e, x - wartosc jaka przyjmuje nasze e
def countE(i,x):
    if(x < 0 or x > 2 ): return 0
    start = h * (i - 1)
    middle = h * i
    end = h * (i + 1)
    if(start <= x <= middle): return (x - start) / h
    if(middle < x <= end): return (end - x) / h
    return 0

def derivativeE(i,x):
    start = h * (i - 1)
    middle = h * i
    end = h * (i + 1)
    if(start <= x <= middle): return 1/h
    if(middle < x <= end): return (-1)/h
    return 0

def countB(i,j):
    result = (-2) * countE(i,0) * countE(j,0) # obliczenie -2ei(0)ej(0)

    def countIntegrate(x):
        if(x <= 1):E = 2
        else: E = 6
        return E * derivativeE(i,x) * derivativeE(j,x)

    result += scipy.integrate.quadrature(countIntegrate,max((i-1) * h,(j-1) * h,0),min((i+1) * h,(j+1) * h,2),vec_func=False)[0]
    return result


bArray = [[0for i in range(n)]for i in range(n)]
for i in range(n):
    for j in range(n):
        if(abs(i-j) <=1):
            bArray[i][j] = countB(j,i)

lArray = [0 for i in range(n)]
for i in range(n):
    lArray[i] = (-20) * countE(i,0) + 2 * 3 * countE(i,0) #pomijamy całkę, ponieważ wiemy, że pochodna u z daszkiem jest rowna 0

uArray = numpy.linalg.inv(bArray).dot(lArray)
xArray = [i*h for i in range(n)]
for i in range(n):
    sumY = 0
    for j in range(n):
        sumY += uArray[i] * countE(j,xArray[i])
    sumY += 3
    uArray[i] = sumY
uArray = numpy.append(uArray,3)
xArray.append(2)
plt.plot(xArray,uArray)
plt.show()
#print(uArray[0])

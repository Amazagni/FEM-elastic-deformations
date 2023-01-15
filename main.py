import numpy
import matplotlib.pyplot as plt
#n = 3

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
    numberOfNodes = 80
    nodes, weights = numpy.polynomial.legendre.leggauss(numberOfNodes)
    for nId in range(numberOfNodes):nodes[nId] += 1 # przejście z (-1,1) na (0,2)
    # nodesI = nodes.copy()
    # nodesJ = nodes.copy()
    # h2 = h/2
    # for nId in range(numberOfNodes):
    #     nodesI[nId] = h2 * i + h2 * nodes[nId]
    #     nodesJ[nId] = h2 * j + h2 * nodes[nId]
    # print(i)
    # print(nodesI)
    for nId in range(numberOfNodes):
        if(nodes[nId] <= 1):E = 2
        else: E = 6
    #      result += weights[nId] * E * derivativeE(i,nodesI[nId]) * derivativeE(j,nodesJ[nId])
        result += weights[nId] * E * derivativeE(i,nodes[nId]) * derivativeE(j,nodes[nId])
    return result


bArray = [[0for i in range(n)]for i in range(n)]
for i in range(n):
    for j in range(n):
        if(abs(i-j) <=1):
            bArray[i][j] = countB(i,j)

lArray = [0 for i in range(n)]
for i in range(n):
    lArray[i] = (-20) * countE(i,0) + 2 * 3 * countE(i,0) #pomijamy całkę, ponieważ wiemy, że pochodna u z daszkiem jest rowna 0
#print(bArray)
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
#print(xArray)
#print(uArray)
plt.plot(xArray,uArray)
plt.show()
print(uArray[0])

#print(countE(0,0.5))

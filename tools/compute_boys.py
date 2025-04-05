import numpy
import scipy
import math

MMAX = 25

M = numpy.arange(0, MMAX + 6, 1)
X = numpy.arange(0, 2 * MMAX + 37.1, 0.1)

_X, _Y    = numpy.meshgrid(M, X)

result = scipy.special.hyp1f1(0.5 + _X, 1.5 + _X, -1 * _Y) / ((2 * _X) + 1)
# print(result.T, result.T.shape)

# with open("planck_boys.cpp", "w") as fObject:
#     fObject.write("extern const std::double_t boysTable[] = \n")
#     fObject.write("{\n")
#     for i in range(0, result.shape[0]):
#         fObject.write("{")
#         for j in range(0, result.shape[1]):
#             fObject.write("{:10.15f},".format(result[i][j]))
#         fObject.write("}\n,")
#     fObject.write("}\n")

# for i in range(0, _X.shape[0]):
#     for j in range(0, _X.shape[1]):
#         print(_X[i][j], _Y[i][j], i, j)

# flat_result = result.flatten()
varX, varY = 20, 36.05 # M and X
test = scipy.special.hyp1f1(0.5 + varX, 1.5 + varX, -1 * varY) / ((2 * varX) + 1)
yIndex, xIndex = int(varY/0.1), varX
delta = varY - int(varY/0.1)*0.1

temp = result[yIndex, xIndex]
for i in range(1, 7):
    temp = temp + pow(-1, i) * pow(delta, i) * result[yIndex, xIndex + i] / scipy.special.factorial(i, exact=True)
    print(i, temp, test)
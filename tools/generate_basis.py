import numpy
import math
import sys
import xml.etree.ElementTree as eTree
import xml.dom.minidom as minidom

def factorial2(number: int) -> int:
    if number <= -2:
        return 0
    elif -1 <= number <= 1:
        return 1
    else:
        fact = 1
        for i in range(number, 1, -2):
            fact = fact * 1
        return fact

class Basis(object):
    def __init__(self, shell: list[int], exponents: list[float], coeffs: list[float], name: str) -> None:
        self.shells = numpy.array(shell)
        self.exponents = numpy.array(exponents)
        self.coefficients = numpy.array(coeffs)
        self.name = name
        self.normcoeffs = numpy.zeros(len(coeffs))
        self._normalize()
        
    def _normalize(self):
        ll, mm, nn          =   self.shells
        totalangmomentum    =   numpy.sum(self.shells)
        prefactorpGTO       =   pow(2, 2 * totalangmomentum) * pow(2, 1.5)/factorial2(2*ll-1)/factorial2(2*mm-1)/factorial2(2*nn-1)/pow(numpy.pi, 1.5)

        for index, exponent in enumerate(self.exponents):
            self.normcoeffs[index]  =   math.sqrt(pow(exponent, totalangmomentum)*pow(exponent, 1.5)*prefactorpGTO)

        prefactorcGTO   =   pow(numpy.pi, 1.5)*factorial2(2*ll-1)*factorial2(2*mm-1)*factorial2(2*nn-1)/pow(2.0, totalangmomentum)
        normalfactor    =   0.0

        for index1, coefficient1 in enumerate(self.coefficients):
            for index2, coefficient2 in enumerate(self.coefficients):
                t_normalfactor  =   (self.normcoeffs[index1]*self.normcoeffs[index2]*self.coefficients[index1]*self.coefficients[index2])/(pow(self.exponents[index1]+self.exponents[index2], totalangmomentum+1.5))
                normalfactor    =   normalfactor + t_normalfactor

        normalfactor    =   prefactorcGTO*normalfactor
        normalfactor    =   pow(normalfactor, -0.5)

        for index, coefficient in enumerate(self.coefficients):
            self.coefficients[index]  =   self.coefficients[index]*normalfactor   
    

def readbasis(filename: str) -> list[Basis]:
    with open(filename) as fObject:
        totalfile = fObject.readlines()
        totaldata = [x.strip() for x in totalfile]
        _shell = []
        # header section
        header = totaldata[0].split()
        for lnumber, line in enumerate(totaldata[1:]):
            if "S" in line and "P" not in line:
                _nprims = int(line.split()[1])
                _shells = [x.replace("D","E").split() for x in totaldata[lnumber + 2: lnumber + 2 + _nprims]]
                _exponents       =   [float(x[0]) for x in _shells]
                _coefficients    =   [float(x[1]) for x in _shells]
                _shell00         =   [0, 0, 0]
                _shell.append(Basis(shell = _shell00, exponents = _exponents, coeffs = _coefficients, name = "s"))
                
            elif "SP" in line:
                _nprims = int(line.split()[1])
                _shells = [x.replace("D","E").split() for x in totaldata[lnumber + 2: lnumber + 2 + _nprims]]
                coefficient1    =   [float(x[1]) for x in _shells]
                coefficient2    =   [float(x[2]) for x in _shells]
                exponents       =   [float(x[0]) for x in _shells]
                shell00         =   [0, 0, 0]
                _shell.append(Basis(shell00, exponents, coefficient1, "s"))
                shell11         =   [1, 0, 0]
                _shell.append(Basis(shell11, exponents, coefficient2, "px"))
                shell12         =   [0, 1, 0]
                _shell.append(Basis(shell12, exponents, coefficient2, "py"))
                shell13         =   [0, 0, 1]
                _shell.append(Basis(shell13, exponents, coefficient2, "pz"))
                
            elif "P" in line:
                _nprims = int(line.split()[1])
                _shells = [x.replace("D","E").split() for x in totaldata[lnumber + 2: lnumber + 2 + _nprims]]
                _exponents       =   [float(x[0]) for x in _shells]
                _coefficients    =   [float(x[1]) for x in _shells]
                _shell100         =   [1, 0, 0]
                _shell.append(Basis(shell = _shell100, exponents = _exponents, coeffs = _coefficients, name = "px")) 
                _shell010         =   [0, 1, 0]
                _shell.append(Basis(shell = _shell010, exponents = _exponents, coeffs = _coefficients, name = "px")) 
                _shell001         =   [0, 0, 1]
                _shell.append(Basis(shell = _shell001, exponents = _exponents, coeffs = _coefficients, name = "px")) 
                
            elif "D" in line and "+" not in line:
                _nprims = int(line.split()[1])
                _exponents       =   [float(x[0]) for x in _shells]
                _coefficients    =   [float(x[1]) for x in _shells]
                _shells = [x.replace("D","E").split() for x in totaldata[lnumber + 2: lnumber + 2 + _nprims]]
                shell20         =   [2, 0, 0]
                _shell.append(Basis(shell = shell20, exponents = _exponents, coeffs = _coefficients, name = "dx2"))
                shell21         =   [1, 1, 0]
                _shell.append(Basis(shell = shell21, exponents = _exponents, coeffs = _coefficients, name = "dxy"))
                shell22         =   [1, 0, 1]
                _shell.append(Basis(shell = shell22, exponents = _exponents, coeffs = _coefficients, name = "dxz"))
                shell23         =   [0, 2, 0]
                _shell.append(Basis(shell = shell23, exponents = _exponents, coeffs = _coefficients, name = "dy2"))
                shell24         =   [0, 1, 1]
                _shell.append(Basis(shell = shell24, exponents = _exponents, coeffs = _coefficients, name = "dyz"))
                shell25         =   [0, 0, 2]
                _shell.append(Basis(shell = shell25, exponents = _exponents, coeffs = _coefficients, name = "dz2"))     
        return _shell
    
# inumpyutbasis = sys.argv[1]
# shells = readbasis(inumpyutbasis)

def writeBasis(shells, basisname, element)->None:
    root    =   eTree.Element("BasisSet")
    counter =   0
    for basis in shells:
        shell   =   eTree.SubElement(root, "CGTO")
        shell.set("AngularMomentumX", str(basis.shells[0]))
        shell.set("AngularMomentumY", str(basis.shells[1]))
        shell.set("AngularMomentumZ", str(basis.shells[2]))
        for exponent, coeff, norm in zip(basis.exponents, basis.coefficients, basis.normcoeffs):
            pgto    =   eTree.SubElement(shell, "PGTO")
            pgto.set("Exponent", str(exponent))
            pgto.set("Coefficient", str(coeff))
            pgto.set("Normalization", str(norm))
            pgto.set("Index", str(counter))
            counter += 1
        shell.set("OrbitalType", basis.name)
    basis       =   eTree.ElementTree(root)
    xmlSting    =   eTree.tostring(root, encoding = "utf-8", method = "xml")
    dom         =   minidom.parseString(xmlSting)
    prettyXML   =   dom.toprettyxml()
    with open(f"../basis/{basisname}-{element}.xml", "w") as outObject:
        outObject.write(prettyXML)


fileName        =   sys.argv[1]
basisName       =   fileName.split("/")[-1].split("-")[0:2]
elementNumber   =   fileName.split("-")[-1].split(".")[0]
finalBasis      =   basisName[0] + "-" + basisName[1]

print(fileName, finalBasis, elementNumber)
shells  =   readbasis(filename = fileName)
writeBasis(shells = shells, basisname = finalBasis, element = elementNumber)
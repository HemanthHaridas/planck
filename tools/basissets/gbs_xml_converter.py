import numpy as np
import scipy as scp
import math as mt
import sys as sys
import typing as tp
import xml.etree.ElementTree as eTree
import xml.dom.minidom as minidom

class Basis(object):
    def __init__(self, shell, exponents, coefficients, name):
        self.shell          =   np.array(shell)
        self.exponents      =   np.array(exponents)
        self.coefficients   =   np.array(coefficients)
        self.normcoeffs     =   np.zeros(self.coefficients.size)
        self.name           =   name

        self.normalizepGTO()

    def normalizepGTO(self)->None:
        ll, mm, nn          =   self.shell
        totalangmomentum    =   np.sum(self.shell)
        prefactorpGTO       =   pow(2, 2*totalangmomentum)*pow(2, 1.5)/scp.special.factorial2(2*ll-1)/scp.special.factorial2(2*mm-1)/scp.special.factorial2(2*nn-1)/pow(np.pi, 1.5)

        for index, exponent in enumerate(self.exponents):
            self.normcoeffs[index]  =   mt.sqrt(pow(exponent, totalangmomentum)*pow(exponent, 1.5)*prefactorpGTO)

        prefactorcGTO   =   pow(np.pi, 1.5)*scp.special.factorial2(2*ll-1)*scp.special.factorial2(2*mm-1)*scp.special.factorial2(2*nn-1)/pow(2.0, totalangmomentum)
        normalfactor    =   0.0

        for index1, coefficient1 in enumerate(self.coefficients):
            for index2, coefficient2 in enumerate(self.coefficients):
                t_normalfactor  =   (self.normcoeffs[index1]*self.normcoeffs[index2]*self.coefficients[index1]*self.coefficients[index2])/(pow(self.exponents[index1]+self.exponents[index2], totalangmomentum+1.5))
                normalfactor    =   normalfactor + t_normalfactor

        normalfactor    =   prefactorcGTO*normalfactor
        normalfactor    =   pow(normalfactor, -0.5)

        for index, coefficient in enumerate(self.coefficients):
            self.coefficients[index]  =   self.coefficients[index]*normalfactor
    
def readBasis(filename: str) -> tp.List[Basis]:
    with open(filename) as target:
        basisdata   =   target.readlines()
        headerline  =   basisdata[0].split()[0]
        _atomname   =   headerline
        shells      =   []
        for lnumber, line in enumerate(basisdata[1:]):
            if "S" in line and "P" not in line:
                nprims          =   int(line.split()[1])
                pgtodata        =   [x.replace('D', 'E').split() for x in basisdata[lnumber+2:lnumber+2+nprims]]
                exponents       =   [float(x[0]) for x in pgtodata]
                coefficients    =   [float(x[1]) for x in pgtodata]
                shell00         =   [0, 0, 0]
                shells.append(Basis(shell00, exponents, coefficients, "s"))

            if "P" in line:
                nprims          =   int(line.split()[1])
                pgtodata        =   [x.replace('D', 'E').split() for x in basisdata[lnumber+2:lnumber+2+nprims]]
                coefficient1    =   [float(x[1]) for x in pgtodata]
                coefficient2    =   [float(x[2]) for x in pgtodata]
                exponents       =   [float(x[0]) for x in pgtodata]
                shell00         =   [0, 0, 0]
                shells.append(Basis(shell00, exponents, coefficient1, "s"))
                shell11         =   [1, 0, 0]
                shells.append(Basis(shell11, exponents, coefficient2, "px"))
                shell12         =   [0, 1, 0]
                shells.append(Basis(shell12, exponents, coefficient2, "py"))
                shell13         =   [0, 0, 1]
                shells.append(Basis(shell13, exponents, coefficient2, "pz"))

            if "D" in line and "+" not in line:
                nprims          =   int(line.split()[1])
                pgtodata        =   [x.replace('D', 'E').split() for x in basisdata[lnumber+2:lnumber+2+nprims]]
                exponents       =   [float(x[0]) for x in pgtodata]
                coefficients    =   [float(x[1]) for x in pgtodata]
                shell20         =   [2, 0, 0]
                shells.append(Basis(shell20, exponents, coefficients, "dx2"))
                shell21         =   [1, 1, 0]
                shells.append(Basis(shell21, exponents, coefficients, "dxy"))
                shell22         =   [1, 0, 1]
                shells.append(Basis(shell22, exponents, coefficients, "dxz"))
                shell23         =   [0, 2, 0]
                shells.append(Basis(shell23, exponents, coefficients, "dy2"))
                shell24         =   [0, 1, 1]
                shells.append(Basis(shell24, exponents, coefficients, "dyz"))
                shell25         =   [0, 0, 2]
                shells.append(Basis(shell25, exponents, coefficients, "dz2"))

    return shells

basisFile   =   sys.argv[1]
basisSet    =   readBasis(basisFile)

def writeBasis(shells)->None:
    root    =   eTree.Element("BasisSet")
    for basis in basisSet:
        shell           =   eTree.SubElement(root, "Shell")
        shell.set("Shell", basis.shell)
        shell.set("Exponents", basis.exponents)
        shell.set("Coefficients", basis.coefficients)
        shell.set("Normalizations", basis.normcoeffs)
        shell.set("ShellType", basis.name)
    basis       =   eTree.ElementTree(root)
    xmlSting    =   eTree.tostring(root, encoding = "utf-8", method = "xml")
    dom         =   minidom.parseString(xmlSting)
    prettyXML   =   dom.toprettyxml()
    with open(f"../basis/{basisFile[:-6]}.xml", "w") as outObject:
        outObject.write(prettyXML)

writeBasis(basisSet)
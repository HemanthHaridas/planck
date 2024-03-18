import numpy as np
import scipy as scp
import math as mt
import sys as sys
import typing as tp
import xml.etree.ElementTree as eTree
import xml.dom.minidom as minidom

"""
 Planck
 Copyright (C) 2024 Hemanth Haridas, University of Utah
 Contact: hemanthhari23@gmail.com
 
 This program is free software: you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation, either version 3 of the License, or a later version.
 
 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along with
 this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
    print(filename)
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

def writeBasis(shells, basisname, element)->None:
    root    =   eTree.Element("BasisSet")
    counter =   0
    for basis in shells:
        shell   =   eTree.SubElement(root, "CGTO")
        shell.set("AngularMomentumX", str(basis.shell[0]))
        shell.set("AngularMomentumY", str(basis.shell[1]))
        shell.set("AngularMomentumZ", str(basis.shell[2]))
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
shells  =   readBasis(filename = fileName)
writeBasis(shells = shells, basisname = finalBasis, element = elementNumber)
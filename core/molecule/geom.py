import  numpy
import  typing  
import  scipy.linalg
import  math
import  scipy
import  os
from    . import periodic_table_helper

class Molecule:
  def __init__(self, atoms : typing.List[typing.Tuple[str, typing.Tuple[float, float, float]]], charge : int, multiplicity : int, basis: str, symmetry = True) -> None:
    self.n_atoms      = len(atoms)
    self.geometry     = atoms
    self.charge       = charge
    self.muliplicity  = multiplicity
    self.basis_set    = basis
    self.exponents    = []
    self.centers      = []
    self.coefficients = []
    self.normalize    = []
    self.shells       = []
    self.charges      = []
    self.point_group  = None
    self.reorient     = None
    
  def read_basis(self):
    self.charges      = [periodic_table_helper.get_element(row[0]) for row in self.geometry]
    object_path       = os.path.dirname(os.path.dirname(__file__))
    for atomindex, atom in enumerate(self.charges):
      basis_file      = "/basissets/{}/{}-{}.txt".format(self.basis_set, self.basis_set, atom)
      basis_path      = object_path+basis_file
      print(basis_path)
      with open(basis_path) as basis_object:
        basis_data  = basis_object.readlines()
        for lnumber, line in enumerate(basis_data[1:]):
          if "S" in line and "P" not in line:
              nprims          =   int(line.split()[1])
              pgtodata        =   [x.replace('D', 'E').split() for x in basis_data[lnumber+2:lnumber+2+nprims]]
              exponents       =   [float(x[0]) for x in pgtodata]
              coefficients    =   [float(x[1]) for x in pgtodata]
              shell00         =   [0, 0, 0]

              self.shells.append(shell00)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

          if "P" in line:
              nprims          =   int(line.split()[1])
              pgtodata        =   [x.replace('D', 'E').split() for x in basis_data[lnumber+2:lnumber+2+nprims]]
              coefficients    =    [float(x[1]) for x in pgtodata]
              exponents       =   [float(x[0]) for x in pgtodata]
              shell11         =   [1, 0, 0]
              shell12         =   [0, 1, 0]
              shell13         =   [0, 0, 1]
            
              self.shells.append(shell11)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

              self.shells.append(shell12)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

              self.shells.append(shell13)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

          if "D" in line and "+" not in line:
              nprims          =   int(line.split()[1])
              pgtodata        =   [x.replace('D', 'E').split() for x in basis_data[lnumber+2:lnumber+2+nprims]]
              exponents       =   [float(x[0]) for x in pgtodata]
              coefficients    =   [float(x[1]) for x in pgtodata]
              shell20         =   [2, 0, 0]
              shell21         =   [1, 1, 0]
              shell22         =   [1, 0, 1]
              shell23         =   [0, 2, 0]
              shell24         =   [0, 1, 1]
              shell25         =   [0, 0, 2]

              self.shells.append(shell20)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

              self.shells.append(shell21)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

              self.shells.append(shell22)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

              self.shells.append(shell23)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

              self.shells.append(shell24)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

              self.shells.append(shell25)
              self.exponents.append(exponents)
              self.coefficients.append(coefficients)
              self.centers.append(numpy.array(self.geometry[atomindex][1])*1.8897259886)

  def normalize_basis(self):
     for basis_object in zip(self.shells, self.coefficients, self.exponents):
        print(basis_object)
        total_momentum  = sum(basis_object[0])
        prefactor_pgto  = pow(2, 2*total_momentum)*pow(2, 1.5)/scipy.special.factorial2(2*basis_object[0][0]-1)/scipy.special.factorial2(2*basis_object[0][1]-1)/scipy.special.factorial2(2*basis_object[0][2]-1)/pow(numpy.pi, 1.5)
        self.normalize.append([math.sqrt(pow(exponent, total_momentum)*pow(exponent, 1.5)*prefactor_pgto) for exponent in basis_object[2]])
        
import sys
from planck.core.molecule.geom import Molecule

class RHF:
    def __init__(self, molecule_object: Molecule, scf_cycles: int = 100, diis: bool = True) -> None:
        self.molecule   = molecule_object
        self.scf_cycles = scf_cycles
        self.diis       = diis
        self.n_alpha    = 0
        self.n_beta     = 0
    
        if self.check_charge_multiplicity():
            pass
        else:
            sys.exit(-1)

    def check_charge_multiplicity(self) -> bool:
        total_electrons = sum(self.molecule.charges) + self.molecule.charge
        try:
            assert total_electrons%2 == 0, f"[ERROR]    Total Number of {total_electrons} electrons and Multiplicity of 0 do not match" 
        except AssertionError as Message:
            print(Message)
            return False

        self.n_alpha    =   total_electrons // 2
        self.n_beta     =   total_electrons - self.n_alpha

        try:
            assert self.n_alpha == self.n_beta, f"[ERROR]   Number of alpha electrons {self.n_alpha} is not equal to Number of beta electrons {self.n_beta}"
        except AssertionError as Message:
            print(Message)
            return False
        return True
    
import typing
import periodictable

def get_element(atomname: str) -> int:
    return periodictable.elements.symbol(atomname).number

def get_mass(atomname: str) -> float:
    return periodictable.elements.symbol(atomname).mass


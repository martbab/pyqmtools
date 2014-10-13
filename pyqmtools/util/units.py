#!/usr/bin/env python
# collection of helper functions for conversion of various units
# used in QM software
# according to CODATA 2010
import numpy as np

# definition of Avogadro's constant
AVOGADRO = 6.02214129e23

# definition of 1 Hartree in terms of J
HARTREE2J = 4.35974434e-18

# conversion of 1 Hartree to kJ
HARTREE2KJ = HARTREE2J / 1.00e3

# conversion of Hartree to kJ/mole
HARTREE2KJMOL = HARTREE2KJ * AVOGADRO

# conversion from kJ/mole to kcal/mole
KCALMOL2KJMOL = 4.184

# conversion from Hartrees to kcal/mole
HARTREE2KCALMOL = HARTREE2KJMOL / KCALMOL2KJMOL


# conversion of 1 Bohr to Angstrom
BOHR2ANGSTROM = 0.52917721092


def _convert(arg, units, reverse = False):
    factor = units
    result = None
    if reverse:
        factor = 1.0 / units
    

    try:
        if type(arg) == str:
            raise TypeError
        elif type(arg) == list:
            result = []

            for a in arg:
                if type(arg) == str:
                    raise TypeError
                
                result.append(a * factor)
        else:
            result = arg * factor

        return result
    except TypeError, ValueError:
        raise TypeError("Cannot convert units for type: \"\"" % type(arg))


# convenience functions
def bohr2angstrom(arg, reverse = False):
    return _convert(arg, units = BOHR2ANGSTROM, reverse = reverse)

def hartree2kj(arg, reverse = False):
    """
    converts data from hartrees to kilojoules and reverse
    """
    return _convert(arg, units = HARTREE2KJ, reverse = reverse)
           

def hartree2kjmol(arg, reverse = False):
    """
    converts data from hartrees to kilojoules per mole and reverse
    """
    return _convert(arg, units = HARTREE2KJMOL, reverse = reverse)


def hartree2kcalmol(arg, reverse = False):
    """
    converts data from hartrees to kcal per mole and reverse
    """
    return _convert(arg, units = HARTREE2KCALMOL, reverse = reverse)


def gradqm2mm(arg, reverse = False):
    """
    converts units of cartesian gradient used in QM software (usually 
    hartree/bohr, check software documentation) to units used in MM 
    forcefields (usually kilocalorie/mol/Angstrom)
    """
    partial = _convert(arg, units = HARTREE2KCALMOL, reverse = reverse)
    full_result = _convert(
        partial, 
        units = 1.0 / BOHR2ANGSTROM, 
        reverse = reverse
    )
    return full_result


def hessqm2mm(arg, reverse = False):
    """
    converts units of cartesian force constants used in QM software (usually
    hartree / bohr**2) to units used in MM forcefields (usually 
    kcal/mol/Angstrom**2)
    """
    partial = _convert(arg, units = HARTREE2KCALMOL, reverse = reverse)
    full_result = _convert(
        partial, 
        units = 1.0 / (BOHR2ANGSTROM**2), 
        reverse = reverse
    )
    return full_result

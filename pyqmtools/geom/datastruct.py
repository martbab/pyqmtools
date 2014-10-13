#!/usr/bin/env python
#TODO: implement a detection of close contacts between voxels and molecule
"""
This file contains definitions of these classes:


(c) 2012 Martin Babinsky
"""

import math
from pyqmtools.util import units as u
from pyqmtools.util import elements as e

_BOHR_NAME = 'bohr'
_ANG_NAME = 'ang'

class Atom(object):
    def __init__(self,
        element = "",
        coords = [0.0, 0.0, 0.0],
        units = _ANG_NAME
    ):
        self.units = units
        self.coords = coords
        self.element = element
        self.charge = None

    @property
    def element(self):
        return self.__element
        
    @element.setter
    def element(self, el):
        self.__element = el

    @property
    def units(self):
        return self.__units
    
    @units.setter
    def units(self, units):
        if units.lower() == self.__units:
            pass
        elif units.lower() == _ANG_NAME:
            self.__units = _ANG_NAME
            self.coords = u.bohr2angstrom(
                self.coords
            )

        elif units.lower() == _BOHR_NAME:
            self.__units = _BOHR_NAME
            self.coords = u.bohr2angstrom(
                self.coords,
                reverse = True
            )
        else:
            raise ValueError(
                "Units must be specified as either '%s' or '%s'!" %\
                    (_ANG_NAME, _BOHR_NAME)
            )

    def read_xyz_record(self, record = ""):
        self.units = _ANG_NAME
        (at_name, x, y, z) = record.strip().split()

        self.name = at_name
        self.coords = map(float, [x, y, z])
        
    def read_cube_record(self, record = ""):
        self.units = _BOHR_NAME
        (at_num, ch, x, y, z) = record.strip().split()

        at_n = int(at_num)

        for k in _ATOM_LOOKUP:
            if _ATOM_LOOKUP[k] == at_n:
                self.name = k

        self.coords = map(float, [x, y, z])

    def read_turbomole_record(self, record = ""):
        self.units = _BOHR_NAME
        (x, y, z, at_name) = record.strip().split()[:4]
        self.name = at_name.capitalize()
        self.coords = map(float, [x, y, z])

    def string_repr(self, format_ = 'xyz'):
        if format_ == 'xyz':
            return self._string_as_xyz()
        elif format_ == 'cube':
            return self._string_as_cube_coords()
        else:
            pass
        
    def _string_as_xyz(self, formatting = "%3s%15.7f%15.7f%15.7f\n"):
        self.units = _ANG_NAME
        return formatting % (self.name, self.coords[0], self.coords[1], \
            self.coords[2])

    def _string_as_cube_coords(self, 
        formatting = "%5d%12.6f%12.6f%12.6f%12.6f\n"
    ):
        self.units = _BOHR_NAME
        return formatting % (self.at_num, self.charge, self.coords[0], \
            self.coords[1], self.coords[2])
        

        

class Coordinates(object):
    def __init__(self, 
        units = _ANG_NAME
    ):
        self.units = units
        self.atoms = []
    
    @property
    def units(self):
        return self.__units
    
    @units.setter
    def units(self, units):
        if units.lower() == self.__units:
            pass
        elif units.lower() == _ANG_NAME:
            self.__units = _ANG_NAME

        elif units.lower() == _BOHR_NAME:
            self.__units = _BOHR_NAME
        else:
            raise ValueError(
                "Units must be specified as either '%s' or '%s'!" %\
                    (_ANG_NAME, _BOHR_NAME)
            )

        for a in self.atoms:
            a.units = units.lower()

    def get_units(self):
        return self.__units
        
    def add_atom(self, atom = None):
        if isinstance(atom, Atom):
            self.atoms.append(atom)
    
    def count_atoms(self):
        return len(self.atoms)

    def dump_as_xyz(self, outp_file = None, header = True, comment = ""):
        if header:
            outp_file.write("%d\n" % self.count_atoms())
            outp_file.write(comment)

        for atom in self.atoms:
            outp_file.write(atom.string_repr())

    def dump_as_cube(self, outp_file = None):
        for atom in self.atoms:
            outp_file.write(atom.string_repr(format_ = 'cube'))

    def read_from_xyz(self, inp_file, header = True):
        if header:
            n_atoms = int(inp_file.readline().strip())
            comment = inp_file.readline()

        for line in inp_file:
            helper_atom = Atom()
            helper_atom.read_xyz_record(line)
            self.add_atom(
                helper_atom
            )

    def read_from_cube(self, inp_file, header = True, units = _BOHR_NAME):
        header_count = 6
        lc = 0

        self.units = units

        if header:
            while lc < header_count:
                inp_file.readline()
                lc += 1

        for line in inp_file:
            line_record = line.strip().split()

            if len(line_record) > 5:
                break

            else:
                helper_atom = Atom(units = self.units)
                helper_atom.read_cube_record(line)
                self.add_atom(
                    helper_atom
                )

    def read_from_turbomole(self, inp_file):
        blk_begin = "$coord"
        blk_end = "$"
        blk_found = False

        for line in inp_file:
            if blk_begin in line:
                blk_found = True
                continue

            if blk_found:
                if blk_end in line:
                    break
                else:
                    helper_atom = Atom()
                    helper_atom.read_turbomole_record(line)
                    self.add_atom(
                        helper_atom
                    )



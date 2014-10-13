#!/usr/bin/env/python

"""
Set of classes to be used for the reading and writing molecular geometry between
various sources, such as XYZ files, Gaussian Cubes, and output files generated
by various QM software.
"""
from pyqmtools.util import units as u
from pyqmtools.util import elements as e
import datastruct as d

class GeomReadException(Exception):
    pass
    
    
class XYZIO(object):
    _atom_line_fmt = """{a.element:3s} {a.coords[0]:10.5f} {a.coords[1]:10.5f} 
{a.coords[2]:10.5f}\n"""
    def __init__(self,
    ):
        self.per_table = e.PeriodicTable()
        
    def read_atom(self, line, units = 'ang'):
        fields = line.strip().split()
        atom = d.Atom(
            units = units
        )
        
        if fields[0].isdigit():
            atom.element = self.per_table.lookup_element(fields[0])
            atom.charge = int(fields[0])
        else:
            atom.charge = self.per_table.lookup_zcharge(fields[0])
            atom.element = fields[0]
            
        atom.coords = map(float, fields[1:4])
        return Atom
        
    def write_atom(self, atom, outp):
        outp.write(
            self.__class__._atom_line_fmt.format(
                a = atom
            )
        )
        
    def read_coord(self, inp, units = 'ang', header = True):
        if header:
            n_atoms = inp.readline()
            comment = inp.readline()
        
        coord = d.Coordinates(
            units = units
        )

        for line in inp:
            coord.add_atom(
                self.read_atom(line)
            )
            
        return coord
        
    def write_coord(self, mol, outp):
        pass
        
    def read(self, inp):
        pass
        
class TurbomoleIO(object):
    pass

    
class GaussianOutputIO(object):
    pass
    
    
class ADFOutputIO(object):
    pass

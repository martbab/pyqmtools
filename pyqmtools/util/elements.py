from collections import OrderedDict

_ELEMS= (
    'H',
    'He',
    'Li',
    'Be',
    'B',
    'C',
    'N',
    'O',
    'F',
    'Ne',
    'Na',
    'Mg',
    'Al',
    'Si',
    'P',
    'S',
    'Cl',
    'Ar',
    'K',
    'Ca',
    'Sc',
    'Ti',
    'V',
    'Cr',
    'Mn',
    'Fe',
    'Co',
    'Ni',
    'Cu',
    'Zn',
    'Ga',
    'Ge',
    'As',
    'Se',
    'Br',
    'Kr',
    'Rb',
    'Sr',
    'Y',
    'Zr',
    'Nb',
    'Mo',
    'Tc',
    'Ru',
    'Rh',
    'Pd',
    'Ag',
    'Cd',
    'In',
    'Sn',
    'Sb',
    'Te',
    'I',
    'Xe',
    'Cs',
    'Ba',
    'La',
    'Ce',
    'Pr',
    'Nd',
    'Pm',
    'Sm',
    'Eu',
    'Gd',
    'Tb',
    'Dy',
    'Ho',
    'Er',
    'Tm',
    'Yb',
    'Lu',
    'Hf',
    'Ta',
    'W',
    'Re',
    'Os',
    'Ir',
    'Pt',
    'Au',
    'Hg',
    'Tl',
    'Pb',
    'Bi',
    'Po',
    'At',
    'Rn',
    'Fr',
    'Ra',
    'Ac',
    'Th',
    'Pa',
    'U',
    'Np',
    'Pu',
    'Am',
    'Cm',
    'Bk',
    'Cf',
    'Es',
    'Fm',
    'Md',
    'No',
    'Lr',
    'Rf',
    'Db',
    'Sg',
    'Bh',
    'Hs',
    'Mt',
    'Ds',
    'Rg',
    'Cn'
)

_UNK_ELEM = 'Xx'

class PeriodicTable(object):
    def __init__(self,
    ):
        self._init_dicts()
        
    def _init_dicts(self):
        self.elem_zcharge = {}
        self.zcharge_elem = {}
                
        for (z, el) in enumerate(_ELEMS):
            self.elem_zcharge[el] = z + 1
            self.zcharge_elem[z + 1] = el
    
    def lookup_elem(self, z):
        if isinstance(z, int):
            if z in self.zcharge_elem:
                return self.zcharge_elem[z]
            else:
                return _UNK_ELEM
        else:
            raise ValueError(
                "Nuclear charge must be of type <int>!"
            )
    
    def lookup_zcharge(self, elem):
        if isinstance(elem, str) and len(elem) < 3:
            if elem in self.elem_zcharge:
                return self.elem_zcharge[elem]
            else:
                return 0
        else:
            raise ValueError(
                "Element symbol must be a string of maximum 2 characters!"
            )
            
    def lookup(self, arg):
        if isinstance(arg, int):
            return self.lookup_elem(arg)
        elif isinstance(arg, str):
            return self.lookup_zcharge(arg)
        else:
            raise ValueError(
                "Argument must be of type <int> or <str>!"
            )

    
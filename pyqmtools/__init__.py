__version__ = '0.1'

from os.path import splitext
import re
from .geom.datastruct import Atom, Coordinates
from .geom.io import GeomReadException, XYZIO, TurbomoleIO, GaussianOutputIO, ADFOutputIO
from .nmr.datastruct import SigmaTensor, SigmaReference, TensorList, TensorList, TensorStats
from .nmr.parsers import NMRTensorReadError, NMRFinishReadException, GaussianOutputParser, ADFOutputParser
from .util.elements import PeriodicTable
from .util.units import *

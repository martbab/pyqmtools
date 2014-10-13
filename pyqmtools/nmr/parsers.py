from datastruct import *
import re

class NMRTensorReadError(Exception):
    pass

class NMRFinishReadException(Exception):
    pass
    
class GaussianOutputParser(object):
    '''
    class for parsing Gaussian logfile and loading data to TensorList
    '''
    _section_begin = "SCF GIAO Magnetic shielding tensor (ppm):"
    _section_end = "End of Minotr Frequency-dependent properties file"
    _tensor_begin = "Isotropic ="
    _eigenvalues  = "Eigenvalues:"
    _eigenvectors = "Eigenvectors:"

    def __init__(
        self,
        filename = '',
        max_index = 0,
        **kwargs
    ):
        self.filename = filename
        self.shielding_type = 'total'
        self.max_index = max_index
        
    def _g0x_float(self, v):
        if v == "*" * len(v):
            return 0.0
        else:
            return float(v)
    
    def _process_block(
        self,
        block,
        check_index = False
    ):
        #print block
        first_line = block[0].split()
        result = SigmaTensor()
        result.index = int(first_line[0])
        result.element = first_line[1]
        
        if check_index and (result.index > self.max_index):
            raise NMRFinishReadException

        if self.__class__._eigenvalues in block[4][:15]:
            eigvals = (

            )
            result.eigenvalues = array(
                [
                    self._g0x_float(block[4][15:26]), 
                    self._g0x_float(block[4][26:37]), 
                    self._g0x_float(block[4][37:48])
                ]
            )
        else:
            raise NMRTensorReadError(
                "Invalid data format"
            )

        if (len(block) > 5):    
            if self.__class__._eigenvectors in block[5][:17]:
                for i in xrange(6,9):
                    result.eigenvectors = array(
                        [
                            map(
                                self._g0x_float,
                                [
                                    block[i][9:20],
                                    block[i][20:31],
                                    block[i][31:42]
                                ]
                            ) for i in xrange(6,9)
                        ]
                    )

            else:
                raise NMRTensorReadError(
                    "Invalid data format"
                )

        return result
    
    def read(
        self
    ):
        section_found = False

        result = TensorList()
        result.shielding_type = self.shielding_type
        result.filename = self.filename
        result.file_type = "Gaussian 0X output"
        tensor_block = None
        check_index = False
        if self.max_index > 0:
            check_index = True
            

        with open(self.filename, 'r') as f:
            try:
                for line in f:
                    if self.__class__._section_begin in line:
                        section_found = True
                        continue

                    if section_found:
                        if self.__class__._section_end in line:
                            result.append(
                                self._process_block(
                                    tensor_block
                                )
                            )
                            break

                        elif self.__class__._tensor_begin in line: 
                            if tensor_block is None:
                                tensor_block = []
                            else:
                                result.append(
                                    self._process_block(
                                        tensor_block,
                                        check_index = check_index
                                    )
                                )
                                tensor_block = []

                        tensor_block.append(
                            line
                        )
            except NMRFinishReadException:
                return result
            except (TypeError, ValueError, IndexError):
                raise NMRTensorReadError("Failed to read Gaussian NMR tensor")

        return result


class ADFOutputParser(object):
    _file_type = "ADF NMR output"
    _nucleus_blk_begin = "****  N U C L E U S : "
    _nucleus_blk_end = "*" * 79
    _atom_numbering = {
        'input' : "Atom input number in the ADF calculation:",
        'internal' : "Internal NMR numbering of atoms:"
    }
    _shielding_blk_delim = "=== SCALED:"
    _shielding_types = {
        'total' : "TOTAL",
        'spin-orbit' : "SPIN-ORBIT",
        'diamagnetic' : "DIAMAGNETIC",
        'paramagnetic' : "PARAMAGNETIC"
    }
    _job_type_blk_begin = '(INPUT FILE)'
    _job_type_blk_end = 'end'
    _outp_type_token = 'out'
    _iso_total_shielding = 'total'
    _principal_axis_rep = "PRINCIPAL AXIS REPRESENTATION"
    _principal_components = "==== Principal components:"
    _pas = "==== Principal Axis System:"
    _principal_end = "-" * 35
    _atom_number_regexp = re.compile(
        r'\s*(?P<elem>[A-Za-z]{1,2})\((?P<index>\d+)\)\s*'
    )
    _nmr_end = 'N M R   E X I T'

    def __init__(
        self,
        filename,
        shielding_type = 'total',
        atom_numbering = 'input',
    ):
        self.filename = filename

        self.shielding_type = shielding_type
        self.atom_numbering = atom_numbering
        self.output_types = {
            'iso' : self._parse_iso_block,
            'tens' : self._parse_tens_block
        }
        self.outp_type = 'iso'

    def set_shielding_type(self, arg):
        if arg in self.__class__._shielding_types:
            self.__shielding_type = arg
        else:
            raise ValueError(
                "Unrecognized shielding type \"%s\"" % arg
            )

    def get_shielding_type(self):
        return self.__shielding_type

    def process_block(
        self,
        block,
    ):
        result = SigmaTensor()

        (result.index, result.element) = ADFOutputParser.parse_element_index(
            block[0][1]
        )
        if self.max_index > 0 and result.index == self.max_index:
            return None
        result.shielding_type = self.shielding_type

        result.eigenvalues = array(
            map(
                float,
                block[1][0:4]
            )
        )
            
            
        principal_axes = array(
            [
                map(float, i)
                for i in block[3:7]
            ]
        )

        result.eigenvectors = principal_axes.T

        return result

    def parse_element_index(
        self,
        arg
    ):
        match = ADFOutputParser._atom_number_regexp.match(arg)

        if match:
            return (
                int(match.group('index')), 
                match.group('elem')
            )
        else:
            raise NMRTensorReadError(
                "Invalid atom numbering format in ADF output"
            )

    def read(
        self,
    ):

        result = TensorList(
            filename = self.filename,
            file_type = ADFOutputParser._file_type,
        )      

        with open(self.filename, 'r') as inp_file:

            line = '\n'

            while line != '':
                line = inp_file.readline()

                if self.__class__._nmr_end in line:
                    break

                if ADFOutputParser._job_type_blk_begin in line:
                    self._check_outp_type(inp_file)

                if ADFOutputParser._nucleus_blk_begin in line:
                    result.append(
                        self._process_block(
                            inp_file, 
                        )
                    )
            
        result.shielding_type = self.shielding_type
        return result
     
    def _check_outp_type(
        self,
        f
    ):
        line = ""

        fields = []

        while line != ADFOutputParser._job_type_blk_end:
            line = f.readline().strip().lower()

            fields = line.split()

            if fields[0] == ADFOutputParser._outp_type_token:
                if fields[1] in self.output_types:
                    self.outp_type = fields[1]



    def _process_block(
        self,
        f,
    ):
        line = ''
        block = []
        shield_blk_found = False

        result = SigmaTensor()
        shielding_type_line = ' '.join(
            [
                self.__class__._shielding_blk_delim,
                self.__class__._shielding_types[self.shielding_type]
            ]
        )

        while not ADFOutputParser._nucleus_blk_end in line:
            line = f.readline().strip()

            if ADFOutputParser._atom_numbering[self.atom_numbering] in line:
                (result.index, result.element) = \
                    self.parse_element_index(
                        line.split(':')[1]
                    )
                continue

            if shielding_type_line in line:
                shield_blk_found = True
                continue

            if shield_blk_found:
                if ADFOutputParser._shielding_blk_delim in line:
                    break
                elif len(line) != 0:
                    block.append(
                        line
                    )
                else:
                    continue

        #print block
        self.output_types[self.outp_type](
            block,
            result
        )
        return result
            

    def _parse_iso_block(
        self,
        block,
        tensor
    ):
            
        for l in block:
            if self.__class__._iso_total_shielding in l:
                tensor.eigenvalues = array(
                    [float(l.split()[-1])] * 3
                )

    def _parse_tens_block(
        self,
        block,
        tensor
    ):
        try:
            pcomp_idx = block.index(
                self.__class__._principal_components
            ) + 1
                
            tensor.eigenvalues = array(
                map(float, block[pcomp_idx].split())
            )

            
            pas_idx = block.index(
                self.__class__._pas
            ) + 1

            tensor.eigenvectors = array(
                [
                    map(float, i.split()) for i in block[pas_idx:pas_idx + 3]
                ]
            ).T

        except (ValueError, IndexError):
            raise NMRTensorReadError(
                "Failed to read ADF tensor"
            )

    shielding_type = property(get_shielding_type, set_shielding_type)

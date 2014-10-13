from numpy import array, sum, dot, zeros, std, mean, sqrt
from numpy.linalg import norm, eig
from collections import MutableSequence

class SigmaTensor(object):
    '''
    TODO:
    '''

    _index_format = "%6d"
    _element_format = "%3s"
    _eigenvalue_format = "%10.3f"
    _eigenvector_format = "%10.6f"

    def __init__(self, 
        element = "", 
        index = 0, 
        shielding_type = 'total',
    ):

        self.shielding_type = shielding_type
        self.cartesian_rep = zeros((3,3))
        self.symm_tensor = zeros((3,3))

        self.eigenvalues = zeros(3)
        self.eigenvectors = zeros((3,3))

        self.element = element
        self.index = index 

    def __eq__(self, other):
        return self.index == other.index and self.element == other.element

    def __ne__(self, other):
        return self.index != other.index or self.element != other.element

    def __lt__(self, other):
        return self.index < other.index

    def __le__(self, other):
        return self.index <= other.index

    def __gt__(self, other):
        return self.index > other.index

    def __ge__(self, other):
        return self.index >= other.index

    def _calc_symmetric_tensor(self):
        symm_elem = 0.0
        for i in xrange(3):
            for j in xrange(3):
                if i == j:
                    self.symm_tensor[i][j] == self.cartesian_rep 
                else:
                    symm_element = (self.cartesian_rep[i][j] \
                        + self.cartesian_rep[j][i]) * 0.5

                    self.symm_tensor[i][j] = symm_element
                    self.symm_tensor[j][i] = symm_element

    def _calc_eigvals_eigvecs(self):
        (vals, vecs) = eig(self.symm_tensor)
        self.eigenvalues = vals
        # transpose eigenvectors to obtain them in row form
        # eig returns them in column matrix
        self.eigenvectors = vecs.T
        


    def __str__(self):
        '''
        prints out a short summary about the tensor parameters as string
        '''
        message = "Index: " + repr(self.index) + "\n"
        message += "Element symbol: "
        message += self.element
        message += "\n"
        message += "Type of shielding tensor: " + self.shielding_type
        message += "\n"
        message += "Eigenvalues: " + repr(self.eigenvalues) + "\n"
        message += "Isotropic value: " + repr(self.eigenvalues.mean()) + "\n"
        message += "Eigenvectors: " + repr(self.eigenvectors) + "\n"
        return message

    def get_element(self):
        return self.__element

    def set_element(self, arg):
        if (type(arg) is str) and (len(arg) < 3):
            self.__element = arg.capitalize()
        else:
            raise ValueError(
                "Invalid argument for element symbol %s" % repr(arg)
            )

    def write_to_file(self, 
        f, 
        level = 1
    ):
        '''write tensor data to file with various levels of verbosity:
'level = 1' prints out atom index, element symbol and isotropic, value of 
shielding, 'level = 2' acts same as 'level 1', but additionaly prints out the 
eigenvalues of shielding tensor, and, finally, 'level = 3' additionaly prints 
out the principal axis system'''

        if level < 1 or level > 3:
            raise TypeError("Invalid verbosity level, should be 1 to 3")

        # information common to all levels
        output_format = ' '.join(
            [
                self.__class__._index_format,
                self.__class__._element_format,
                self.__class__._eigenvalue_format
            ]
        )
        f.write(
            output_format % (self.index, self.element, self.eigenvalues.mean())
        )

        if level == 1:
            f.write('\n')

        if level > 1:
            outp_format = ' '.join(
                [self.__class__._eigenvalue_format] * 3
            )
            f.write(
                outp_format % tuple(self.eigenvalues)
            )
            f.write('\n')

        if level == 3:
            f.write(
                "# Principal axis system (eigenvectors in rows):\n"
            )
            outp_format = self.__class__._eigenvector_format * 3

            for v in self.eigenvectors:
                f.write("    ")
                f.write(
                    outp_format % tuple(v)
                )
                f.write('\n')

            f.write('#' * 78 + '\n')

    def get_sigma_iso(self):
        return self.eigenvalues.mean()
    
    
    element = property(get_element, set_element)

class SigmaReference():
    '''
    class used for representation of secondary references used for 
    conversion of shielding to chemical shift
    '''

    def __init__(self):
        self.refs = {}

    def __str__(self):
        msg = '\n'.join(
            [
                self.__class__.__name__,
                'Secondary chemical shift references:',
            ]
        )

        for r in sorted(self.refs):
            msg = '\n'.join(
                [
                    msg,
                    'Nucleus: %s sigma: %10.3f delta: %10.3f' % \
                        (r, self.refs[r][0], self.refs[r][1])
                ]
            )

        return msg

    def insert_reference(
        self,
        element,
        sigma,
        delta
    ):
        '''
        insert new reference as a pair of (shielding, shift) values
        '''
        self.refs[element] = (sigma, delta)

    def transform_tensor(self, tensor):
        if tensor.element in self.refs:
            tensor.eigenvalues = \
                self.refs[tensor.element][0] - tensor.eigenvalues + \
                self.refs[tensor.element][1]

    def transform_tensor_list(self, tens_list):
        if isinstance(tens_list, TensorList):
            if tens_list.referenced == False:
                for d in tens_list:
                    self.transform_tensor(d)
                tens_list.referenced = True
        else:
            raise TypeError(
                "Argument must be of type \'TensorList\'"
            )
   
    def read_from_file(self, f):
        '''read secondary reference values from external file in simple form:
        '<element symbol> <calc'd shielding constant> <exp. chemical shift>'
        Lines beginning with '#' are treated as comment and ignored
        '''

        for line in f:
            if line[0] == '#':
                continue

            fields = line.split()
            try:
                self.insert_reference(
                    fields[0],
                    float(fields[1]),
                    float(fields[2])
                )
            except (IndexError, ValueError):
                raise ValueError("Invalid record: \"%s\"" % line)


class TensorList(MutableSequence):
    '''
    TODO:
    '''

    def __init__(self, 
        filename = "", 
        file_type = "", 
        shielding_type = 'total',
    ):
        self.data = list()
        # the name of the file the data are from
        self.filename = filename
        # type of the file the data are from
        self.file_type = file_type
        # shielding type for future use
        self.shielding_type = shielding_type
        # flag whether data were already referenced
        self.referenced = False

    def check_type(self, v):
        if not isinstance(v, SigmaTensor):
            raise TypeError("unsupported type %s" % type(v))

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        self.check_type(value)
        self.data[index] = value

    def __delitem__(self, index):
        del self.data[index]

    def __len__(self):
        return len(self.data)

    def get_shielding_type(self):
        return self.__shielding_type

    def set_shielding_type(self, arg):
        self.__shielding_type = arg

        if len(self.data) != 0:
            for d in self.data:
                d.shielding_type = arg

    def insert(self, index, value):
        self.check_type(value)
        self.data[index] = value

    def append(self, value):
        self.check_type(value)
        self.data.append(value)

    def sort(self):
        self.data.sort()

    def write_to_file(self,
        f,
        level = 1,
        reference = None
    ):
        if level < 0 or level > 3:
            raise TypeError("Invalid verbosity level, should be 1 to 3")

        file_info ='''# NMR shielding tensors from file \"%s\"
# source file type: %s
''' % (self.filename, self.file_type)       

        if reference is not None and isinstance(reference, SigmaReference):
            file_info = ''.join(
                [
                    file_info,
                    '# The shielding values of the following nuclei were \n',
                    '# converted to chemical shifts using \n'
                    '# these secondary references:\n',
                ]
            )

            for r in reference.refs:
                file_info = ''.join(
                    [
                        file_info,
                        "# nucleus: %2s sigma: %10.3f delta: %10.3f\n" %\
                            (r, reference.refs[r][0], reference.refs[r][1])
                    ]
                )

            reference.transform_tensor_list(self)

        f.write(file_info)

        list_info = "# list of %s shielding tensor parameters\n" % \
            self.shielding_type

        header = "# %8s %10s" % ("atom", "isotropic")

        if level > 1:
            header = ''.join(
                [
                    header,
                    "%10s %10s %10s" % ('value 11', 'value 22', 'value 33')
                ]
            )

        if level == 3:
            list_info = '\n'.join(
                [
                    list_info,
                    '# principal axes are also printed',
                    '\n'
                ]
            )

        header = ''.join(
            [
                header,
                '\n'
            ]
        )

        f.write(
            list_info
        )

        f.write(
            header
        )

        for d in self.data:
            d.write_to_file(
                f,
                level = level
            )

    shielding_type = property(get_shielding_type, set_shielding_type)

    
class TensorStats(object):
    """
    class holding multiple instances of SigmaTensor for single atom.
    useful for calculating statistics from multiple calculations
    """
    _stat_header = "#%9s %6s %8s %8s %18s %20s \n" %(
        "Atom",
        "count",
        "mean",
        "sigma",
        "std. err. mean",
        "95% conf. int. (+/-)"
    )
    _stat_line_fmt = "%6d %2s %6d %8.3f %8.3f %18.3f %20.3f\n"
    def __init__(self,
        filenames = [],
        file_type = '',
        shield_type = 'total'
    ):
        self.filenames = filenames
        self.file_type = file_type
        self.shield_type = shield_type
        self.data = {}
        self.stats = {}
        
    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        self.check_type(value)
        self.data[index] = value

    def __delitem__(self, index):
        del self.data[index]

    def __len__(self):
        return len(self.data)
    
    def check_type(self, v):
        if not isinstance(v, SigmaTensor):
            raise TypeError("unsupported type %s" % type(v))
            
    def add_tensor(self, tensor):
        self.check_type(tensor)
        if tensor.index in self.data:
            if tensor.element == self.data[tensor.index][-1].element:
                self.data[tensor.index].append(
                    tensor
                )
            else:
                raise ValueError(
                    "Element mismatch between tensors"
                )
        else:
            self.data[tensor.index] = [tensor]
    
    def add_tensors_from_list(self, tens_list):
        if isinstance(tens_list, TensorList):
            for t in tens_list:
                self.add_tensor(t)
    
    def write_tensors(self, index, outp_file, verb_level = 1):
        outp_file.write(
            "# Array NMR Shielding/shift tensors for atom No. %d\n" % index
        )
        outp_file.write(
            "# %6s %10s %8s %8s %8s %8s\n" % (
                'Number',
                'atom',
                'd_iso',
                'd_11',
                'd_22',
                'd_33'
            )
        )
        
        for (i, t) in enumerate(self.data[index]):
            outp_file.write(
                '%6d ' % (i + 1)
            )
            t.write_to_file(outp_file, level = verb_level)
    
    def write_header(self, outp_file):
        outp_file.write(
            self.__class__._stat_header
        )
        
    def write_tensor_stats(self, index, outp_file):
        iso_values = array([t.get_sigma_iso() for t in self.data[index]])
        
        mean = iso_values.mean()
        std_dev = iso_values.std(ddof = 1)
        std_err_mean = std_dev / sqrt(iso_values.size)
        confidence_int = std_err_mean * 1.96
        outp_file.write(
            self.__class__._stat_line_fmt % (
                index,
                self.data[index][0].element,
                iso_values.size,
                mean,
                std_dev,
                std_err_mean,
                confidence_int
            )
        )
        
    def write_stats(self, outp_file):
        self.write_header(outp_file)
        
        for k in sorted(self.data):
            self.write_tensor_stats(k, outp_file)
        
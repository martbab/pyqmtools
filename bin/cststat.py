#!/home/martbab/devel/python/virtual/pyqmtools/bin/python

import pyqmtools as qmt
from optparse import OptionParser
import sys
from time import strftime, localtime
import os
from numpy import mean, std

def print_time():
    return strftime("%a %d %b %Y %H:%M:%S", localtime())

def print_info(msg, quiet = False):
    if not quiet:
        message = ''.join(
            [
                '<', 
                'CSTStat',
                ' ',
                print_time(),
                ' INFO',
                '>  ',
                str(msg),
                '\n'
            ]
        )

        sys.stdout.write(
            message
        )

def print_error(msg):
    message = ''.join(
        [
            '<', 
            'CSTStat',
            ' ',
            print_time(),
            ' ERROR',
            '> ',
            str(msg),
            '\n'
        ]
    )
    sys.exit(
        message
    )

def extract_csts(
    inp_filenames,
    series_dir,
    cst_parser,
    stat_filename = 'cststat.txt',
    verbose_level = 1,
    reference = None,
    quiet = False,
    max_index = 0
):
    tens_list = None
    tens_stat = qmt.nmr.datastruct.TensorStats(
        filenames = inp_filenames,
    )

    for fn in inp_filenames:
        try:
            print_info(
                "Attempting to read file \"%s\"..." % fn,
                quiet
            )
            cst_parser.filename = fn
            cst_parser.max_index = max_index
            tens_list = cst_parser.read()
            if not quiet:
                print_info(
                    "Processed %d entries." % len(tens_list),
                    quiet
                )
            tens_list.sort()

            if reference is not None: 
                print_info(
                    "Referencing entries...",
                    quiet
                )
                reference.transform_tensor_list(tens_list)

                print_info(
                    "Done.",
                    quiet
                )

            tens_stat.add_tensors_from_list(
                tens_list
            )

        except qmt.nmr.parsers.NMRTensorReadError, e:
            print_error(e)

        except IOError:
            print_error(
                "I/O error during writing."
            )
    if not quiet:
        print_info(
            "Calculating statistics and writing entries to file \"%s\"..." % \
                stat_filename,
             quiet
        )
    with open(stat_filename, 'w') as f:
        tens_stat.write_stats(
            f
        )
    print_info(
        "Success.",
        quiet
    )
    if series_dir != '':
        print_info(
            "Writing individual samples to files in directory \"%s\"" % \
                series_dir,
             quiet
        )
        fname = ""
        if not os.path.isdir(series_dir):
            os.mkdir(series_dir)
            
        for k in tens_stat.data:

            fname = "%06d%s_samples.txt" % (
                k, 
                tens_stat.data[k][0].element
            )

            full_path_name = os.path.join(
                series_dir,
                fname
            )
            with open(full_path_name, 'w') as f:
                tens_stat.write_tensors(
                    k,
                    f
                )
        print_info(
            "Success.",
            quiet
        )
                
        

def extract_reference(
    reference_filename,
    quiet = False
):
        print_info(
            "Reading file \"%s\" containing secondary references..." %\
                reference_filename,
            quiet
        )
        reference = qmt.nmr.datastruct.SigmaReference()

        with open(reference_filename, 'r') as ref_file:
            try:
                reference.read_from_file(ref_file)
            except ValueError:
                print_error("Invalid format of reference file. Aborting...")

        print_info(
            "Found these reference values.",
            quiet
        )
        print_info(reference)

        print_info(
            "Other nuclei will NOT be referenced.",
            quiet
        )
        return reference


def main():
    desc='''A script to process a series of Gaussian or ADF output files and 
print a number statistical descriptors (sample mean, sample standard deviation,
 standard error of the mean) for each nucleus.'''

    usage = '''Usage: %prog [options] input_files'''

    file_types_parsers = {
        'gaussian' : qmt.nmr.parsers.GaussianOutputParser,
        'adf' : qmt.nmr.parsers.ADFOutputParser
    }

    shielding_types = qmt.nmr.parsers.ADFOutputParser._shielding_types
    numbering_types = qmt.nmr.parsers.ADFOutputParser._atom_numbering

    opt_parser = OptionParser(
        usage = usage,
        description = desc
    )

    opt_parser.add_option(
        '-o',
        '--output-filename',
        dest = 'outp_file',
        help = '''Name of the main output file containing statistics. 
Defaults to \'cststat.txt\'.''',
        default = 'cststat.txt',
        metavar = 'FILENAME',
    )
    opt_parser.add_option(
        '-d',
        '--series-directory',
        dest = 'series_dir',
        help = '''Optional directory which can store values of all samples for 
each nucleus as a separate file. Can potentially result in a large number of 
files.''',
        default = '',
        metavar = 'SERIES_DIR'
    )

    opt_parser.add_option(
        '-r',
        '--reference-file',
        dest = 'ref_filename',
        help = qmt.nmr.datastruct.SigmaReference.read_from_file.__doc__,
        default = None,
        metavar = 'FILENAME'
    )

    opt_parser.add_option(
        '-t',
        '--file-type',
        dest = 'file_type',
        type = 'choice',
        choices = file_types_parsers.keys(),
        help = '''type of output file to be processed, currently implemented are 
ADF and Gaussian. Default is ADF.''',
        default = 'adf'
    )

    opt_parser.add_option(
        '-s',
        '--shielding-type',
        dest = 'shield_type',
        type = 'choice',
        choices = shielding_types.keys(),
        help = '''In case of ADF NMR output, sets whether total, diamagnetic, 
paramagnetic or spin-orbit shielding tensor is printed out. 
Defaults to total shielding tensor''',
        default = 'total'
    )
    opt_parser.add_option(
        '-m',
        '--max-index',
        dest = 'max_index',
        type = 'int',
        help = '''Read only atoms up to certain atomic index. Useful if the 
system contains molecules which are of no interest (like solvent for example).
''',
        default = 0,
        metavar = 'INDEX'
    )
#    opt_parser.add_option(
#        '-l',
#        '--verbosity-level',
#        dest = 'verbosity_level',
#        type = 'choice',
#        choices = ('1', '2', '3'),
#        help = qmtools.nmr.datastruct.SigmaTensor.write_to_file.__doc__,
#        default = '1'
#    )

    opt_parser.add_option(
        '-q',
        '--quiet',
        dest = 'quiet',
        action = 'store_true',
        help = "Suppress the amount of output from program.",
        default = False
    )

    opt_parser.add_option(
        '-n',
        '--numbering-type',
        dest = 'numbering_type',
        type = 'choice',
        choices = numbering_types.keys(),
        help = '''In case of ADF NMR output 
(yeah, ADF output files are major pain in the ass) specify whether use atom 
numbering specified by user with input structure or the numbering assigned 
by ADF during fragment generation.
Default is to use input numbering''',
        default = 'input'
    )

    (options, args) = opt_parser.parse_args()


    if len(args) == 0:
        opt_parser.error('No input file name specified!')
    
    print_info(
        'CSTStat v 0.1:'
    )
    print_info(
        'A script to calculate statistics from Gaussian or ADF NMR outputs.' 
    )

    print_info("invoked with following command-line options:")

    for o in options.__dict__.items():
        print_info(
            "%s" % str(o)
        )
        
    reference = None
    if options.ref_filename is not None:
        reference = extract_reference(options.ref_filename, options.quiet)

    cst_parser = file_types_parsers[options.file_type](
        '',
        shielding_type = options.shield_type,
        atom_numbering = options.numbering_type,
    )

    extract_csts(
        args,
        options.series_dir,
        cst_parser,
        stat_filename = options.outp_file,
        reference = reference,
        quiet = options.quiet,
        max_index = options.max_index
    )

    print_info("Finished.")

if __name__ == '__main__':
    main()

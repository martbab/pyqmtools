#!/home/martbab/devel/python/virtual/pyqmtools/bin/python

import pyqmtools as qmt
from optparse import OptionParser
import sys
from time import strftime, localtime
import os

def print_time():
    return strftime("%a %d %b %Y %H:%M:%S", localtime())

def print_info(msg):
    message = ''.join(
        [
            '<', 
            'CSTExtract',
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
            'CSTExtract',
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


def main():
    desc='''
A script to extract shielding tensor information from output of Gaussian 
and ADF QM calculations.'''

    usage = '''Usage: %prog [options] input_file [output_file] 
Writes data to standard output when output file is not specified'''

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
        '-l',
        '--verbosity-level',
        dest = 'verbosity_level',
        type = 'choice',
        choices = ('1', '2', '3'),
        help = qmt.nmr.datastruct.SigmaTensor.write_to_file.__doc__,
        default = '1'
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

    outp_file = None
    tens_list = None

    if len(args) == 0:
        opt_parser.error('No input file name specified!')
    
    print_info(
        'CSTExtract v 0.1:'
    )
    print_info(
        'A script to process results of Gaussian or ADF NMR calculations' 
    )

    print_info("invoked with following command-line options:")

    for o in options.__dict__.items():
        print_info(
            "%s" % str(o)
        )
        
    print_info(
        "Attempting to read file \"%s\"..." % args[0]
    )
    try:
        p = file_types_parsers[options.file_type](
            args[0],
            shielding_type = options.shield_type,
            atom_numbering = options.numbering_type
        )
        tens_list = p.read()
        tens_list.sort()

    except qmt.nmr.parsers.NMRTensorReadError, e:
        print_error(e)
        

    print_info(
        "Processed %d entries." % len(tens_list)
    )
    reference = None
    if options.ref_filename is not None:
        print_info(
            "Reading file \"%s\" containing secondary references..." %\
            options.ref_filename
        )
        reference = qmt.nmr.datastruct.SigmaReference()

        with open(options.ref_filename, 'r') as ref_file:
            try:
                reference.read_from_file(ref_file)
            except ValueError:
                print_error("Invalid format of reference file. Aborting...")

        print_info(
            "Found these reference values."
        )
        print_info(reference)

        print_info(
            "Other nuclei will NOT be referenced."
        )

        print_info(
            "Referencing entries..."
        )
        reference.transform_tensor_list(tens_list)

        print_info(
            "Done."
        )

    try:
        if len(args) >= 2:
            outp_file = open(args[1], 'w')
            print_info(
                "Writing entries to file \"%s\"..." % args[1]
            )
        else:
            outp_file = sys.stdout
            print_info(
                "Writing entries to standard output..."
            )

        tens_list.write_to_file(
            outp_file,
            level = int(options.verbosity_level),
            reference = reference
        )


        print_info("Successfully written data to file/std. output.")
    except IOError:
        print_error(
            "I/O error during writing."
        )
    finally:
        print_info("Finished.")
        outp_file.close()


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import pyqmtools
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

def extract_csts(
    inp_filenames,
    outp_dir,
    cst_parser,
    suffix = '_cstext.txt',
    verbose_level = 1,
    reference = None
):
    tens_list = None
    if outp_dir != '':
        if not os.path.isdir(outp_dir):
            os.mkdir(outp_dir)

    for fn in inp_filenames:
        try:
            print_info(
                "Attempting to read file \"%s\"..." % fn
            )
            cst_parser.filename = fn
            tens_list = cst_parser.read()
            print_info(
                "Processed %d entries." % len(tens_list)
            )
            tens_list.sort()

            if reference is not None: 
                print_info(
                    "Referencing entries..."
                )
                reference.transform_tensor_list(tens_list)

                print_info(
                    "Done."
                )

            prefix = os.path.splitext(fn)[0]
            outp_filename = prefix + suffix

            if outp_dir != '':
                prefix = os.path.splitext(
                    os.path.basename(fn)
                )[0]
                outp_filename = os.path.join(outp_dir, prefix) + suffix

            print_info(
                "Writing entries to file \"%s\"..." % outp_filename
            )

            with open(outp_filename, 'w') as outp_file:
                tens_list.write_to_file(
                    outp_file,
                    level = verbose_level,
                    reference = reference
                )

                print_info("Success.")

        except qmtools.nmr.parsers.NMRTensorReadError, e:
            print_error(e)

        except IOError:
            print_error(
                "I/O error during writing."
            )

def extract_reference(
    reference_filename
):
        print_info(
            "Reading file \"%s\" containing secondary references..." %\
            reference_filename
        )
        reference = qmtools.nmr.datastruct.SigmaReference()

        with open(reference_filename, 'r') as ref_file:
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
        return reference


def main():
    desc='''A script to extract shielding tensor information from output of 
Gaussian and ADF QM calculations. This version is suited for processing of 
multiple files.'''

    usage = '''Usage: %prog [options] input_files'''

    file_types_parsers = {
        'gaussian' : qmtools.nmr.parsers.GaussianOutputParser,
        'adf' : qmtools.nmr.parsers.ADFOutputParser
    }

    shielding_types = qmtools.nmr.parsers.ADFOutputParser._shielding_types
    numbering_types = qmtools.nmr.parsers.ADFOutputParser._atom_numbering

    opt_parser = OptionParser(
        usage = usage,
        description = desc
    )

    opt_parser.add_option(
        '-o',
        '--output-suffix',
        dest = 'suffix',
        help = '''suffix for output filename(s). Defaults to \'_cstext.txt\'.''',
        default = '_cstext.txt',
        metavar = 'SUFFIX',
    )
    opt_parser.add_option(
        '-d',
        '--output-directory',
        dest = 'outp_dir',
        help = 'optional output directory to which to store results. \
Defaults to empty string, which means that the results will be written \
to the same directories as the input files.',
        default = '',
        metavar = 'OUTP_DIR'
    )

    opt_parser.add_option(
        '-r',
        '--reference-file',
        dest = 'ref_filename',
        help = qmtools.nmr.datastruct.SigmaReference.read_from_file.__doc__,
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
        help = qmtools.nmr.datastruct.SigmaTensor.write_to_file.__doc__,
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
        
    reference = None
    if options.ref_filename is not None:
        reference = extract_reference(options.ref_filename)

    cst_parser = file_types_parsers[options.file_type](
        '',
        shielding_type = options.shield_type,
        atom_numbering = options.numbering_type
    )

    extract_csts(
        args,
        options.outp_dir,
        cst_parser,
        suffix = options.suffix,
        verbose_level = int(options.verbosity_level),
        reference = reference
    )

    print_info("Finished.")

if __name__ == '__main__':
    main()

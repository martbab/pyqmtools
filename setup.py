from distutils.core import setup

setup(
    # Application name:
    name="pyqmtools",

    # Version number (initial):
    version="0.1.0",

    # Application author details:
    author="Martin Babinsky",
    author_email="martbab@gmail.com",

    # Packages
    packages=[
        "pyqmtools",
        "pyqmtools.util",
        "pyqmtools.geom",
        "pyqmtools.nmr",
        "bin", 
        "test"
    ],
    # scripts
    scripts = [
        'bin/cstextract.py',
        'bin/cstextract-m.py',
        'bin/cststat.py',
    ],

    # Include additional files into the package
    include_package_data=True,

    # Details
    url="",

    #
    # license="LICENSE.txt",
    description="""A collection of module and scripts for the extraction of 
molecular properties form the results of QM calculations""",

    # long_description=open("README.txt").read(),

    # Dependent packages (distributions)
    install_requires=[
        "scipy",
    ],
)
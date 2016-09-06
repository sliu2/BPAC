__author__ = 'Sheng Liu'

import imp

current_version = imp.load_source('BPACVersion', 'BPAC/_version.py').__version__

try:
    import numpy
except ImportError:
    raise ImportError("BPAC requires numpy to be installed before starting setup (pip install numpy)")
    
try:
    import pyBigWig
except ImportError:
    raise ImportError("BPAC requires pyBigWig to be installed before starting setup (pip install pyBigWig)")

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

setup(
    name='BPAC',
    version=current_version,
    description='Predict TFBS from sequence and accessibility data',
    long_description=open('README.rst',"rt").read(),
    author='Sheng Liu',
    author_email='sliu96@jhmi.edu',
    url='http://bioinfo.wilmer.jhu.edu/BPAC',
    license='GPLv3',
    packages= [
        'BPAC',
    ],

    install_requires=[
        "numpy", 
        "matplotlib", 
        "pysam >= 0.7.5",
        "pyDNase",
        "pybedtools",
        "pyBigWig",
        "pandas",
        "sklearn",
    ],
    package_data = {'BPAC':["data/*"]},
    
    scripts=[
        "BPAC/scripts/generateAttributes",
        "BPAC/scripts/constructModel",
        "BPAC/scripts/predictWithModel",
        "BPAC/scripts/evalPrediction"
        ],
    
)


__author__ = 'Sheng Liu'

import imp

current_version = imp.load_source('pyTFBS_predictVersion', 'pyTFBS_predict/_version.py').__version__

try:
    import numpy
except ImportError:
    raise ImportError("pyTFBS_predict requires numpy to be installed before starting setup (pip install numpy)")
    
try:
    import pyBigWig
except ImportError:
    raise ImportError("pyTFBS_predict requires pyBigWig to be installed before starting setup (pip install pyBigWig)")

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

setup(
    name='pyTFBS_predict',
    version=current_version,
    description='Predict TFBS from sequence and DNase-seq/ATAC-seq data',
    long_description=open('README.rst',"rt").read(),
    author='Sheng Liu',
    author_email='sliu96@jhmi.edu',
    url='http://bioinfo.wilmer.jhu.edu/pyTFBS_predict',
    license='GPLv3',
    packages= [
        'pyTFBS_predict',
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
    package_data = {'pyTFBS_predict':["data/*"]},
    
    scripts=[
        "pyTFBS_predict/scripts/generateAttributes",
        "pyTFBS_predict/scripts/constructModel",
        "pyTFBS_predict/scripts/predictWithModel",
        "pyTFBS_predict/scripts/evalPrediction"
        ],
    
)


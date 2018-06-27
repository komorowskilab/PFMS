#!/usr/bin/env python
# Time-stamp: <2009-10-16 14:46:53 Tao Liu>

"""Description

Setup script for MACS -- Model Based Analysis for ChIP-Seq data

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  beta
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

import os
import sys
from distutils.core import setup, Extension
try:
    import py2exe
except ImportError:
    pass
try:
    import py2app
except ImportError:
    pass

def main():
    if not float(sys.version[:3])>=2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4! python 2.6.2 is recommended!\n")
        sys.exit(1)
    setup(name="MACS",
          version="1.3.7.1",
          description="Model Based Analysis for ChIP-Seq data",
          author='Yong Zhang; Tao (Foo) Liu',
          author_email='zy@jimmy.harvard.edu; taoliu@jimmy.harvard.edu',
          url='http://liulab.dfci.harvard.edu/MACS/',
          package_dir={'MACS' : 'lib'},
          packages=['MACS', 'MACS.IO'],
          scripts=['bin/macs','bin/elandmulti2bed.py','bin/elandresult2bed.py'],
          console=['bin/macs'],
          app    =['bin/macs'],
          classifiers=[
            'Development Status :: 5 - productive',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            ],
          )

if __name__ == '__main__':
    main()

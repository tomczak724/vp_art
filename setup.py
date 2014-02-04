try:
        from setuptools import setup
except ImportError:
        from distutils.core import setup

setup(
        name = 'vp_art',
        version = 'alpha',
        author = 'Adam, Steven',
        classifiers = ['Development Status :: 3 - Alpha',
                    'Environment :: Console',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: '+\
                        'GNU Library or Lesser General Public License (LGPL)',
                    'Natural Language :: English',
                    'Operating System :: POSIX',
                    'Programming Language :: Python',
                    'Topic :: Scientific/Engineering :: Astronomy',
                    'Topic :: Software Development :: Libraries'],
        description = 'A python module for the reduction of virus-p spectra',
        install_requires = ['numpy >= 1.5.1',
                            'scipy >= 0.8.0',
                            'pyfits >= 2.4.0',
                            'matplotlib >= 1.0.1',
                            'astLib >= 0.4.0'],
        packages = ['vp_art']
        )

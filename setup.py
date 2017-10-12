from setuptools import setup
import os
import versioneer

module_dir = os.path.dirname(os.path.abspath(__file__))


def readme():
    with open(os.path.join(module_dir, 'README.rst')) as f:
        return f.read()


setup(
    name='prlworkflows',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=['prlworkflows'],
    description='atomate workflows and scripts for the Phases Research Lab',
    long_description=readme(),
    install_requires=['atomate'],
    extras_require={
        'dev': [
            'sphinx',
            'sphinx_rtd_theme',
            'pytest',
            'twine',
        ]
    },
    author='Brandon Bocklund',
    author_email='brandonbocklund@gmail.com',
    url='https://github.com/phasesresearchlab/prlworkflows',
    license='MIT',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
    ],
)

from setuptools import setup
import os

module_dir = os.path.dirname(os.path.abspath(__file__))


def readme():
    with open(os.path.join(module_dir, 'README.rst')) as f:
        return f.read()


setup(
    name='PRLWorkflows',
    packages=['PRLWorkflows'],
    description='atomate workflows and scripts for the Phases Research Lab',
    long_description=readme(),
    version='0.0.1',
    install_requires=['atomate'],
    author='Brandon Bocklund',
    author_email='brandonbocklund@gmail.com',
    url='https://gitlab.com/bocklund/prl-workflows',
    )

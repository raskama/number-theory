"""
Simple setup script for the project, using setuptools
"""
from setuptools import setup
from setuptools import find_packages

def get_readme():
    """
    Read from README file
    """
    with open('README.md') as readme:
        return readme.read()

setup(
    name='biquad',
    version='1.0',
    description='Python scripts for computing Pythagoras number of biquadratic fields',
    long_description=get_readme(),
    classifiers=[
        'Programming Language :: Python :: 3.8.5',
    ],
    keywords='number_theory biquadratic sum_of_squares',
    include_package_data=True,
    zip_safe=False,
    packages=['biquad','testing'],
)


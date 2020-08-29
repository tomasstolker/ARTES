#!/usr/bin/env python

from setuptools import setup

from pip._internal.network.session import PipSession
from pip._internal.req import parse_requirements

reqs = parse_requirements('requirements.txt', session=PipSession())
reqs = [str(req.requirement) for req in reqs]

setup(
    name='artes',
    version='0.0.1',
    description='Three-dimensional scattering radiative transfer for planetary atmospheres',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    author='Tomas Stolker',
    author_email='tomas.stolker@phys.ethz.ch',
    url='https://github.com/tomasstolker/ARTES',
    project_urls={'Documentation': 'https://species.readthedocs.io'},
    packages=['artes'],
    package_dir={'artes': 'artes'},
    include_package_data=True,
    install_requires=reqs,
    license='MIT',
    zip_safe=False,
    keywords='artes',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)

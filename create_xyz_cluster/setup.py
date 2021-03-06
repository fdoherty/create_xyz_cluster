import io
import os
import re

# from setuptools import find_packages
from setuptools import setup


def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding='utf-8') as fd:
        return re.sub(text_type(r':[a-z]+:`~?(.*?)`'), text_type(r'``\1``'), fd.read())


setup(
    name="create_xyz_cluster",
    version="1.0.0",
    url="https://github.com/fdoherty/create_xyz_cluster",
    license='MIT',

    author="Frank Doherty",
    author_email="fdoherty@umich.edu",

    description="A simple python script to build a xyz coordinate list of atoms from a genetic algorithm output",
    long_description=read("README.rst"),

    packages=['create_xyz_cluster', "tests"],

    package_data={'create_xyz_cluster': ["data/*.xyz"]
                  },

    entry_points={'console_scripts': ['create_xyz_cluster = create_xyz_cluster.create_xyz_cluster:main',
                                      ],
                  },
    package_dir={'create_xyz_cluster': 'create_xyz_cluster'},
#    packages=find_packages(exclude=('tests',)),

    install_requires=['numpy'],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)

from setuptools import setup


import os
setup(
    name='ocelot',
    version='0.01',
    description="",
    author="qai222",
    author_email='qai222@uky.edu',
    url='',
    packages=['ocelot', 'ocelot.schema', 'ocelot.routines', 'ocelot.task'],
    install_requires=['pymatgen', 'numpy', 'shapely', 'scipy', 'rdkit', 'networkx', 'matplotlib', 'numba'],
    keywords='ocelot',
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)

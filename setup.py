try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
setup(
    name='ocelot',
    version='0.2',
    description="",
    author="Ai, Qianxiang",
    author_email='qai222@uky.edu',
    packages=['ocelot', 'ocelot.schema', 'ocelot.routines', 'ocelot.task', 'ocelot.curator'],
    keywords='ocelot',
)

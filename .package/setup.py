from setuptools import setup, find_packages
from os.path import dirname, abspath

where = abspath(__file__)
where = dirname(dirname(where))

setup(
    name="phenocellpy",
    version="0.0.9.alpha",
    packages=find_packages(where=where),
    package_dir={'': '..'},
    install_requires=[
        "numpy >= 1.20",
        "scipy >= 1.6",
    ],
)

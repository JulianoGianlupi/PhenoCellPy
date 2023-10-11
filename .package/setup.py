from setuptools import setup, find_packages
from os.path import dirname, abspath, join

where = abspath(__file__)
where = dirname(dirname(where))

with open(join(where, "README.md"), "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="phenocellpy",
    version="0.0.9a0",
    packages=find_packages(where=where),
    package_dir={'': where},
    long_description=long_description,
    long_description_content_type='text/markdown',
    # install_requires=[
    #     "numpy >= 1.20",
    #     "scipy >= 1.6",
    # ],
)

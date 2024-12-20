from setuptools import setup

INSTALL_REQUIRES = ["MDAnalysis"]

TEST_REQUIRES = [
    # testing and coverage
    "pytest",
    "coverage",
    "pytest-cov",
    # to be able to run `python setup.py checkdocs`
    "collective.checkdocs",
    "pygments",
]


with open("README.md", "r") as f:
    long_description = f.read()

with open("HyresRebuilder/__init__.py", "r") as f:
    init = f.readlines()

for line in init:
    if "__version__" in line:
        __version__ = line.split('"')[-2]

setup(
    name="HyresRebuilder",
    version=__version__,
    author="Shanlong Li",
    author_email="shanlongli@umass.edu",
    description="Rebuild atomistic model from HyRes & iConRNA model.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wayuer19/HyresRebuilder",
    download_url="https://github.com/wayuer19/HyresRebuilder/releases",
    platforms="Tested on Ubuntu 22.04",
    packages=["HyresRebuilder"],
    package_dir={'HyresRebuilder':'HyresRebuilder'},
    package_data={
        "HyresRebuilder":["map/*.pdb"]
    },
    install_requires=INSTALL_REQUIRES,
    extras_require={"test": TEST_REQUIRES + INSTALL_REQUIRES,},
    classifiers=[
        # Trove classifiers
        # (https://pypi.python.org/pypi?%3Aaction=list_classifiers)
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Intended Audience :: Science/Research",
    ],
    entry_points={
        'console_scripts': [
            'RNArebuilder = HyresRebuilder.RNA_Rebuilder:RNArebuild',
            'hyresrebuilder = HyresRebuilder.Rebuilder:rebuild',
        ]
    }
)

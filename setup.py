from setuptools import setup  # , find_packages
from io import open
import versioneer

VERSION = versioneer.get_version()
DESCRIPTION = (
    "ANANSE: Prediction of key transcription factors in cell fate "
    "determination using enhancer networks"
)

with open("README.md", encoding="utf-8") as f:
    long_description = f.read().strip("\n")

setup(
    name="ananse",
    version=VERSION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    description=DESCRIPTION,
    author="Quan Xu, Georgios Georgiou, Siebren Frölich, Maarten van der Sande, Jos Smits, Simon van Heeringen",
    author_email="simon.vanheeringen@gmail.com",
    url="https://github.com/vanheeringen-lab/ananse/",
    download_url=f"https://github.com/vanheeringen-lab/ananse/{VERSION}",
    license="MIT",
    packages=["ananse", "ananse.commands"],  # find_packages(),
    scripts=["scripts/ananse"],
    include_package_data=True,
    zip_safe=False,  # This is necessary, otherwise files won't be installed
    classifiers=[
        "Development Status :: 4 Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "adjusttext",
        "dask",
        "genomepy >=0.14",
        "gimmemotifs >=0.17",
        "loguru",
        "matplotlib >=3.3",
        "networkx",
        "numpy >=1.6",
        "openpyxl",
        "pandas",
        "pybedtools",
        "pydot >=1.4.1",
        "pygraphviz >=1.7",
        "pyranges",
        "tables",  # pytables on conda
        "scikit-learn",
        "scipy >=1.9",
        "seaborn",
        "tqdm",
    ],
)

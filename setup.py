from setuptools import setup, find_packages
from io import open
import versioneer

DESCRIPTION = (
    "ANANSE: Prediction of key transcription factors in cell fate "
    "determination using enhancer networks"
)

with open("README.md", encoding="utf-8") as f:
    long_description = f.read().strip("\n")

setup(
    name="ananse",
    version=versioneer.get_version(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    description=DESCRIPTION,
    author="Quan Xu, Georgios Georgiou, Siebren Frölich, Maarten van der Sande, Jos Smits, Simon van Heeringen",
    author_email="simon.vanheeringen@gmail.com",
    url="https://github.com/vanheeringen-lab/ananse/",
    download_url="https://github.com/vanheeringen-lab/ananse/"
    + versioneer.get_version(),
    license="MIT",
    packages=find_packages(),
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
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=[
        "setuptools >=0.7",
        "adjusttext",
        "dask",
        "gimmemotifs >=0.15.3",
        "loguru",
        "networkx",
        "numpy",
        "openpyxl",
        "pandas <1.3",  # tmp version limit. will be fixed in gimmemotifs >0.16.1
        "scipy",
        "scikit-learn",
        "tables",
        "genomepy >=0.9.3",
        "pyranges",
    ],
)

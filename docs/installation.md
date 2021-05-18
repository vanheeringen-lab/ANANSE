## Installation

ANANSE runs on Linux and Windows 10 using the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10). Mac OSX should work, however, as we don't use it ourselves so unexpected issues might pop up. Let us know, so we can try to fix it.

### The easiest way to install

The recommended way to install ANANSE is by using [conda](https://docs.continuum.io/anaconda). Activate the [bioconda](https://bioconda.github.io/) channel if you haven't used bioconda before.
You only have to do this once.

``` bash
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

You can install ANANSE by creating a specific environment:

``` bash
$ conda config --set use_only_tar_bz2 True
$ conda create -n ananse ananse
```

Don't forget to activate the environment with `conda activate ananse` whenever you want to use ANANSE.

``` bash
# Activate the environment before you use ANANSE
$ conda activate ananse
```

### Alternative installation

You can also install ANANSE with `pip`. 

``` bash
$ pip install git+https://github.com/vanheeringen-lab/ANANSE
``` 

### Developers installation

ANANSE uses setuptools for installations. If you wish to develop in ANANSE, create a conda environment with all dependencies (can be achieved by install ANANSE from conda). Next, clone the repository from [git](https://github.com/vanheeringen-lab/ANANSE). After changing directory to the ANANSE repository, run:

``` bash
$ git checkout develop
$ python setup.py develop
$ python setup.py build
```

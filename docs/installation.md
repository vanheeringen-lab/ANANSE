## Installation

ANANSE runs on Linux and Windows 10 using the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10). 
Mac OSX should work as well, but we don't use it ourselves, so unexpected issues might pop up. 
If you have any issues let us know, so we can try to fix it!

### The easiest way to install

The recommended way to install ANANSE is by using [conda](https://docs.continuum.io/anaconda). 
Activate the [bioconda](https://bioconda.github.io/) channel if you haven't used bioconda before.
You only have to do this once.

``` bash
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

You can install ANANSE by creating a specific environment.
You only have to do this once.

``` bash
$ conda create -n ananse ananse
```

Don't forget to activate the environment whenever you want to use ANANSE.

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

If you wish to try the developmental build, or want to help improve ANANSE, you can!
ANANSE development occurs in the branch of the same name.
You can work and run the development code with the following steps:

``` bash
$ git clone https://github.com/vanheeringen-lab/ANANSE.git
$ cd ANANSE
$ git checkout develop
$ conda env create -n ananse_dev -f requirements.yaml
$ conda activate ananse_dev
$ python setup.py develop
```

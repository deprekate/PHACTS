# PHACTS

A simple program to classify the lifestyle of phages.

Introduction
------------
Phage Classification Tool Set (PHACTS) utilizes a novel similarity algorithm and
a supervised Random Forest classifier to make a prediction whether the lifestyle
of a phage, described by its proteome, is virulent or temperate. The similarity 
algorithm creates a training set from phages with known lifestyles and along with
the lifestyle annotation, trains a Random Forest to classify the lifestyle of a
phage. PHACTS predictions are shown to have a 99% precision rate. 

To install `PHACTS`,
```
pip install phacts
```

or

```sh
 git clone https://github.com/deprekate/PHACTS.git
 pip install PHACTS/ --user
```

This is the new python based recoding of my original perl code, and is in an alpha state for now.
* I have included fasta35 binaries for Linux and OSX, that should install automatically.  For 
Windows, you will need to download a windows fasta35 binary, and make it visible on the PATH

I have tried to replicate exactly the functionality of the original published version.  

PHACTS Example
--------------

Run on included sample data, the output is the predicted lifestyle and the class probabilities:
```sh
$ phacts.py tests/NC_001416.faa 
Class       probability  standard deviation
Temperate   0.626273726  0.101584648
```

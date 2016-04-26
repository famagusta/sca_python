### What is this repository for? ###

* This repository is for holding the Command Line Interface in python for the SCA toolbox (http://systems.swmed.edu/rr_lab/sca.html)
We have recreated the toolbox in python for further use in a Pymol plugin for visualization
* Beta - test version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)


### What is SCA?
* SCA is a computation technique to determine potential allosteric sites in protein families
* based on evolutionary coupling between amino acid sequences in proteins. 

* Note that there is a lower bound on number of sequences used in an MSA determined by binomial
* distribution of frequencies of amino acids and KL entropy function. Around 100 is what the
* notes on SCA calculations suggest. 

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies: Procedure for Linux machines:
    1. SciPy has the following dependencies which need to be installed globally.
       
    2. Requirements.txt contains the list of the libraries/packages required. 

The following commands can be used to install all the dependencies:
        
```
#!
sudo apt-get install gfortran libopenblas-dev liblapack-dev
pip install -r <location of requirements.txt file>

```

* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Contact - robinphilip1989@gmail.com

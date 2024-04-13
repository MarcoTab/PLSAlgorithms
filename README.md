# PLSLib v0.1

**PLSLib** is a collection of Python and R scripts demonstrating the different algorithms detailed in *Partial Least Squares Regression and Related Dimension Reduction Methods* by R. Dennis Cook and Liliana Forzani, available [here](https://lforzani.github.io/PLSR-book/).

## Usage

### Download the Files
To download the files, either download the ZIP file by clicking on the green `< > Code` button: ![image](https://github.com/MarcoTab/PLSLib/assets/64563061/0f69178c-fd4c-4c31-b42d-cd411ed09788)

Or use git to clone the repository to your local: 
- cd into the directory where you want to put the scripts
- run `git clone https://github.com/MarcoTab/PLSLib.git`
- the scripts to reproduce tables and figures are found in the directory `examples`.


### Python
* Python must be version 3.8 or newer.
* Install all required packages using the provided file in the top level directory: `pip3 install -r requirements.txt`.
* Run scripts from the top level directory (i.e. the directory where this README is located). 
Example:

```user@machine: ~/path/to/PLSLib$ python3 example/chapter4/mussels.py```

### R
* R must be version 4.2.2 or newer
* Set your working directory to the top level directory in the repository before running scripts. I.e. `setwd("~/path/to/PLSLib")`.

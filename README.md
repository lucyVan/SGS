# Using HMBFL with SGS：

## HMBFL：

The HMBFL folder stores the original script for the higher-order mutation-based fault localization method and two original benchmark


### Using HMBFL：

#### Program Under Test  
* Put the program under test into 
"./HMBFLtool/test_data/defect_root/source/".

*  Modify the value of "defect_source_path" in 
 "HMBFLtool/start.sh".


#### test cases
* Put the test cases into "./HMBFLtool/test_data/inputs/"

* if the test cases are read by the program under test by file redirection, the value of the "input_type" in "start.sh" need to be changed to "file", not "args"

#### MBFL-guided HOMs generation strategy
* The specific execution flow is viewed in "HMBFLtool/start_MBFL_randomization_higher_order_MBFL.sh"

#### SBFL-guided HOMs generation strategy
* The specific execution flow is viewed in "HMBFLtool/start_SBFL_randomization_higher_order_MBFL.sh"

#### Entire random HOMs generation strategy
* The specific execution flow is viewed in "HMBFLtool/start_entire_randomization_higher_order_MBFL.sh"

The folder of the respective HOMS generation strategies contains the implementation code for the three suspiciousness statistics methods


## mbfl sbfl result_get:
These are some tools for reduction of mutants.


## SGS.py:

The SGS folder stores the script for SGS mutant reduction method.

### Using SGS.py：
After getting the set of higher-order mutants by the code in the HMBFL directory, the table of experimental results can be obtained by configuring the directory where the higher-order mutants are located in the CHMBFL function.

## Remarks：

* This tool passed experiment on the "Linux version 3.10.0-957.el7.x86_64".
* python version : 3.6.8

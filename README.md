# sugarMassesPredict
command line tool to calculate all possible glycan molecules and the *m/z* values of their ions given a set of input parameters

**NOTE: I have now also added an 'R' version - if you use the reticulate package in R, you can source the file (sugarMassesPredict-r.py) and then use the function 'predict_sugars' directly in R to generate output there. **

this tool is being frequently updated :) if you have any questions or issue please feel free to contact me at mbligh@mpi-bremen.de

_word of caution: please note that this tool will predict sugars that are not possible as the nature of sugar chemistry means that it would take a long time to add in all the constraints!_


e.g.
```
library(reticulate)
py_install("pandas", "numpy")
source_python("sugarMassesPredict-r.py")
dp1 = as.integer(1)
dp2 = as.integer(3)
ESI_mode = 'pos'
scan_range1 = as.integer(100)
scan_range2 = as.integer(800)
pent_option = as.integer(1)
modifications = list('sulphate', 'deoxy')
label = "procainamide"
df <- predict_sugars(dp1 = dp1, dp2 = dp2, ESI_mode = ESI_mode, scan_range1 = scan_range1, scan_range2 = scan_range2, pent_option = pent_option, modifications = modifications, label = label) 
```

# dependencies
* pandas
* numpy
* python 3

# input parameters 
## required
* dp (degree of polymerisation) range
* whether pentose monomers should be used in addition to hexose
* modifications - possible options are none, all, or any combination of:
    * sulphate
    * carboxyl
    * phosphate
    * deoxy
    * N-acetyl
    * O-acetyl
    * O-methyl
    * anhydrobridge
    * unsaturated
    * alditol
    * amino
    * dehydrated
* maximum number of modifications per monomer on average
* ionisation mode
* scan range (*m/z*)
## optional
* label - current options are procainamide (added by reductive amination) and benzoic acid (added on free alcohol groups, will calculate glycans with no label to the maximum number of labels possible)
* output file path - defaults to "predicted_sugars.txt"
* options to do with the calculation of the possible number of structural isomers (but this section needs to be fixed)

# output
tab delimited text file with one row per molecule. *m/z* values outside the scan range as shown as "NA", and molecules with no ions with *m/z* values within the scan are not returned. columns are as follows:
* degree of polymerisation (dp)
* name
* monoisotopic mass (Da)
* sum formula
* columns with *m/* values of possible ions given the input parameters.
    * positive mode: 
        * \[M+H\]<sup>+</sup>
        * \[M+Na\]<sup>+</sup>
    * negative mode: 
      * \[M+Cl\]<sup>-</sup>
      * \[M+CHOO\]<sup>-</sup>
      * \[M+2Cl\]<sup>2-</sup>
      * \[M+2CHOO\]<sup>2-</sup>
      * \[M+Cl-H\]<sup>2-</sup>
      * \[M+CHOO-H\]<sup>2-</sup>
      * \[M+CHOO+Cl\]<sup>2-</sup>
      * \[M-*n*H\]<sup>*n*-</sup>, where n is 1 to the maximum number anionic groups that any single molecule in the table has
* if parameters related to isomers were specified in the input, there is an additional column for the number of possible isomers per molecule

# how to run
the help menu accessed with:
```
sugarMassesPredict.py -h
```

returns the following:
```
usage: sugarMassesPredict.py [-h] -dp int int [-p int] -m str [str ...]
                             [-n int] [-ds int] [-ld str [str ...]] [-oh int]
                             [-b int] -i str [str ...] -s int int [-l label]
                             [-o filepath]

Script to predict possible masses of unknown sugars. Written by Margot Bligh.

optional arguments:
  -h, --help            show this help message and exit
  -dp int int, --dp_range int int
                        DP range to predict within: two space separated
                        numbers required (lower first)
  -p int, --pent_option int
                        should pentose monomers be considered as well as
                        hexose: 0 for no {default}, 1 for yes
  -m str [str ...], --modifications str [str ...]
                        space separated list of modifications to consider.
                        note that alditol and unsaturated are max once per
                        saccharide. allowed values: none OR all OR any
                        combination of carboxyl, phosphate, deoxy, nacetyl,
                        omethyl, anhydrobridge, oacetyl, unsaturated, alditol,
                        sulphate
  -n int, --nmod_max int
                        max no. of modifications per monomer on average
                        {default 1}. does not take into account unsaturated or
                        alditol.
  -ds int, --double_sulphate int
                        can monomers be double-sulphated: 0 for no {default},
                        1 for yes. for this you MUST give a value of at least
                        2 to -n/--nmod_max
  -ld str [str ...], --LorD_isomers str [str ...]
                        isomers calculated for L and/or D enantiomers {default
                        D only}. write space separated if both
  -oh int, --OH_stereo int
                        stereochem of OH groups considered when calculating
                        no. of isomers: 0 for no {default}, 1 for yes
  -b int, --bond_stereo int
                        stereochem of glycosidic bonds and reducing end
                        anomeric carbons considered when calculating no. of
                        isomers: 0 for no {default}, 1 for yes
  -i str [str ...], --ESI_mode str [str ...]
                        neg and/or pos mode for ionisation (space separated if
                        both)
  -s int int, --scan_range int int
                        mass spec scan range to predict within: two space
                        separated numbers required (lower first)
  -l label, --label label
                        name a label added to the oligosaccharide. if not
                        labelled do not include. options: procainamide OR
                        benzoic_acid.
  -o filepath, --output filepath
                        filepath to .txt file for output table {default:
                        predicted_sugars.txt}

```





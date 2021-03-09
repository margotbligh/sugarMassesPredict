# sugarMassesPredict
Command line executable tool to calculate all possible glycan molecules and the m/z values of their ions given a set up input parameters.

# dependencies
* pandas
* numpy
* python 3

# input parameters 
## required
* dp (degree of polymerisation) range
* whether pentose monomers should be used in addition to hexose
* modifications - possible options are none, all, or on any combination of:
  * sulphate
  * carboxyl
  * phosphate
  * deoxy
  * N-acetyl
  * O-acetyl
  * O-methyl
  * anhydrobridge
* maximum number of modifications per monomer on average
* ionisatione mode
* scan range (*m/z*)
## optional
* label - current options are procainamide (added by reductive amination) and benzoic acid (added on free alcohol groups, will calculate glycans with no label to the maximum number of labels possible)
* output file path - defaults to "predicted_sugars.txt"
* options to do with the calculation of the possible number of structural isomers (but this section needs to be fixed)

# how to run
The help menu can be accessed with:
````
sugarMassesPredict.py -h

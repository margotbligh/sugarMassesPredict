#!/usr/bin/env python
# import modules
import pandas as pd
import numpy as np
import itertools
import argparse
import time
import gc
# start timer
start_time = time.time()
# suppress warnings
pd.options.mode.chained_assignment = None  # default='warn'
# DEFINE INPUT ARGUMENTS
# ----------------------
msg = "Script to predict possible masses of unknown sugars. Written by Margot Bligh."

# initialise parser
parser = argparse.ArgumentParser(description=msg)

# add arguments
possible_modifications = ['carboxyl',
                          'phosphate',
                          'deoxy',
                          'nacetyl',
                          'omethyl',
                          'anhydrobridge',
                          'oacetyl',
                          'unsaturated',
                          'alditol',
                          'amino',
                          'dehydrated',
                          'sulphate']
parser.add_argument('-dp',
                    '--dp_range',
                    help='DP range to predict within: two space separated numbers required (lower first)',
                    nargs=2,
                    required=True,
                    type=int,
                    dest='dp_range',
                    metavar='int')

parser.add_argument('-p',
                    '--pent_option',
                    help='should pentose monomers be considered as well as hexose: 0 for no {default}, 1 for yes',
                    nargs=1,
                    type=int,
                    dest='pent_option',
                    metavar='int',
                    default=0,
                    choices=[0, 1])

parser.add_argument('-m',
                    '--modifications',
                    help='space separated list of modifications to consider. note that alditol, dehydrated and unsaturated are max once per saccharide. allowed values: none OR all OR any combination of ' + ', '.join(
                        possible_modifications),
                    nargs='+',
                    dest='modifications',
                    required=True,
                    metavar='str',
                    choices=possible_modifications + ['none', 'all'])

parser.add_argument('-n',
                    '--nmod_max',
                    help='max no. of modifications per monomer on average {default 1}. does not take into account unsaturated, dehydrated or alditol.',
                    nargs=1,
                    type=int,
                    default=1,
                    dest='nmod_max',
                    metavar='int')

parser.add_argument('-ds',
                    '--double_sulphate',
                    help='can monomers be double-sulphated: 0 for no {default}, 1 for yes. for this you MUST give a value of at least 2 to -n/--nmod_max',
                    nargs=1,
                    type=int,
                    default=0,
                    dest='double_sulphate',
                    metavar='int',
                    choices=[0, 1])

parser.add_argument('-ld',
                    '--LorD_isomers',
                    help='isomers calculated for L and/or D enantiomers {default D only}. write space separated if both',
                    nargs='+',
                    metavar='str',
                    dest='LorD_isomers',
                    choices=['L', 'D'],
                    default='D')

parser.add_argument('-oh',
                    '--OH_stereo',
                    help='stereochem of OH groups considered when calculating no. of isomers: 0 for no {default}, 1 for yes',
                    nargs=1,
                    dest='OH_stereo',
                    metavar='int',
                    choices=[0, 1],
                    default=0,
                    type=int)
parser.add_argument('-b',
                    '--bond_stereo',
                    help='stereochem of glycosidic bonds and reducing end anomeric carbons considered when calculating no. of isomers: 0 for no {default}, 1 for yes',
                    nargs=1,
                    dest='bond_stereo',
                    metavar='int',
                    choices=[0, 1],
                    default=0,
                    type=int)
parser.add_argument('-i',
                    '--ESI_mode',
                    help='neg and/or pos mode for ionisation (space separated if both)',
                    nargs='+',
                    required=True,
                    dest='ESI_mode',
                    metavar='str',
                    choices=["neg", "pos"])
parser.add_argument('-s',
                    '--scan_range',
                    help='mass spec scan range to predict within: two space separated numbers required (lower first)',
                    nargs=2,
                    required=True,
                    type=int,
                    dest='scan_range',
                    metavar='int')

parser.add_argument('-l',
                    '--label',
                    help='name a label added to the oligosaccharide. if not labelled do not include. options: procainamide OR benzoic_acid.',
                    metavar='label',
                    dest='label',
                    default='none')

parser.add_argument('-o',
                    '--output',
                    help='filepath to .txt file for output table {default: predicted_sugars.txt}',
                    nargs=1,
                    metavar='filepath',
                    dest='filepath',
                    default='predicted_sugars.txt')

# read arguments from command line
args = parser.parse_args()

dp_range = args.dp_range
pent_option = args.pent_option
if isinstance(pent_option, list):
    pent_option = pent_option[0]
modifications = args.modifications
nmod_max = args.nmod_max
if isinstance(nmod_max, list):
    nmod_max = nmod_max[0]
double_sulphate = args.double_sulphate
if isinstance(double_sulphate, list):
    double_sulphate = double_sulphate[0]
LorD_isomers = args.LorD_isomers
OH_stereo = args.OH_stereo
if isinstance(OH_stereo, list):
    OH_stereo = OH_stereo[0]
bond_stereo = args.bond_stereo
if isinstance(bond_stereo, list):
    bond_stereo = bond_stereo[0]
ESI_mode = args.ESI_mode
scan_range = args.scan_range
label = args.label
outfile = args.filepath
if isinstance(outfile, list):
    outfile = outfile[0]

if "all" in modifications:
    modifications = possible_modifications

if "sulphate" in modifications:
    modifications.append(modifications.pop(modifications.index('sulphate')))

if "alditol" in modifications:
    alditol_option = 'y'
    modifications.remove('alditol')
elif "alditol" not in modifications:
    alditol_option = 'n'

if "unsaturated" in modifications:
    unsaturated_option = 'y'
    modifications.remove('unsaturated')
elif "unsaturated" not in modifications:
    unsaturated_option = 'n'

if "dehydrated" in modifications:
    dehydrated_option = 'y'
    modifications.remove('dehydrated')
elif "dehydrated" not in modifications:
    dehydrated_option = 'n'


# 1: DEFINE MASSES / FORMULAS / ISOMERS / MODIFICATIONS VARIABLES / FUNCTIONS
# ----------------------

print("step #1: defining mass, formula, isomer and modification variables, and functions")
print("-------------------------------------------------------------------------------\n")

# hexose and water masses to build molecule base
hex_mass = 180.06339
water_mass = 18.010565
# mass differences for modifications
pent_mdiff = -30.010566
modifications_mdiff = {
    "sulphate": 79.956817,
    "anhydrobridge": -water_mass,
    "omethyl": 14.01565,
    "carboxyl": 13.979265,
    "nacetyl": 41.026549,
    "oacetyl": 42.010565,
    "phosphate": 79.966333,
    "deoxy": -15.994915,
    "unsaturated": -2.015650,
    "alditol": 2.015650,
    "amino": -0.984016,
    "dehydrated": -water_mass
}

# mass differences for labels
procainamide_mdiff = 219.173546
benzoic_acid_mdiff = 104.026215

# mass differences for ions
ion_mdiff = {
    "H": 1.00782500000003,
    "Na": 22.98977,
    "Cl": 34.968853,
    "CHOO": 44.997655,
    "NH4": 18.034374,
    "K": 38.963708
}

e_mdiff = 0.000548579909

# formulas

formulas = {
    "hex": [6, 12, 0, 6, 0, 0],
    "pent": [5, 10, 0, 5, 0, 0],
    "water": [0, -2, 0, -1, 0, 0],
    "sulphate": [0, 0, 0, 3, 1, 0],
    "anhydrobridge": [0, -2, 0, -1, 0, 0],
    "omethyl": [1, 2, 0, 0, 0, 0],
    "carboxyl": [0, -2, 0, 1, 0, 0],
    "nacetyl": [2, 3, 1, 0, 0, 0],
    "oacetyl": [2, 2, 0, 1, 0, 0],
    "phosphate": [0, 1, 0, 3, 0, 1],
    "deoxy": [0, 0, 0, -1, 0, 0],
    "procainamide": [13, 21, 3, 0, 0, 0],
    "benzoic_acid": [7, 4, 0, 1, 0, 0],
    "unsaturated": [0, -2, 0, 0, 0, 0],
    "alditol": [0, +2, 0, 0, 0, 0],
    "amino": [0, +1, +1, -1, 0, 0],
    "dehydrated": [0, -2, 0, -1, 0, 0]
}
# modification types
modifications_anionic = {"sulphate",
                         "phosphate",
                         "carboxyl"}
modifications_neutral = {"anhydrobridge",
                         "omethyl",
                         "nacetyl",
                         "oacetyl",
                         "deoxy",
                         "unsaturated",
                         "amino",
                         "dehydrated"}

# isomers
isomers_OHdiff = {"anhydrobridge",
                  "omethyl",
                  "nacetyl",
                  "oacetyl",
                  "deoxy"}

# set up general variables and functions
dp_range_list = list(range(dp_range[0], dp_range[1] + 1))


def dpRepeats(dp_range_list):
    repeats_list = []
    for i in dp_range_list:
        repeats_list = repeats_list + list(range(0, i + 1))
    return repeats_list


def bcRepeats(nmax_bc):
    repeats_list = []
    for i in nmax_bc:
        repeats_list = repeats_list + list(range(0, i + 1))
    return repeats_list


def getHexMasses(dp_range_list):
    dp = pd.Series(dp_range_list)
    name = "hex-" + dp.astype(str)
    hex = dp
    mass = dp * hex_mass - (dp - 1) * water_mass
    masses = pd.DataFrame({'dp': dp,
                           'name': name,
                           'hex': hex.astype(int),
                           'mass': mass})
    return masses


def getPentMasses(masses):
    dp = masses.dp.repeat(masses.dp.array + 1).reset_index(drop=True)
    pent = pd.Series(dpRepeats(dp_range_list))
    hex = dp - pent
    name = "hex-" + hex.astype(str) + "-pent-" + pent.astype(str)
    mass = masses.mass.repeat(masses.dp.array + 1).reset_index(drop=True)
    mass = mass + pent * pent_mdiff
    masses = pd.DataFrame({'dp': dp,
                           'name': name,
                           'hex': hex,
                           'pent': pent,
                           'mass': mass})
    return masses


def getModificationNumbers(dp_range_list, m, pent_option, modifications):
    modification_numbers = []
    for i in dp_range_list:
        a = list(range(0, i + 1))
        if pent_option == 1:
            modification_numbers = modification_numbers + \
                                   list(itertools.product(a, repeat=m)) * (i + 1)
        elif pent_option == 0:
            modification_numbers = modification_numbers + \
                                   list(itertools.product(a, repeat=m))
    modification_numbers = pd.DataFrame(modification_numbers)
    modification_numbers.columns = modifications
    return modification_numbers


elapsed_time = time.time() - start_time
print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

# 2: CALCULATE ALL POSSIBLE MASSES
# ----------------------

print("\nstep #2: calculating all possible masses")
print("----------------------------------------\n")

# build hexose molecules
print("--> getting hexose masses")
masses = getHexMasses(dp_range_list)
elapsed_time = time.time() - start_time
print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if "none" in modifications and pent_option == 0 and 'benzoic_acid' in label:
    print("--> adding benzoic acid")
    masses['maxn_benzoic_acid'] = (masses.hex * 5) - (masses.dp * 2 - 2)
    nmax_bc = list(masses.maxn_benzoic_acid)
    masses_array = np.array(masses)
    masses_array = masses_array.repeat(masses.maxn_benzoic_acid.array + 1, axis=0)
    colNames = masses.columns
    masses = pd.DataFrame(masses_array)
    masses.columns = colNames
    n_bc = bcRepeats(nmax_bc)
    masses['benzoic_acid'] = n_bc
    masses.mass = masses.mass + benzoic_acid_mdiff * masses.benzoic_acid
    masses = masses.drop(columns='maxn_benzoic_acid')
    masses.name = masses.name + "-benzoic_acid-" + masses.benzoic_acid.astype(str)
    masses.name = masses.name.str.replace("-benzoic_acid-0", "")
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    del masses_array

# calculate masses for pentose molecules if selected
if pent_option == 1:
    print("--> getting pentose masses")
    masses = getPentMasses(masses)
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if "none" in modifications and pent_option == 1 and 'benzoic_acid' in label:
    print("--> adding benzoic acid")
    masses['maxn_benzoic_acid'] = (masses.hex * 5) + (masses.pent * 4) - (masses.dp * 2 - 2)
    nmax_bc = list(masses.maxn_benzoic_acid)
    masses_array = np.array(masses)
    masses_array = masses_array.repeat(masses.maxn_benzoic_acid.array + 1, axis=0)
    colNames = masses.columns
    masses = pd.DataFrame(masses_array)
    masses.columns = colNames
    n_bc = bcRepeats(nmax_bc)
    masses['benzoic_acid'] = n_bc
    masses.mass = masses.mass + benzoic_acid_mdiff * masses.benzoic_acid
    masses = masses.drop(columns='maxn_benzoic_acid')
    masses.name = masses.name + "-benzoic_acid-" + masses.benzoic_acid.astype(str)
    masses.name = masses.name.str.replace("-benzoic_acid-0", "")
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    del masses_array

# add modifications
if "none" not in modifications and pent_option == 1:
    print("--> adding modifications")
    m = len(modifications)
    dp = masses.dp.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
    hex = masses.hex.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
    pent = masses.pent.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
    modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
    name = "hex-" + hex.astype(str) + "-pent-" + pent.astype(str)
    for i in range(m):
        name = name + "-" + modifications[i] + "-" + modification_numbers[modifications[i]].astype(str)
    name = name.str.replace("-\D+-0", "")
    name = name.str.replace("hex-0-", "")
    mass = masses.mass.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
    for i in range(m):
        mass = mass + modifications_mdiff[modifications[i]] * modification_numbers[modifications[i]]
    masses = pd.DataFrame({'dp': dp,
                           'name': name,
                           'hex': hex,
                           'pent': pent})
    masses = pd.concat([masses, modification_numbers], axis=1)
    masses['mass'] = mass
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    # add benzoic acid
    if "benzoic_acid" in label and "anhydrobridge" not in modifications:
        print("--> adding benzoic acid")
        modification_numbers_sum = modification_numbers.sum(axis=1)
        masses['maxn_benzoic_acid'] = (masses.hex * 5) + (masses.pent * 4) - (
                masses.dp * 2 - 2) - modification_numbers_sum
        nmax_bc = list(masses.maxn_benzoic_acid)
        masses_array = np.array(masses)
        masses_array = masses_array.repeat(masses.maxn_benzoic_acid.array + 1, axis=0)
        colNames = masses.columns
        masses = pd.DataFrame(masses_array)
        masses.columns = colNames
        n_bc = bcRepeats(nmax_bc)
        masses['benzoic_acid'] = n_bc
        masses.mass = masses.mass + benzoic_acid_mdiff * masses.benzoic_acid
        masses = masses.drop(columns='maxn_benzoic_acid')
        masses.name = masses.name + "-benzoic_acid-" + masses.benzoic_acid.astype(str)
        masses.name = masses.name.str.replace("-benzoic_acid-0", "")
        elapsed_time = time.time() - start_time
        print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    if "benzoic_acid" in label and "anhydrobridge" in modifications:
        print("--> adding benzoic acid")
        modification_numbers_sum = modification_numbers.drop(columns="anhydrobridge").sum(
            axis=1) + modification_numbers.anhydrobridge * 2
        nmax_bc = (masses.hex * 5) + (masses.pent * 4) - (masses.dp * 2 - 2) - modification_numbers_sum
        nmax_bc = np.array(nmax_bc)
        nmax_bc[nmax_bc < 0] = 0
        nmax_bc = pd.Series(nmax_bc)
        masses['maxn_benzoic_acid'] = nmax_bc
        nmax_bc = list(masses.maxn_benzoic_acid)
        masses_array = np.array(masses)
        masses_array = masses_array.repeat(masses.maxn_benzoic_acid.array + 1, axis=0)
        colNames = masses.columns
        masses = pd.DataFrame(masses_array)
        masses.columns = colNames
        n_bc = bcRepeats(nmax_bc)
        masses['benzoic_acid'] = n_bc
        masses.mass = masses.mass + benzoic_acid_mdiff * masses.benzoic_acid
        masses = masses.drop(columns='maxn_benzoic_acid')
        masses.name = masses.name + "-benzoic_acid-" + masses.benzoic_acid.astype(str)
        masses.name = masses.name.str.replace("-benzoic_acid-0", "")
        elapsed_time = time.time() - start_time
        print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if "none" not in modifications and pent_option == 0:
    print("--> adding modifications")
    m = len(modifications)
    dp = masses.dp.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
    hex = masses.hex.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
    modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
    name = "hex-" + hex.astype(str)
    for i in range(m):
        name = name + "-" + modifications[i] + "-" + modification_numbers[modifications[i]].astype(str)
    name = name.str.replace("-\D+-0", "")
    name = name.str.replace("hex-0-", "")
    mass = masses.mass.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
    for i in range(m):
        mass = mass + modifications_mdiff[modifications[i]] * modification_numbers[modifications[i]]
    masses = pd.DataFrame({'dp': dp,
                           'name': name,
                           'hex': hex})
    masses = pd.concat([masses, modification_numbers], axis=1)
    masses['mass'] = mass
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    # add benzoic acid
    if "benzoic_acid" in label and "anhydrobridge" not in modifications:
        print("--> adding benzoic acid")
        modification_numbers_sum = modification_numbers.sum(axis=1)
        masses['maxn_benzoic_acid'] = (masses.hex * 5) - (masses.dp * 2 - 2) - modification_numbers_sum
        nmax_bc = list(masses.maxn_benzoic_acid)
        masses_array = np.array(masses)
        masses_array = masses_array.repeat(masses.maxn_benzoic_acid.array + 1, axis=0)
        colNames = masses.columns
        masses = pd.DataFrame(masses_array)
        masses.columns = colNames
        n_bc = bcRepeats(nmax_bc)
        masses['benzoic_acid'] = n_bc
        masses.mass = masses.mass + benzoic_acid_mdiff * masses.benzoic_acid
        masses = masses.drop(columns='maxn_benzoic_acid')
        masses.name = masses.name + "-benzoic_acid-" + masses.benzoic_acid.astype(str)
        masses.name = masses.name.str.replace("-benzoic_acid-0", "")
        elapsed_time = time.time() - start_time
        print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    if "benzoic_acid" in label and "anhydrobridge" in modifications:
        print("--> adding benzoic acid")
        modification_numbers_sum = modification_numbers.drop(columns="anhydrobridge").sum(
            axis=1) + modification_numbers.anhydrobridge * 2
        masses['maxn_benzoic_acid'] = (masses.hex * 5) - (masses.dp * 2 - 2) - modification_numbers_sum
        nmax_bc = list(masses.maxn_benzoic_acid)
        masses_array = np.array(masses)
        masses_array = masses_array.repeat(masses.maxn_benzoic_acid.array + 1, axis=0)
        colNames = masses.columns
        masses = pd.DataFrame(masses_array)
        masses.columns = colNames
        n_bc = bcRepeats(nmax_bc)
        masses['benzoic_acid'] = n_bc
        masses.mass = masses.mass + benzoic_acid_mdiff * masses.benzoic_acid
        masses = masses.drop(columns='maxn_benzoic_acid')
        masses.name = masses.name + "-benzoic_acid-" + masses.benzoic_acid.astype(str)
        masses.name = masses.name.str.replace("-benzoic_acid-0", "")
        elapsed_time = time.time() - start_time
        print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if "sulphate" in modifications and double_sulphate == 1:
    print("--> adding extra sulphate groups")
    masses_s1 = masses.loc[masses['sulphate'] >=1]
    masses_s2 = masses_s1
    masses_s2.sulphate = masses_s1.sulphate + masses_s1.dp
    masses_s2.name = masses_s2.name.str.replace("-sulphate-\d{1,2}", "")
    masses_s2.name = masses_s2.name + '-sulphate-' + masses_s2.sulphate.astype(str)
    masses_s2.mass = masses_s2.mass + modifications_mdiff['sulphate'] * masses_s2.dp
    masses = masses.append(masses_s2).reset_index()
    del masses_s1
    del masses_s2
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if "procainamide" in label:
    print("--> adding procainamide label")
    masses['name'] = masses.name + '-procA'
    masses['mass'] = masses.mass + procainamide_mdiff
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if unsaturated_option == 'y':
    print("--> adding unsaturated sugars")
    masses_a = masses.copy()
    masses_a.name = "unsaturated-" + masses.name
    masses_a['unsaturated'] = 1
    masses['unsaturated'] = 0
    masses_a.mass = masses.mass + modifications_mdiff['unsaturated']
    masses = masses.append(masses_a).reset_index()
    del masses_a
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if alditol_option == 'y':
    print("--> adding alditol sugars")
    masses_a = masses.copy()
    masses_a.name = "alditol-" + masses_a.name
    masses_a['alditol'] = 1
    masses['alditol'] = 0
    masses_a.mass = masses_a.mass + modifications_mdiff['alditol']
    masses = masses.append(masses_a).reset_index(drop = True)
    del masses_a
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

if dehydrated_option == 'y':
    print("--> adding dehydration")
    masses_a = masses.copy()
    masses_a.name = "dehydrated-" + masses_a.name
    masses_a['dehydrated'] = 1
    masses['dehydrated'] = 0
    masses_a.mass = masses_a.mass + modifications_mdiff['dehydrated']
    masses = masses.append(masses_a).reset_index(drop = True)
    del masses_a
    elapsed_time = time.time() - start_time
    print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

gc.collect()

# 3: GET FORMULAS
# ----------------------

print("\nstep #3: building formulas")
print("----------------------------------------\n")

if "none" in modifications and pent_option == 1:
    if "benzoic_acid" in label:
        dp = masses.dp
        hex = masses.hex
        pent = masses.pent
        benzoic_acid = masses.benzoic_acid
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex,
                                         'pent': pent,
                                         'benzoic_acid': benzoic_acid})
    else:
        dp = masses.dp
        hex = masses.hex
        pent = masses.pent
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex,
                                         'pent': pent})

    molecules = list(molecule_numbers.drop('dp', axis=1).columns)
    atom_names = ["C", "H", "N", "O", "S", "P"]
    atom_list = []
    for i in range(len(atom_names)):
        n = np.array([0] * len(masses.index))
        for j in range(len(molecules)):
            form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers[molecules[j]])
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
        if "procainamide" in label:
            p = np.array([formulas['procainamide'][i]] * len(masses.index))
            n = n + p
        atom_list.append(list(n))
    # remove molecules from formula for glycosidic bonds
    atom_list_2 = []
    for i in range(len(atom_names)):
        n = np.array(atom_list[i])
        form_n = np.array([formulas['water'][i]] * len(masses.index))
        mol_n = np.array(molecule_numbers['dp'] - 1)
        form_mol_n = form_n * mol_n
        n = n + form_mol_n
        atom_list_2.append(list(n))
    # concatenate to build formulas
    for i in range(len(atom_names)):
        if i == 0:
            formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        else:
            formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
    # fix to remove atoms with zero
    formulas_final = formulas_final.str.replace("\D0", "")
    masses['formula'] = formulas_final

if "none" in modifications and pent_option == 0:
    if "benzoic_acid" in label:
        dp = masses.dp
        hex = masses.hex
        benzoic_acid = masses.benzoic_acid
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex,
                                         'benzoic_acid': benzoic_acid})
    else:
        dp = masses.dp
        hex = masses.hex
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex})
    molecules = list(molecule_numbers.drop('dp', axis=1).columns)
    atom_names = ["C", "H", "N", "O", "S", "P"]
    atom_list = []
    for i in range(len(atom_names)):
        n = np.array([0] * len(masses.index))
        for j in range(len(molecules)):
            form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers[molecules[j]])
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
        if "procainamide" in label:
            p = np.array([formulas['procainamide'][i]] * len(masses.index))
            n = n + p
        atom_list.append(list(n))
    # remove molecules from formula for glycosidic bonds
    atom_list_2 = []
    for i in range(len(atom_names)):
        n = np.array(atom_list[i])
        form_n = np.array([formulas['water'][i]] * len(masses.index))
        mol_n = np.array(molecule_numbers['dp'] - 1)
        form_mol_n = form_n * mol_n
        n = n + form_mol_n
        atom_list_2.append(list(n))
    # concatenate to build formulas
    for i in range(len(atom_names)):
        if i == 0:
            formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        else:
            formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
    # fix to remove atoms with zero
    formulas_final = formulas_final.str.replace("\D0", "")
    masses['formula'] = formulas_final

if "none" not in modifications and pent_option == 1:
    if unsaturated_option == 'y':
        modifications.append('unsaturated')
    if alditol_option == 'y':
        modifications.append('alditol')
    if dehydrated_option == 'y':
        modifications.append('dehydrated')
    if "benzoic_acid" in label:
        dp = masses.dp
        hex = masses.hex
        pent = masses.pent
        benzoic_acid = masses.benzoic_acid
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex,
                                         'pent': pent,
                                         'benzoic_acid': benzoic_acid})
        modification_numbers = masses[modifications]
        modification_numbers_array = np.array(modification_numbers)
        modification_numbers_array = modification_numbers_array.repeat(np.array(nmax_bc) + 1, axis=0)
        colNames = modification_numbers.columns
        modification_numbers = pd.DataFrame(modification_numbers_array)
        modification_numbers.columns = colNames
        molecule_numbers = pd.concat([molecule_numbers, modification_numbers], axis=1)
    else:
        dp = masses.dp
        hex = masses.hex
        pent = masses.pent
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex,
                                         'pent': pent})
        modification_numbers = masses[modifications]
        molecule_numbers = pd.concat([molecule_numbers, modification_numbers], axis=1)
    molecules = list(molecule_numbers.drop('dp', axis=1).columns)
    atom_names = ["C", "H", "N", "O", "S", "P"]
    atom_list = []
    for i in range(len(atom_names)):
        n = np.array([0] * len(masses.index))
        for j in range(len(molecules)):
            form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers[molecules[j]])
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
        if "procainamide" in label:
            p = np.array([formulas['procainamide'][i]] * len(masses.index))
            n = n + p
        atom_list.append(list(n))
    # remove molecules from formula for glycosidic bonds
    atom_list_2 = []
    for i in range(len(atom_names)):
        n = np.array(atom_list[i])
        form_n = np.array([formulas['water'][i]] * len(masses.index))
        mol_n = np.array(molecule_numbers['dp'] - 1)
        form_mol_n = form_n * mol_n
        n = n + form_mol_n
        atom_list_2.append(list(n))
    # concatenate to build formulas
    for i in range(len(atom_names)):
        if i == 0:
            formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        else:
            formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
    # fix to remove atoms with zero
    formulas_final = formulas_final.str.replace("\D0", "")
    masses['formula'] = formulas_final

if "none" not in modifications and pent_option == 0:
    if unsaturated_option == 'y':
        modifications.append('unsaturated')
    if alditol_option == 'y':
        modifications.append('alditol')
    if dehydrated_option == 'y':
        modifications.append('dehydrated')
    if "benzoic_acid" in label:
        dp = masses.dp
        hex = masses.hex
        benzoic_acid = masses.benzoic_acid
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex,
                                         'benzoic_acid': benzoic_acid})
        modification_numbers = masses[modifications]
        modification_numbers_array = np.array(modification_numbers)
        modification_numbers_array = modification_numbers_array.repeat(np.array(nmax_bc) + 1, axis=0)
        colNames = modification_numbers.columns
        modification_numbers = pd.DataFrame(modification_numbers_array)
        modification_numbers.columns = colNames
        molecule_numbers = pd.concat([molecule_numbers, modification_numbers], axis=1)
    else:
        dp = masses.dp
        hex = masses.hex
        molecule_numbers = pd.DataFrame({'dp': dp,
                                         'hex': hex})
        modification_numbers = masses[modifications]
        molecule_numbers = pd.concat([molecule_numbers, modification_numbers], axis=1)
    molecules = list(molecule_numbers.drop('dp', axis=1).columns)
    atom_names = ["C", "H", "N", "O", "S", "P"]
    atom_list = []
    for i in range(len(atom_names)):
        n = np.array([0] * len(masses.index))
        for j in range(len(molecules)):
            form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers[molecules[j]])
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
        if "procainamide" in label:
            p = np.array([formulas['procainamide'][i]] * len(masses.index))
            n = n + p
        atom_list.append(list(n))
    # remove molecules from formula for glycosidic bonds
    atom_list_2 = []
    for i in range(len(atom_names)):
        n = np.array(atom_list[i])
        form_n = np.array([formulas['water'][i]] * len(masses.index))
        mol_n = np.array(molecule_numbers['dp'] - 1)
        form_mol_n = form_n * mol_n
        n = n + form_mol_n
        atom_list_2.append(list(n))
    # concatenate to build formulas
    for i in range(len(atom_names)):
        if i == 0:
            formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        else:
            formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
    # fix to remove atoms with zero
    formulas_final = formulas_final.str.replace("\D0", "")
    masses['formula'] = formulas_final

gc.collect()

elapsed_time = time.time() - start_time
print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

# 4: FILTER BASED ON NUMBER OF POSSIBLE MODIFICATIONS
# ----------------------

print("\nstep #4: filtering based on number of modifications per monomer")
print("----------------------------------------------------------------\n")

if "none" not in modifications:
    if unsaturated_option == 'y':
        modifications.remove('unsaturated')
    if alditol_option == 'y':
        modifications.remove('alditol')
    if dehydrated_option == 'y':
        modifications.remove('dehydrated')
    masses['nmod'] = masses[modifications].sum(axis=1)
    masses['nmod_avg'] = masses.nmod / masses.dp
    masses = masses.drop(masses[masses.nmod_avg > nmod_max].index)

if 'anhydrobridge' in modifications and pent_option == 1:
    indexDelete = masses[masses.hex < masses.anhydrobridge].index
    masses.drop(indexDelete, inplace=True)
    masses = masses.reset_index()

gc.collect()
elapsed_time = time.time() - start_time
print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

# 5: ESTIMATE NUMBERS OF POSSIBLE ISOMERS
# ----------------------

print("\nstep #5: estimating number of possible isomers")
print("----------------------------------------------------------------\n")

# calculate for L or D forms
if len(LorD_isomers) == 2:
    masses['isomers'] = 2 ** masses.dp

if OH_stereo == 1:
    # hexose
    masses['OH'] = 0
    mask1 = (masses['hex'] != 0)
    masses_temp = masses[mask1]
    masses.loc[mask1, 'OH'] = ((masses_temp.hex - 1) * 2) + 3
    # pentose
    if pent_option == 1:
        mask2 = (masses['hex'] != 0) & (masses['pent'] != 0)
        masses_temp = masses[mask2]
        masses.loc[mask2, 'OH'] = masses_temp.OH + masses_temp.pent * 1
        mask3 = (masses['hex'] == 0)
        masses_temp = masses[mask3]
        masses.loc[mask3, 'OH'] = masses_temp.OH + masses_temp.pent * 1 + 2
    # modifications
    if "none" not in modifications:
        modifications_OHdiff = [val for i, val in enumerate(modifications) if val in isomers_OHdiff]
        for i in range(len(modifications_OHdiff)):
            masses['OH'] = masses.OH - masses[modifications_OHdiff[i]]
    if len(LorD_isomers) == 2:
        masses['isomers'] = masses.isomers * (2 ** masses.OH)
    elif len(LorD_isomers) == 1:
        masses['isomers'] = 2 ** masses.OH
    masses = masses.drop(columns=['OH'])

if bond_stereo == 1:
    if len(LorD_isomers) == 2 or OH_stereo == 1:
        masses['isomers'] = masses.isomers * (2 ** masses.dp)
    elif len(LorD_isomers) == 1 and OH_stereo == 0:
        masses['isomers'] = 2 ** masses.dp

gc.collect()
elapsed_time = time.time() - start_time
print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

# 6: CALCULATE M/Z VALUES
# ----------------------

print("\nstep #6: calculating m/z values of ions")
print("----------------------------------------------------------------\n")

if len(list(set(modifications).intersection(modifications_anionic))) >= 1:
    # create separate tables of sugars with (any) anionic modifications, and with (only) neutral modifications
    anionic_mod_used = list(set(modifications).intersection(modifications_anionic))
    masses_anionic = masses[masses['name'].str.contains('|'.join(anionic_mod_used))]
    masses_all = masses.merge(masses_anionic.drop_duplicates(), how='left', indicator=True)
    masses_neutral = masses_all[masses_all._merge == 'left_only']

    # calculate m/z values for neutral molecules
    if "neg" in ESI_mode:
        masses_neutral['[M-H]-'] = masses_neutral.mass - ion_mdiff['H'] + e_mdiff
        masses_neutral['[M+Cl]-'] = masses_neutral.mass + ion_mdiff['Cl'] + e_mdiff
        masses_neutral['[M+CHOO]-'] = masses_neutral.mass + ion_mdiff['CHOO'] + e_mdiff
        masses_neutral['[M-2H]-2'] = (masses_neutral.mass - 2 * ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses_neutral['[M+2Cl]-2'] = (masses_neutral.mass + 2 * ion_mdiff['Cl'] + 2 * e_mdiff) / 2
        masses_neutral['[M+2CHOO]-2'] = (masses_neutral.mass + 2 * ion_mdiff['CHOO'] + 2 * e_mdiff) / 2
        masses_neutral['[M+Cl-H]-2'] = (masses_neutral.mass + ion_mdiff['Cl'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses_neutral['[M+CHOO-H]-2'] = (masses_neutral.mass + ion_mdiff['CHOO'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses_neutral['[M+CHOO+Cl]-2'] = (masses_neutral.mass + ion_mdiff['CHOO'] + ion_mdiff['Cl'] + 2 * e_mdiff) / 2
    if "pos" in ESI_mode:
        masses_neutral['[M+H]+'] = masses_neutral.mass + ion_mdiff['H'] - e_mdiff
        masses_neutral['[M+Na]+'] = masses_neutral.mass + ion_mdiff['Na'] - e_mdiff
        masses_neutral['[M+NH4]+'] = masses_neutral.mass + ion_mdiff['NH4'] - e_mdiff
        masses_neutral['[M+K]+'] = masses_neutral.mass + ion_mdiff['K'] - e_mdiff

    # filter neutral molecules based on scan range
    # set values outside range to NaN
    # remove rows where all ions are outside range
    my_cols = list(masses_neutral.filter(like='[M', axis=1).columns)
    masses_neutral[my_cols] = masses_neutral[my_cols].where(masses_neutral[my_cols] >= scan_range[0])
    masses_neutral[my_cols] = masses_neutral[my_cols].where(masses_neutral[my_cols] <= scan_range[1])
    masses_neutral = masses_neutral.dropna(subset=my_cols, how='all')

    # calculate m/z values for anionic molecules
    if len(anionic_mod_used) > 1:
        masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].sum(axis=1)
        masses_anionic['nmod_anionic'] = masses_anionic.nmod_anionic.astype(int)
    elif len(anionic_mod_used) == 1:
        masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].astype(int)
    if "neg" in ESI_mode:
        ions = list(range(1, masses_anionic.nmod_anionic.max() + 1))
        ions = list("[M-" + pd.Series(ions).astype(str) + "H]-" + pd.Series(ions).astype(str))
        for i in range(len(ions)):
            masses_anionic[ions[i]] = (masses_anionic.mass - ion_mdiff['H'] * (i + 1) + e_mdiff * (i + 1)) / (i + 1)
            masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i + 1))
        masses_anionic = masses_anionic.rename({'[M-1H]-1': '[M-H]-'}, axis=1)
        masses_anionic['[M+Cl]-'] = masses_anionic.mass + ion_mdiff['Cl'] + e_mdiff
        masses_anionic['[M+CHOO]-'] = masses_anionic.mass + ion_mdiff['CHOO'] + e_mdiff
        masses_anionic['[M+2Cl]-2'] = (masses_anionic.mass + 2 * ion_mdiff['Cl'] + 2 * e_mdiff) / 2
        masses_anionic['[M+2CHOO]-2'] = (masses_anionic.mass + 2 * ion_mdiff['CHOO'] + 2 * e_mdiff) / 2
        masses_anionic['[M+Cl-H]-2'] = (masses_anionic.mass + ion_mdiff['Cl'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses_anionic['[M+CHOO-H]-2'] = (masses_anionic.mass + ion_mdiff['CHOO'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses_anionic['[M+CHOO+Cl]-2'] = (masses_anionic.mass + ion_mdiff['CHOO'] + ion_mdiff['Cl'] + 2 * e_mdiff) / 2
    if "pos" in ESI_mode:
        masses_anionic['[M+H]+'] = masses_anionic.mass + ion_mdiff['H'] - e_mdiff
        masses_anionic['[M+Na]+'] = masses_anionic.mass + ion_mdiff['Na'] - e_mdiff
        masses_anionic['[M+NH4]+'] = masses_anionic.mass + ion_mdiff['NH4'] - e_mdiff
        masses_anionic['[M+K]+'] = masses_anionic.mass + ion_mdiff['K'] - e_mdiff

    # filter anionic molecules based on scan range
    # set values outside range to NaN
    # remove rows where all ions are outside range
    my_cols = list(masses_anionic.filter(like='[M', axis=1).columns)
    masses_anionic[my_cols] = masses_anionic[my_cols].where(masses_anionic[my_cols] >= scan_range[0])
    masses_anionic[my_cols] = masses_anionic[my_cols].where(masses_anionic[my_cols] <= scan_range[1])
    masses_anionic = masses_anionic.dropna(subset=my_cols, how='all')

    # concatenate dataframes and format nicely to only have useful columns
    masses_final = pd.concat([masses_anionic, masses_neutral])
    if "benzoic_acid" in label:
        bad_cols = {'level_0',
                    'index',
                    'hex',
                    'alditol',
                    'pent',
                    'dehydrated',
                    'nmod',
                    'nmod_avg',
                    'nmod_anionic',
                    '_merge',
                    'benzoic_acid'}
    else:
        bad_cols = {'level_0',
                    'index',
                    'hex',
                    'pent',
                    'alditol',
                    'dehydrated',
                    'nmod',
                    'nmod_avg',
                    'nmod_anionic',
                    '_merge'}
    bad_cols.update(modifications_anionic)
    bad_cols.update(modifications_neutral)
    cols_del = list(set(masses_final.columns).intersection(bad_cols))
    masses_final = masses_final.drop(columns=cols_del)

if len(list(set(modifications).intersection(modifications_anionic))) == 0:
    # calculate m/z values for neutral molecules
    if "neg" in ESI_mode:
        masses['[M-H]-'] = masses.mass - ion_mdiff['H'] + e_mdiff
        masses['[M+Cl]-'] = masses.mass + ion_mdiff['Cl'] + e_mdiff
        masses['[M+CHOO]-'] = masses.mass + ion_mdiff['CHOO'] + e_mdiff
        masses['[M-2H]-2'] = (masses.mass - 2 * ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses['[M+2Cl]-2'] = (masses.mass + 2 * ion_mdiff['Cl'] + 2 * e_mdiff) / 2
        masses['[M+2CHOO]-2'] = (masses.mass + 2 * ion_mdiff['CHOO'] + 2 * e_mdiff) / 2
        masses['[M+Cl-H]-2'] = (masses.mass + ion_mdiff['Cl'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses['[M+CHOO-H]-2'] = (masses.mass + ion_mdiff['CHOO'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
        masses['[M+CHOO+Cl]-2'] = (masses.mass + ion_mdiff['CHOO'] + ion_mdiff['Cl'] + 2 * e_mdiff) / 2
    if "pos" in ESI_mode:
        masses['[M+H]+'] = masses.mass + ion_mdiff['H']
        masses['[M+Na]+'] = masses.mass + ion_mdiff['Na']
        masses['[M+NH4]+'] = masses.mass + ion_mdiff['NH4'] - e_mdiff
        masses['[M+K]+'] = masses.mass + ion_mdiff['K'] - e_mdiff


    # filter neutral molecules based on scan range
    # set values outside range to NaN
    # remove rows where all ions are outside range
    my_cols = list(masses.filter(like='[M', axis=1).columns)
    masses[my_cols] = masses[my_cols].where(masses[my_cols] >= scan_range[0])
    masses[my_cols] = masses[my_cols].where(masses[my_cols] <= scan_range[1])
    masses = masses.dropna(subset=my_cols, how='all')

    # format nicely to only have useful columns
    masses_final = masses
    if "benzoic_acid" in label:
        bad_cols = {'level_0',
                    'index',
                    'alditol',
                    'dehydrated',
                    'hex',
                    'pent',
                    'nmod',
                    'nmod_avg',
                    'nmod_anionic',
                    '_merge',
                    'benzoic_acid'}
    else:
        bad_cols = {'level_0',
                    'index',
                    'alditol',
                    'hex',
                    'dehydrated',
                    'pent',
                    'nmod',
                    'nmod_avg',
                    'nmod_anionic',
                    '_merge'}
    bad_cols.update(modifications_neutral)
    cols_del = list(set(masses_final.columns).intersection(bad_cols))
    masses_final = masses_final.drop(columns=cols_del)

elapsed_time = time.time() - start_time
print("finished. elapsed time = " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

# 7: WRITE OUTPUT TO FILE
# ----------------------

print("\nstep #7: writing output to file % s" % outfile)
print("----------------------------------------------------------------\n")


masses_final.to_csv(outfile,
                    sep='\t',
                    header=True,
                    index=False,
                    na_rep="NA")



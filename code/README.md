# Analysis walk-through: Search for a Higgs boson in the 4 muon final state.

## How to Use
Download the file with the fake data-sets from the book's website [terascale.de](https://www.terascale.de/e211693/index_eng.html#terascale_e211716). The files should be located under Chapter 11.
Open `analysis_functions.py` and make sure the path to the file `Histograms_fake.root` is correct on line 38.  
Run the script using the command:
```
$ source run.sh
```
If they do not exist already, the following folders are created:
- output/histograms
- output/figures

NOTE: The first time the script is run the test-statistic distributions for both the signal+background hypothesis and the background-only hypothesis is calculated using 10 000 toys for 4 different luminosities, so the first run may take some time.

If the files below have already been created (i.e this is not your first time running the script), you will be given an option to simply do a quick run where these existing distributions are used.

```
output/histograms/test_statistic_distribution_sfb1_sfs1_toys10000_bin1.root exists
output/histograms/test_statistic_distribution_sfb1.5_sfs1.5_toys10000_bin10.root exists
output/histograms/test_statistic_distribution_sfb2.0_sfs2.0_toys10000_bin10.root exists
output/histograms/test_statistic_distribution_sfb5.0_sfs5.0_toys10000_bin10.root exists

To use these distributions in the analysis type 0 (recommended)
To calculate new distributions press any key:
```
Simply type 0 (zero) to initiate a quick run, or some other key to run the full analysis again.

## Credits
This analysis is based on the chapter "Analysis Walkthrough" in the book Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods, Wiley-VCH (2013) and the code is based on the skeleton code written by Ivo van Vulpen and Aart Heijboer in 2013

Project: Exercises 11.1-11.3
File: Walkthrough_skeleton.C
Author: Ivo van Vulpen, Aart Heijboer
Version (date): 1.0 (23.06.2013)
Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
All rights reserved.

### TODO:
- Should implement in c++, loops in python slow.
- Make functions more general to be applicable on other datasets.
- Make statistics.py a module so it does not have to be in the working dir.
-

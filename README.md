This is an internal use tool for finding orthologous enhancers in various drosophila species.

## create_orth_datasets.pl <input.fa> <other_species_name> <output_file_prefix>

This program will use installed mapping files to go from dm3 (which is release 5) D.Mel coordinates to any of the installed other species. The code is based on inherited code. See the file "default_targets.txt" for some of the available other species.



## testcompare.py

This (unfinished) program is to help in judging the quality/accuracy of extracted orthologs.


## FetchOrthologousRegions.py

This program goes to the FlyBase website for a particular gene name, and cuts out the region surrounding that gene, +- a specified number of bp.
It will look on FlyBase for known orthologs of that gene, and get all orthologs (and the up/downstream sequence data.)

It's up to you, after that, to find the orthologous enhancer.

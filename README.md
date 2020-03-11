# lcMLkin_v2
Updated version of lcMLkin in python. Uses allele frequencies from an outside population (via a plink file), filters for LD on the fly, and allows some FST betweent the reference population and the target individuals.

To run, use:

> lcMLkin_optim_customAF_PL_LD_AWfst_MP.py input.vcf plink_file FST_value nb_threads

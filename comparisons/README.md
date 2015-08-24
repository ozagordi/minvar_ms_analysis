
## Watch out!

Both MinVar and VirVarSeq report the mutations with respect to consensus B and
the csv files do not list the percentage of remaining wild type. For the
scatterplots this is not a problem: a point with 30% mutation will always
mean a 70% wild type. For precision and recall the story changes if we want
to compare with Sanger based calls. Here if the peaks show both wt and mutant
it would be slightly unfair to discard the wt and implicitly assuming that this
is equivalent to the cases where only the mutant was detected.


# What to run
First expand mutations as explained above, then R script `precrec.R`

    [ozagordi@virologysrv04 comparisons]$ expand_mutations.py
    [ozagordi@virologysrv04 comparisons]$ Rscript precrec.R

This makes two figures: bar plots of F1 and, separately, precision/recall.
Further, it outputs a table with all numbers.

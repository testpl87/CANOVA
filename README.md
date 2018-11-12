# CANOVA
This is CANOVA tool for omics data analysis

It can perform constrained ANOVA to continuous expression data such as microarray gene expression data, RNA-seq so on.
This function uses ordered information of phenotypes.

compute_statistic_twoside calculates CANOVA statistic for two-sided test.
compute_statistic_oneside calculates CANOVA statistic for one-sided test.

dist_generate function is used for generating the null distribution. Distributions of CANOVA statistics depend on the numbers of sample sizes for the ordered groups.
Thus, input of dist_generate function is the number of samples for the groups. The more number of iteration gives the more accurate distribution.

test_with_one function is main function of CANOVA. test_with_one function compute the test statistic and p-value of the test.
It also gives the histogram to show the null distribution and the location of the computed statistic.

test_bulk gives the set of test statsitics and p-values.

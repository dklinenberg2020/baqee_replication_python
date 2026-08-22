README file for replication

The code in this folder replicates the results from
"What if? The Economic Effects for Germany of a Stop of Energy Imports from Russia"
by Rüdiger Bachmann, David Baqaee, Christian Bayer, Moritz Kuhn, Andreas Löschel, Benjamin Moll, Andreas Peichl, Karen Pittel, and Moritz Schularick
The policy brief is available at 
https://www.econtribute.de/RePEc/ajk/ajkpbs/ECONtribute_PB_028_2022.pdf

We thank to Lennard Schlattmann (University of Bonn) and Sihwan Yang (UCLA) for outstanding research assistance. 

The replication material comprises two parts. Part I replicates the model results (Section 2). Part II the empirical results on the distribution of energy expenditures (Section 3).

*Part I : Model results (section 2)*
(1) To replicate the results from the Baqaee-Farhi model in column 1 of Table 2 in the main text and Appendix Table 1, go to folder "baqaee_farhi_model". It contains a separate readme file. You will need MATLAB to run the model code.
- Open and run file "main_dlogW_rev_bigshocks_EU_Russian_v2.m"
- The current setting for the size of shocks is set high enough to ensure that the volume of trade between EU and Russia goes very close to zero. 

(2) To replicate the results in columns 2 and 3 of Table 2 in the main text as well as all figures in Appendix A ("Appendix to Section 2"), run the file "elasticity.m". You will need MATLAB to run the model code.


*Part II : Empirical results on distribution of expenditure/income shares (section 3)*

consumption_energy.do
(1) You need access to the SUF of the EVS data from the German Statistical Office (see do file for details)
(2) The working directory for the code on line 20 must be adjusted and the data needs to be put in the data folder of the working directory (line 24).
(3) Running the code will reproduce the figures on expenditure and income shares from the paper.
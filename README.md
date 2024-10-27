This is the code of the model reported in Stem cell homeostasis in the root of Arabidopsis involves cell type specific complex formation of key transcription factors, Vivien I. Strotmann, Monica L. Garci­a-Gomez, and Yvonne Stahl. https://www.biorxiv.org/content/10.1101/2024.04.26.591257v1.abstract

To simulate the formation of protein complexes in the root stem cell niche, first run 1-ParameterEstimation.R, and then run either 2-Parameter_test.R, 3-RSCN_simulations.R, or 4-RSCN_simulations_controls.R. Below is a description of what each one does: 

1-ParameterEstimation.R: Estimation of assoaciation and dissociation parameters that explain relative binding rates determined experimentally
2-Parameter_test.R: Many parameter sets are recovered in (1), and here we show that the protein complexes formed using the protein levels of the cells of the root stem cell niche are recovered independently of the specific parameters used.
3-RSCN_simulations.R: Wt and PLT3dPrD simulations using the association / dissociation rates found in (1) and the protein levels determined experimentally for the stele initials (SI), quiescent center (QC), columella stem cells (CSC) and columella cells (CC).
4-RSCN_simulations_controls.R: Control simulations changing the association and dissociation rates used, and the protein levels simulated.

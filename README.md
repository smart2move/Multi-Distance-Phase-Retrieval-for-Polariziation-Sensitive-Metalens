# Multi-Distance-Phase-Retrieval-for-Polariziation-Sensitive-Metalens
This project contains files for polarization sensitive metalens design and simulation in Ansys lumerical. Furthermore, It contains the main algorithm for MDPR in MATLAB. Finally some other files for data export from ImageJ files from the CCD to matlab, and other files used for the project

Short explaination of all the files:

Multi_distance_phase_retrieval.pdf: Paper on algorithm with deeper explaination

unit_cell.lsf: code for Ansys lumerical unit cell simulation and parametersweep
unit_cell_git.7z: Ansys project inviroment for unit cell simulation and parametersweep
full_lens.lsf: code for Ansys lumerical full lens simulation
full_lens_git.zip: Ansys project inviroment for full lens simulation
stich_nearfield_ZOS_R100um.lsf: Near field stitching for large radius simation of metalens
gds_export.lsf: GDS export for metalens design

Phase_unwrapping_and_interpolation.mat: 2D phase unwrapping and interpolation for data from unitcell sweep 

MDPR_operator.mat: Operator for MDPR algorithm
multidistance_phase_retrieval_algorithm.mat: function for MDPR algorithm to use with MDPR_operator

CCD_to_MATLAB_matrix_converter.mat: Converting images taken by CCD to matlab using imagaj. MATLAB imageJ extention needed
Pixelsize_increaser.mat: increasing pixel size of matlab matrix of CCD image

max_unit_cell.mat: calculator for max unit cell in metalens design
measurement_requirements.mat: calculator for measurement requirements MDPR





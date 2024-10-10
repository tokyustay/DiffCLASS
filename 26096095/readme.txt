Reflection matrix microscopy and CLASS sample data and algorithms

This software package includes essential algorithms for reflection matrix microscope (RMM) and the 
CLASS algorithm. RMM is an advanced deep-tissue imaging technique that measures the reflected 
electric field of thick scattering samples while varying the illumination channels. The package contains 
detailed algorithm implementations, experimental and numerical data to promote the widespread adoption 
of RMM in various deep-tissue imaging applications.

The package comprises RMM and CLASS functions, sample data, and example scripts that demonstrate 
the processing of the sample data. All functions and scripts in the package were developed and tested 
with Matlab R2021a.

All codes and datasets were created by Sungsam Kang from Super-depth imaging group at Korea 
University. (https://www.bioimaging.korea.ac.kr/)

Sungsam Kang
2024. 06. 25


References
*	CASS microscopy: S. Kang, et. al., Nature Photonics 9 (4), 253-258 (2015)
               https://www.nature.com/articles/nphoton.2015.24
*	CLASS microscopy: S. Kang, et. al., Nature Communications 8, 2157 (2017)
               https://www.nature.com/articles/s41467-017-02117-8

Demo scripts
*	Example_numerical_data.m: Demonstrates RM construction, image reconstruction, and CLASS 
functions from the simulation dataset ¡®numerical_data_RMM_simulation.mat¡¯.
*	Example_experimental_data.m: Demonstrate RM construction, image reconstruction, and CLASS 
functions from the experimental dataset ¡®Experimental_data_USAF.mat¡¯.
*	Example_3D_reconstruction.m: Demonstrate 3D image reconstruction from the simulation 
dataset ¡®numerical_data_3D_beads.mat¡¯.
*	Example_RMM_simulation.m: Generate RM simulation data.

Functions (...\RMM_functions)
*	RM_get_image.m: Reconstruct target images from an RM using aperture synthesis. 
*	RM_get_dk_logical_index.m: Retrieves indices for delta-k space mapping for logical indexing.
*	RM_fft.m: Performs matrix Fourier transform from position basis to spatial frequency basis.
*	RM_ifft.m: Performs inverse Fourier transform from spatial frequency basis to position basis.
*	RM_3D_recon.m: Reconstructs 3D images from RM.
*	Indexing_array_mt.mexw64: Mex function for multi-thread logical indexing.
*	Indexing_array_mt.cpp: C++ source code for ¡®Indexing_array_mt.mexw64¡¯. 
Compile this file if mex function do not work properly using following command:
¡®mex -R2018a indexing_array_mt.cpp¡¯
*	CLASS.m: The CLASS function that quantifies and corrects input and output aberrations from 
RM.
*	RM_ab_correction.m: Corrects input and output aberrations from RM using the CLASS result.

Datasets (...\RMM_sample_data)
*	Experimental_data_USAF.mat: Experimentally measured RM data of a USAF target under a 
1.5mm thick glass plate. Detailed description of the data is included in the script 
¡®Example_experimental_data.m¡¯.
*	numerical_data_RMM_simulation.mat: Numerically generated RM data of the ¡®cameraman.tif¡¯ 
image with and without aberrations. Detailed description of the data is included in the script 
¡®Example_numerical_data.m¡¯.
*	numerical_data_3D_beads.mat: Numerically generated RM data of point particles distributed 
within a 3D volume. Detailed description of the data is included in the script 
¡®Example_3D_reconstruction.m¡¯.
*	numerical_ground_truth_aberrations: Numerically generated input and output aberrations for RM 
simulations.




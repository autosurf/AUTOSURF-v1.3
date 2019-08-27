
-----------------------------------------------------------------------------------
                                 A U T O S U R F      
                                 
-----------------------------------------------------------------------------------

   AUTOSURF Program Suite

   Copyright (c) 2019 Missouri University of Science and Technology

   This file is part of AUTOSURF.
   
-----------------------------------------------------------------------------------

 AUTOSURF is a freely distributed suite of codes for the automated construction of potential 
 energy surfaces (PES) for vdW systems. The fitting algorithms implemented in the code are 
 based on the L-IMLS methodology, and have many advanced features such as options for data-
 point placement, flexibility to include gradients in the fit, iterative refinement, and 
 symmetry recognition. The code completely automates all the steps and procedures that go 
 into fitting various classes of PESs and interfaces to popular electronic structure codes 
 such as MOLPRO and GAUSSIAN. 

 For non-commercial purposes, there is no charge to obtain AUTOSURF package for academic 
 users (research institutes, universities, individuals), the user must simply sign a License 
 Agreement. The License form can be obtained from the authors by email. The link for 
 downloading the most recent version of AUTOSURF will be sent to the email address used in 
 the License form.

 The package (v1.3) is composed of three main programs: AUTOSURF-ABI, AUTOSURF-PES, and AUTOSURF-PLOT:
 * AUTOSURF-ABI performs guided surveys of the PES (various cuts), facilitating 
 the benchmarking of electronic structure methods, and the development of composite 
 schemes such as complete basis set (CBS) extrapolation.
 * AUTOSURF-PES carries out the automated construction of the PES to a user-specified     
 accuracy target in a fairly black-box fashion: starting with a sparse 
 set of initial ab initio seed points, the program grows a fitted PES over 
 predefined ranges of energy and coordinates until the desired level of precision 
 is reached. 
 * AUTOSURF-PLOT permits arbitrary evaluations of the PES, the generation of 
 plots of 1D or 2D cuts of the surface (with optional relaxation) for any of the 
 internal variables, and also to perform a variety of fitting error analyses in 
 specified energy and coordinate ranges.

-----------------------------------------------------------------------------------
 In order to install and run AUTOSURF suite on a Linux system (Windows and OS X 
 are currently not supported) you will need the following:
 ========================================================================
 * A Fortran90 compatible compiler.
 * MPICH or MPICH2.
 * The Fortran LAPACK Library.
 * Electronic structure code (MOLPRO and/or Gaussian).

 The standard AUTOSURF (v1.3) package is distributed as a compressed archive file 
 named autosurf-v1.3.tar.gz, which includes the Fortran source codes of the three 
 programs conforming AUTOSURF. Once the package have been downloaded, the user 
 should simply unpack the file in the desired location and execute the corresponding 
 Makefiles to generate the binaries (the Makefiles has to be modified according to 
 the user's system requirements). The former is just a brief description of the 
 installation procedure. Users should read AUTOSURF documentation for more detailed 
 information.

-----------------------------------------------------------------------------------
***********************************************************************************
******                Welcome to the AUTOSURF Program Suite !                ******         
***********************************************************************************

-----------------------------------------------------------------------------------
   Please refer to the following files for further information:
-----------------------------------------------------------------------------------
   * AUTHORS.tst:   The list of AUTOSURF developers.
   * COPYRIGHT.txt: The license under which AUTOSURF can be used.
   * INSTALL.txt:   Tnstructions to compile and install the program.
-----------------------------------------------------------------------------------


   Directory structure:
========================================================================
   * SOURCE:        The (Fortran, 90 & 77) source files + Makefiles.
   * BIN:           The compiled executables.
   * DOC:           All the documentation related to AUTOSURF (manual & tutorial).
   * EXAMPLES:      Documented examples INPUTS/OUTPUTS for benchmark systems that 
                    can be used as a test suite to check the installation.
========================================================================

-----------------------------------------------------------------------------------
   For further information, visit the web page:
========================================================================
   https://equintas1.wixsite.com/autosurf-pes-rigid4d



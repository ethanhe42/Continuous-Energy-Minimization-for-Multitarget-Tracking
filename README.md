###################################################################
#                                                                 #
#     Continuous Energy Minimization For Multitarget Tracking     #
#         Anton Milan and Stefan Roth and Konrad Schindler        #
#                                                                 #
#                Copyright 2010-2014 Anton Milan                  #
#                                                                 #
###################################################################



About
=====

This software implements our approach to multi-target tracking
using continuous energy minimization [1,2,3].

The additional package
 - minimize.m
 
is released under a different license and is included for your convenience.



==========================================================================
DISCLAIMER:
This software has been refactored for the sake of simplifying the
implementation. Therefore, the results produced by the code may differ
from those presented in the papers [1,2,3].
==========================================================================


To use this software, you should cite the following in any resulting publication:
	
    [1]	Continuous Energy Minimization for Multitarget Tracking [(pdf)](http://www.milanton.de/files/pami2014/pami2014-anton.pdf)
	A. Milan, S. Roth, and K. Schindler.
	In IEEE TPAMI 36(1), 2014
	
    [2] Multi-target Tracking by Continuous Energy Minimization [(pdf)](http://www.milanton.de/files/cvpr2011/cvpr2011-anton.pdf)
        A. Andriyenko and K. Schindler. 
        In CVPR, Colorado Springs, USA, June 2011

    [3] An Analytical Formulation of Global Occlusion Reasoning for Multi-Target Tracking [(pdf)](http://www.milanton.de/files/vs2011/vs2011-anton.pdf)
        A. Andriyenko, S. Roth, and K. Schindler. 
        In IEEE International Workshop on Visual Surveillance (in conjunction with ICCV), Barcelona, Spain, November 2011 


Installation
============
This section describes how to get dctracking running under Linux.
Open a terminal window.

Get the code

    hg clone https://bitbucket.org/amilan/dctracking
    cd dctracking
    

Start MATLAB and run compileMex.m to build the utilities binaries.


Running
=======

run cemTrackerDemo.m



CHANGES
	May 25, 2014	Included PAMI code
	May 25, 2012	Initial public release
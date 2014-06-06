Continuous Energy Minimization for Multitarget Tracking
=======================================================
Anton Milan and Stefan Roth and Konrad Schindler
------------------------------------------------




About
=====

This software implements our approach to multi-target tracking
using continuous energy minimization [1,2,3].

The additional package

    minimize.m
 
is released under a different license and is included for your convenience.




Important:
==========
This software has been refactored for the sake of simplifying the
implementation. Therefore, the results produced by the code may differ
from those presented in the papers [1,2,3].


References
==========
To use this software, you should cite a subset of these publications:

    [1] Continuous Energy Minimization for Multitarget Tracking
        A. Milan, S. Roth, and K. Schindler.
        In IEEE TPAMI 36(1), 2014
        
    [2] Multi-target Tracking by Continuous Energy Minimization
        A. Andriyenko and K. Schindler. 
        In CVPR, Colorado Springs, USA, June 2011

    [3] An Analytical Formulation of Global Occlusion Reasoning for Multi-Target Tracking
        A. Andriyenko, S. Roth, and K. Schindler. 
        In IEEE International Workshop on Visual Surveillance (in conjunction with ICCV), Barcelona, Spain, November 2011 

All papers are available [here](http://research.milanton.net)


Installation
============
This section describes how to get dctracking running under Linux.
Open a terminal window.

Get the code

    hg clone https://bitbucket.org/amilan/contracking
    cd contracking
    

Start MATLAB and run compileMex.m to build the utilities binaries.


Running
=======

run cemTrackerDemo.m



CHANGES

	Jun 06, 2014	Included Dynamic Programming [Pirsiavash et al., CVPR '11] as initialization
	May 25, 2014	Included PAMI code
	May 25, 2012	Initial public release
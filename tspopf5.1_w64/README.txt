=======================================================
 TSPOPF - High Performance AC OPF Solvers for MATPOWER
=======================================================

Version:    5.1

Home Page:  http://www.pserc.cornell.edu/tspopf/

Authors:    Hongye Wang <hw41@cornell.edu>
            Ray Zimmerman <rz10@cornell.edu>

            Thu, Mar 12, 2015

$Id: README.txt 2645 2015-03-11 20:53:01Z ray $
Copyright (c) 2007-2014 by Power System Engineering Research Center (PSERC)
See http://www.pserc.cornell.edu/tspopf/ for more info.


--------------
 INTRODUCTION
--------------

TSPOPF is a collection of three high performance AC optimal power flow
solvers for use with MATPOWER <http://www.pserc.cornell.edu/matpower/>,
a MATLAB(R) power system simulation package. The three solvers are:

  PDIPM   - primal/dual interior point method
  SCPDIPM - step-controlled primal/dual interior point method
  TRALM   - trust region based augmented Lagrangian method
   
The algorithms are described in:

    H. Wang, C. E. Murillo-S‡nchez, R. D. Zimmerman, R. J. Thomas,
    "On Computational Issues of Market-Based Optimal Power Flow",
    IEEE Transactions on Power Systems, Vol. 22, No. 3, Aug. 2007,
    pp. 1185-1193.

The PDIPM in particular is significantly faster for large systems than
any previous MATPOWER OPF solver, including MINOPF. When TSPOPF is
installed, the PDIPM solver becomes the default optimal power flow
solver for MATPOWER. Additional options for TSPOPF can be set using
mpoption (see 'help mpoption' for details).

--------------
 TERMS OF USE
--------------

- TSPOPF is free of charge. Anyone may use it.
- We make no warranties, express or implied. Specifically, we make no
  guarantees regarding the correctness TSPOPF's code or its fitness for any
  particular purpose.
- Any publications derived from the use of MATPOWER (and TSPOPF) must
  acknowledge MATPOWER <http://www.pserc.cornell.edu/matpower/>.
- Anyone may modify TSPOPF for their own use as long as the original
  copyright notices remain in place.
- TSPOPF may not be redistributed without written permission.
- Modified versions of TSPOPF, or works derived from TSPOPF, may not be
  distributed without written permission.


---------------------
 SYSTEM REQUIREMENTS
---------------------

    - MATLAB(R) version 7 or later
        - TRALM requires MATLAB 7.3 or later
        - Windows builds may require MATLAB 7.3 or later
    - MATPOWER 5.0 or later <http://www.pserc.cornell.edu/matpower/>
      (for TSPOPF 5.x)


--------------
 INSTALLATION
--------------

    1.  Unzip the downloaded file. 
    2.  Place the files in a location on your MATLAB path. 


---------------------------
 WHAT'S NEW IN VERSION 5.1
---------------------------

Below is a summary of the changes since version 4.1 of TSPOPF. See the
CHANGES file for all the gory details.

* New features:
  - Compatibility with MATPOWER 5.x and new options architecture.
  - Minor modification to make it more robust in cases where some
    parameters are Inf (e.g. limits for dummy gen representing DC line).


-----------
 FILE LIST
-----------

  CHANGES              A detailed change history, with Unix line endings
  CHANGES.txt          A detailed change history, with DOS line endings (for
                       Windows users)
  mpoption_info_tspopf.m  Returns MATPOWER option info for TSPOPF
  pdipmopf.mexglx      PDIPMOPF MEX file for 32-bit Linux
  pdipmopf.mexa64      PDIPMOPF MEX file for 64-bit Linux
  pdipmopf.mexmac      PDIPMOPF MEX file for Mac OS X (PPC)
  pdipmopf.mexmaci     PDIPMOPF MEX file for Mac OS X (Intel 32-bit)
  pdipmopf.mexmaci64   PDIPMOPF MEX file for Mac OS X (Intel 64-bit)
  pdipmopf.mexw32      PDIPMOPF MEX file for 32-bit Windows
  pdipmopf.mexw64      PDIPMOPF MEX file for 64-bit Windows
  pdipmopfver.m        Prints or returns version information for PDIPMOPF
  README               This file, but with Unix line endings
  README.txt           This file
  scpdipmopf.mexglx    SCPDIPMOPF MEX file 32-bit Linux
  scpdipmopf.mexa64    SCPDIPMOPF MEX file 64-bit Linux
  scpdipmopf.mexmac    SCPDIPMOPF MEX file for Mac OS X (PPC)
  scpdipmopf.mexmaci   SCPDIPMOPF MEX file Mac OS X (Intel 32-bit)
  scpdipmopf.mexmaci64 SCPDIPMOPF MEX file Mac OS X (Intel 64-bit)
  scpdipmopf.mexw32    SCPDIPMOPF MEX file for 32-bit Windows
  scpdipmopf.mexw64    SCPDIPMOPF MEX file for 64-bit Windows
  scpdipmopfver.m      Prints or returns version information for SCPDIPMOPF
  tralmopf.mexglx      TRALMOPF MEX file 32-bit Linux
  tralmopf.mexa64      TRALMOPF MEX file 64-bit Linux
  tralmopf.mexmac      TRALMOPF MEX file for Mac OS X (PPC)
  tralmopf.mexmaci     TRALMOPF MEX file Mac OS X (Intel 32-bit)
  tralmopf.mexmaci64   TRALMOPF MEX file Mac OS X (Intel 64-bit)
  tralmopf.mexw32      TRALMOPF MEX file for 32-bit Windows
  tralmopf.mexw64      TRALMOPF MEX file for 64-bit Windows
  tralmopfver.m        Prints or returns version information for TRALMOPF
  tspopf_solver.m      OPF solver invoked by MATPOWER's opf() which calls
                       the appropriate MEX file


---------
 SUPPORT
---------

Questions about TSPOPF can be addressed to the MATPOWER mailing list
(see <http://www.pserc.cornell.edu/matpower/#mailinglist>).

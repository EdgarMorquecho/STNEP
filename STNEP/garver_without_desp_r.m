function mpc = garver_without_desp_r

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenario without dispatchable generation and allowing shunt compensation in the load nodes
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Joe H. Chow's book, p. 70.

%   MATPOWER
%   $Id: case9.m,v 1.11 2010/03/10 18:08:14 ray Exp $

%% MATPOWER Case Format : Version 2
mpc.version = '2';
%% system MVA base
mpc.baseMVA = 100;
%% bus data
%  bus_i type Pd   Qd      Gs   Bs area	Vm	Va baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	  80	16	0	0	1	1	 0	230	    1	1.05	0.95; %80  16
	2	1	  240	48	0	0	1	1	 0	230	    1	1.05	0.95; %240  48
	3	2	  40 	8	0	0	1	1    0	230  	1	1.05	0.95; %40  8  
	4	1	  160	32	0	0	1	1	 0	230	    1	1.05	0.95; %160  32 
	5	1	  240	48	0	0	1	1	 0	230	    1	1.05	0.95; %240  48
	6	1	   0	0	0	0	1	1    0	230	    1	1.05	0.95; %0    0
];
%% generator data
%	bus	Pg	Qg	    Qmax	Qmin	Vg	  mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	 50	 40 	48   -10 	 1.05	 100	  1	    150	    0;   %50
    2     0  0      0      0     1.05    100      1     1000    0;   %ficticio   
    2     0  0      0      0     1.05    100      1     0       0;   %ficticio   
    3	165	 100	101  -10 	 1.05	 100	  1	    165	   165;   %165
    4     0  0      0      0     1.05    100      1     1000    0;  %ficticio
    4     0  0      0      0     1.05    100      1     0       0;  %ficticio
    5     0  0      0      0     1.05    100      1     1000    0;   %ficticio
    5     0  0      0      0     1.05    100      1     0       0;   %ficticio
	6	545	150  	183  -10 	 1.05	 100	  1     545	   545;   %545
	];


%% branch data
%  fbus	tbus	r	  x  	b  rateA rateB rateC	ratio	angle	status
mpc.branch = [
	1	 2	   0.040 	 0.40  0.00  120   120	 120      1     0	  1 -30	30;
	1	 4	   0.060 	 0.60  0.00	 100   100	 100      1     0	  1 -30	30;
	1	 5     0.020 	 0.20  0.00	 120   120	 120   	  1     0	  1 -30	30;
	2	 3	   0.020     0.20  0.00	 120   120	 120      1     0	  1 -30	30;
	2	 4	   0.040     0.40  0.00	 120   120	 120      1     0     1 -30	30;
	3	 5	   0.020     0.20  0.00	 120   120	 120      1     0	  1 -30	30;
    ]; 

mpc.ne_branch = [
	1	 2	   0.040 	 0.40  0.00  120.0   0.0	 0.0    0.0     0.0    0     -30	30    40;
  	1	 3	   0.038	 0.38  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30	30    38;
  	1	 4	   0.060 	 0.60  0.00	 100.0   0.0	 0.0    0.0     0.0    0     -30	30    60;
 	1	 5     0.020 	 0.20  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30	30    20;
 	1	 6	   0.068     0.68  0.00	  90.0   0.0	 0.0    0.0     0.0    0     -30	30    68;
	2	 3	   0.020     0.20  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30 	30	  20;
 	2	 4	   0.040     0.40  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30	30    40;   
	2	 5	   0.031     0.31  0.00	 120.0   0.0     0.0    0.0     0.0    0     -30	30    31;
	2	 6	   0.030     0.30  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30	30    30;    
	3	 4	   0.059     0.59  0.00	 102.0   0.0	 0.0    0.0     0.0    0     -30	30    59;
	3	 5	   0.020     0.20  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30	30    20;        
    3	 6	   0.048     0.48  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30	30    48;
    4	 5	   0.063     0.63  0.00	  95.0   0.0	 0.0    0.0     0.0    0     -30	30    63;        
    4	 6	   0.030     0.30  0.00	 120.0   0.0	 0.0    0.0     0.0    0     -30	30    30;
    5	 6	   0.061     0.61  0.00	  98.0   0.0	 0.0    0.0     0.0    0     -30	30    61;
	];  
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	cn	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0 0 0 ;  %%1
	2	0	0	3	0 1000 0;%%2
	2	0	0	3   0 0 0;%%2
	2	0	0	3   0 0 0;%%3
    2   0   0   3   0 1000 0;%%4
	2	0	0	3   0 0 0;%%4
    2   0   0   3   0 1000 0;%%5
    2   0   0   3   0 0 0;%%5
    2   0   0   3   0 0 0;%%6
    2	0	0	3	0 0 0;%%1
	2	0	0	3	0 0 0;%%2
	2	0	0	3	0 0 0;%%2
	2	0	0	3   0 0 0;%%3
    2   0   0   3   0 0 0;%%4
    2   0   0   3   0 0 0;%%4
    2   0   0   3   0 0 0;%%5
    2   0   0   3   0 0 0;%%5
    2   0   0   3   0 0 0;%%6
];
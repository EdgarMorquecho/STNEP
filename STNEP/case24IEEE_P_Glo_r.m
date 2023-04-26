function mpc = case24IEEE_P_Glo_r

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenario with dispatchable generation and without allowing shunt compensation in the load nodes
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Joe H. Chow's book, p. 70.

%   MATPOWER
%   $Id: case9.m,v 1.11 2010/03/10 18:08:14 ray Exp $

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
        1	3	324	66	0	0	1	0.955	10.67	138	1	1.05	0.95;
        2	2	291	60	0	0	1	0.971	11.22	138	1	1.05	0.95;
        3	2	540	111	0	0	1	0.968	11.56	138	1	1.05	0.95;
        4	1	222	45	0	0	1	0.998	15.28	138	1	1.05	0.95;
        5	1	213	42	0	0	1	1.002	15.73	138	1	1.05	0.95;
        6	2	408	84	0	0	1	0.99	13	138	1	1.05	0.95;
        7	2	375	75	0	0	1	0.989	12.56	138	1	1.05	0.95;
        8	1	513	105	0	0	1	1.015	20.77	345	1	1.05	0.95;
        9	2	525	108	0	0	1	1.043	28.02	345	1	1.05	0.95;
        10	1	585	120	0	0	1	1.05	35.61	345	1	1.05	0.95;
        11	1	0	0	0	0	1	0.985	12.72	138	1	1.05	0.95;
        12	1	0	0	0	0	1	0.99	12.2	138	1	1.05	0.95;
        13	2	795	162	0	0	1	0.968	11.35	138	1	1.05	0.95;
        14	2	582	117	0	0	1	0.984	11.5	138	1	1.05	0.95;
        15	2	951	192	0	0	1	0.97	11.23	138	1	1.05	0.95;
        16	2	300	60	0	0	1	0.984	11.91	138	1	1.05	0.95;
        17	1	0	0	0	0	1	0.995	13.74	138	1	1.05	0.95;
        18	2	999	204	0	0	1	0.973	11.53	138	1	1.05	0.95;
        19	1	543	111	0	0	1	0.963	11.05	138	1	1.05	0.95;
        20	1	384	78	0	0	1	0.958	11.93	138	1	1.05	0.95;
        21	2	0	0	0	0	1	0.959	13.52	138	1	1.05	0.95;
        22	2	0	0	0	0	1	0.97	16.08	138	1	1.05	0.95;
        23	2	0	0	0	0	1	1	    21	    138	1	1.05	0.95;
        24	1	0	0	0	0	1	0.992	20.89	138	1	1.05	0.95;
];

%% generator data

mpc.gen = [
1	0	0	240	-150	1	100	1	576	0	0	0	0	0	0	0	0	0	0	0	0;
2	0	0	240	-150	1	100	1	576	0	0	0	0	0	0	0	0	0	0	0	0;
3	0	0	9999	-9999	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
4   0  0    000         0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
4   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
5   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
5   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
6	0	0	0	-300	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
7	0	0	540	0	1	100	1	900	0	0	0	0	0	0	0	0	0	0	0	0;
8   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
8   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
9	0	0	9999	-9999	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
10   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
10   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
11   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
11   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
12   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
12   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
13	0	0	720	0	1	100	1	1773	0	0	0	0	0	0	0	0	0	0	0	0;
14	0	0	600	-150	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
15	0	0	330	-150	1	100	1	645	0	0	0	0	0	0	0	0	0	0	0	0;
16	0	0	240	-150	1	100	1	465	0	0	0	0	0	0	0	0	0	0	0	0;
17   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
17   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
18	0	0	600	-150	1	100	1	1200	0	0	0	0	0	0	0	0	0	0	0	0;
19   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
19   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
20   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
20   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;
21	0	0	600	-150	1	100	1	1200	0	0	0	0	0	0	0	0	0	0	0	0;
22	0	0	288	-180	1	100	1	900	0	0	0	0	0	0	0	0	0	0	0	0;
23	0	0	930	-375	1	100	1	1980	0	0	0	0	0	0	0	0	0	0	0	0;
24   0  0    000        0   1   100 1   1000   0   0   0   0   0   0   0   0   0   0   0   0;
24   0  0       0     -00   1   100 1   0   0   0   0   0   0   0   0   0   0   0   0   0;    
];



%% branch data
%	fbus	tbus  r	       x	   b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
1	2	0.0026	0.0139	0.4611	200	0	0	0	0	1	-360	360
1	3	0.0546	0.2112	0.0572	220	0	0	0	0	1	-360	360
1	5	0.0218	0.0845	0.0229	220	0	0	0	0	1	-360	360
2	4	0.0328	0.1267	0.0343	220	0	0	0	0	1	-360	360
2	6	0.0497	0.192	0.052	220	0	0	0	0	1	-360	360
3	9	0.0308	0.119	0.0322	220	0	0	0	0	1	-360	360
3	24	0.0023	0.0839	0	600	0	0	0	0	1	-360	360
4	9	0.0268	0.1037	0.0281	220	0	0	0	0	1	-360	360
5	10	0.0228	0.0883	0.0239	220	0	0	0	0	1	-360	360
6  10   0.0139  0.0605  2.4590  200  0 0   1   0   1 -360	360
7   8   0.0159  0.0614  0.0166  220  0 0   1   0   1 -360	360
8	9	0.0427	0.1651	0.0447	220	0	0	0	0	1	-360	360
8	10	0.0427	0.1651	0.0447	220	0	0	0	0	1	-360	360
9	11	0.0023	0.0839	0	600	0	0	0	0	1	-360	360
9	12	0.0023	0.0839	0	600	0	0	0	0	1	-360	360
10	11	0.0023	0.0839	0	600	0	0	0	0	1	-360	360
10	12	0.0023	0.0839	0	600	0	0	0	0	1	-360	360
11	13	0.0061	0.0476	0.0999	625	0	0	0	0	1	-360	360
11	14	0.0054	0.0418	0.0879	625	0	0	0	0	1	-360	360
12	13	0.0061	0.0476	0.0999	625	0	0	0	0	1	-360	360
12	23	0.0124	0.0966	0.203	625	0	0	0	0	1	-360	360
13	23	0.0111	0.0865	0.1818	625	0	0	0	0	1	-360	360
14	16	0.005	0.0389	0.0818	625	0	0	0	0	1	-360	360
15	16	0.0022	0.0173	0.0364	625	0	0	0	0	1	-360	360
15	21	0.0063	0.049	0.1030	625	0	0	0	0	1	-360	360
15	21	0.0063	0.049	0.1030	625	0	0	0	0	1	-360	360
15	24	0.0067	0.0519	0.1091	625	0	0	0	0	1	-360	360
16	17	0.0033	0.0259	0.0545	625	0	0	0	0	1	-360	360
16	19	0.003	0.0231	0.0485	625	0	0	0	0	1	-360	360
17	18	0.0018	0.0144	0.0303	625	0	0	0	0	1	-360	360
17	22	0.0135	0.1053	0.2212	625	0	0	0	0	1	-360	360
18	21	0.0033	0.0259	0.0545	625	0	0	0	0	1	-360	360
18	21	0.0033	0.0259	0.0545	625	0	0	0	0	1	-360	360
19	20	0.0051	0.0396	0.0833	625	0	0	0	0	1	-360	360
19	20	0.0051	0.0396	0.0833	625	0	0	0	0	1	-360	360
20	23	0.0028	0.0216	0.0455	625	0	0	0	0	1	-360	360
20	23	0.0028	0.0216	0.0455	625	0	0	0	0	1	-360	360
21	22	0.0087	0.0678	0.1424	625	0	0	0	0	1	-360	360
];

%% branch data
%	fbus	tbus  r	       x	   b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.ne_branch = [
1	2	0.0026	0.0139	0.4611	200	0	0	0	0	1	-360	360 3
1	3	0.0546	0.2112	0.0572	220	0	0	0	0	1	-360	360 55
1	5	0.0218	0.0845	0.0229	220	0	0	0	0	1	-360	360 22
1	8	0.0348	0.1344	0	220	0	0	0	0	1	-360	360  35
2	4	0.0328	0.1267	0.0343	220	0	0	0	0	1	-360	360 33 
2	6	0.0497	0.192	0.052	220	0	0	0	0	1	-360	360 50
2	8	0.0328	0.1267	0	220	0	0	0	0	1	-360	360 33
3	9	0.0308	0.119	0.0322	220	0	0	0	0	1	-360	360  31
3	24	0.0023	0.0839	0	600	0	0	0	0	1	-360	360 50
4	9	0.0268	0.1037	0.0281	220	0	0	0	0	1	-360	360 27
5	10	0.0228	0.0883	0.0239	220	0	0	0	0	1	-360	360 23
6	7	0.0497	0.192	0	220	0	0	0	0	1	-360	360 50
6  10   0.0139  0.0605  2.4590  200  0 0   1   0   1 -360	360 16
7   8   0.0159  0.0614  0.0166  220  0 0   1   0   1 -360	360 16
8	9	0.0427	0.1651	0.0447	220	0	0	0	0	1	-360	360 43
8	10	0.0427	0.1651	0.0447	220	0	0	0	0	1	-360	360 43
9	11	0.0023	0.0839	0	600	0	0	0	0	1	-360	360 50
9	12	0.0023	0.0839	0	600	0	0	0	0	1	-360	360 50 
10	11	0.0023	0.0839	0	600	0	0	0	0	1	-360	360 50 
10	12	0.0023	0.0839	0	600	0	0	0	0	1	-360	360 50
11	13	0.0061	0.0476	0.0999	625	0	0	0	0	1	-360	360  66
11	14	0.0054	0.0418	0.0879	625	0	0	0	0	1	-360	360  58
12	13	0.0061	0.0476	0.0999	625	0	0	0	0	1	-360	360  66
12	23	0.0124	0.0966	0.203	625	0	0	0	0	1	-360	360 134
13	14	0.0057	0.0447	0	625	0	0	0	0	1	-360	360 62
13	23	0.0111	0.0865	0.1818	625	0	0	0	0	1	-360	360 120
14	16	0.005	0.0389	0.0818	625	0	0	0	0	1	-360	360 54
14	23	0.008	0.062	0	625	0	0	0	0	1	-360	360 86
15	16	0.0022	0.0173	0.0364	625	0	0	0	0	1	-360	360 24
15	21	0.0063	0.049	0.1030	625	0	0	0	0	1	-360	360 68
15	24	0.0067	0.0519	0.1091	625	0	0	0	0	1	-360	360 72
16	17	0.0033	0.0259	0.0545	625	0	0	0	0	1	-360	360 36
16	19	0.003	0.0231	0.0485	625	0	0	0	0	1	-360	360 32
16	23	0.0105	0.0822	0	625	0	0	0	0	1	-360	360 114
17	18	0.0018	0.0144	0.0303	625	0	0	0	0	1	-360	360  20
17	22	0.0135	0.1053	0.2212	625	0	0	0	0	1	-360	360 146
18	21	0.0033	0.0259	0.0545	625	0	0	0	0	1	-360	360 36
19	20	0.0051	0.0396	0.0833	625	0	0	0	0	1	-360	360 55 
19	23	0.0078	0.0606	0	625	0	0	0	0	1	-360	360 84
20	23	0.0028	0.0216	0.0455	625	0	0	0	0	1	-360	360  30
21	22	0.0087	0.0678	0.1424	625	0	0	0	0	1	-360	360 94
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	cn	c(n-1)	...	c0
mpc.gencost =  [2   0	0	3	0 1*10^-2*0 0; %1  termica
                2   0   0   3   0 1*10^-2*0 0;  %2 termica
                2	0	0	3	0 1*10^-2*0 0; %3
                2	0	0	3	0 1000 0; %4
                2	0	0	3	0 0 0; %4
                2	0	0	3	0 1000 0; %5
                2	0	0	3	0 0 0; %5
                2	0	0	3	0 1*10^-2*0 0; %6
                2	0	0	3   0 1*10^-2*0 0; %7 termica
                2   0   0   3   0 1000 0; %8
                2   0   0   3   0 0 0; %8
                2   0   0   3   0 1*10^-2*0 0; %9
                2   0   0   3   0 1000 0; %10
                2   0   0   3   0 0 0; %10
                2   0   0   3   0 1000 0; %11
                2   0   0   3   0 0 0; %11
                2   0   0   3   0 1000 0; %12
                2   0   0   3   0 0 0; %12
                2   0   0   3   0 1*10^-2*0 0; %13
                2   0   0   3   0 1*10^-2*0 0; %14
                2   0   0   3   0 1*10^-2*0 0; %15
                2   0   0   3   0 1*10^-2*0 0; %16
                2   0   0   3   0 1000 0; %17
                2   0   0   3   0 0 0; %17
                2   0   0   3   0 1*10^-2*0 0; %18
                2   0   0   3   0 1000 0; %19
                2   0   0   3   0 0 0; %19
                2   0   0   3   0 1000 0; %20
                2   0   0   3   0 0 0; %20
                2   0   0   3   0 1*10^-2*0 0; %21
                2   0   0   3   0 1*10^-2*0 0; %22
                2   0   0   3   0 1*10^-2*0 0; %23
                2   0   0   3   0 1000 0; %24
                2   0   0   3   0 0 0; %24
                2   0	0	3	0 0 0; %1  termica
                2   0   0   3   0 0 0;  %2 termica
                2	0	0	3	0 0 0; %3
                2	0	0	3	0 0 0; %4
                2	0	0	3	0 0 0; %4
                2	0	0	3	0 0 0; %5
                2	0	0	3	0 0 0; %5
                2	0	0	3	0 0 0; %6
                2	0	0	3   0 0 0; %7 termica
                2   0   0   3   0 0 0; %8
                2   0   0   3   0 0 0; %8
                2   0   0   3   0 0 0; %9
                2   0   0   3   0 0 0; %10
                2   0   0   3   0 0 0; %10
                2   0   0   3   0 0 0; %11
                2   0   0   3   0 0 0; %11
                2   0   0   3   0 0 0; %12
                2   0   0   3   0 0 0; %12
                2   0   0   3   0 0 0; %13
                2   0   0   3   0 0 0; %14
                2   0   0   3   0 0 0; %15
                2   0   0   3   0 0 0; %16
                2   0   0   3   0 0 0; %17
                2   0   0   3   0 0 0; %17
                2   0   0   3   0 0 0; %18
                2   0   0   3   0 0 0; %19
                2   0   0   3   0 0 0; %19
                2   0   0   3   0 0 0; %20
                2   0   0   3   0 0 0; %20
                2   0   0   3   0 0 0; %21
                2   0   0   3   0 0 0; %22
                2   0   0   3   0 0 0; %23
                2   0   0   3   0 0 0; %24
                2   0   0   3   0 0 0; %24
                 ]; 
             
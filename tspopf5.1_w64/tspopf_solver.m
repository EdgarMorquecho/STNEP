function [results, success, raw] = tspopf_solver(om, mpopt)
%TSPOPF_SOLVER  Solves AC optimal power flow using TSPOPF.
%
%   [RESULTS, SUCCESS, RAW] = TSPOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an OPF model object and a MATPOWER options struct.
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   RESULTS is a MATPOWER case struct (mpc) with the usual baseMVA, bus
%   branch, gen, gencost fields, along with the following additional
%   fields:
%       .order      see 'help ext2int' for details of this field
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .nln
%               .l  lower bounds on nonlinear constraints
%               .u  upper bounds on nonlinear constraints
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
%       .info   solver specific termination code
%       .g      constraint values
%       .dg     constraint 1st derivatives
%
%   See also OPF.

%   MATPOWER
%   $Id: tspopf_solver.m 2595 2015-02-10 18:25:06Z ray $
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 2000-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% form vector of indices into table columns; needed by minopf.mex$(PLATFORM)
col = cell(1, 59);
[col{1:21}]  = idx_bus;     %% PQ ... MU_VMIN
[col{22:35}] = idx_gen;     %% GEN_BUS ... MU_QMIN
[col{36:42}] = idx_cost;    %% PW_LINEAR ... COST
[col{43:59}] = idx_brch;    %% F_BUS ... MU_ST
col = cell2mat(col);

%% options
alg = upper(mpopt.opf.ac.solver);

%% initialize some default options
if mpopt.pdipm.feastol == 0
  mpopt.pdipm.feastol = mpopt.opf.violation;  	%% = MPOPT.opf.violation by default
end
if mpopt.tralm.feastol == 0
  mpopt.tralm.feastol = mpopt.opf.violation;  	%% = MPOPT.opf.violation by default
end
if strcmp(alg, 'PDIPM')
  feastol = mpopt.pdipm.feastol;
else                        %% alg == 'TRALM'
  feastol = mpopt.tralm.feastol;
end

%% check for active power or current line flow limit
if upper(mpopt.opf.flow_lim(1)) ~= 'S'
  error('tspopf_solver: MPOPT.opf.flow_lim == ''%s'' not supported by MPOPT.opf.ac.solver == ''%s''', ...
        mpopt.opf.flow_lim, mpopt.opf.ac.solver);
end

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
cp = get_cost_params(om);
[N, H, Cw] = deal(cp.N, cp.H, cp.Cw);
fparm = [cp.dd cp.rh cp.kk cp.mm];
[vv, ll] = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs
nxyz = getN(om, 'var');     %% total number of control vars of all types
nxy = vv.N.Va + vv.N.Vm + vv.N.Pg + vv.N.Qg + ny;   %% number of x, y vars
nz = nxyz - nxy;            %% number of user vars
i1z = nxy + 1;              %% starting index for user vars

%% get user var set names
vs = get(om, 'var', 'order');   %% all var sets
k = 5 + (ny > 0);   %% everything following 'Qg', unless we also have 'y'
uvs = cell(1, length(vs) - k + 1);
if ~isempty(uvs)
  [uvs{:}] = deal(vs(k:end).name);
end

%% get bounds on user vars
[x0, LB, UB] = getv(om);
[z0, zl, zu] = deal(x0(i1z:nxyz), LB(i1z:nxyz), UB(i1z:nxyz));

%% WORKAROUND: TSPOPF MEX files assume N does NOT have columns for y
if ~isempty(N) && ny > 0
  N(:, vv.i1.y:vv.iN.y) = [];   %% strip out y columns
end

%% linear constraints
%% WORKAROUND: Add bounds on z vars to A, l, u
%% WORKAROUND: Move equality var bounds to general linear constraints
%%             where they are handled properly as equality constraints
om2 = [];
if ~isempty(zl) || ~isempty(zu)
  if isempty(zl)
    zl = -Inf(nz,1);
  end
  if isempty(zu)
    zu = Inf(nz,1);
  end
  if isempty(z0)
    z0 = zeros(size(zl));
  end
end
zb  = find(zl ~= -Inf | zu ~= Inf);     %% indices of z w/bounds
veq = find(abs(bus(:, VMAX) - bus(:, VMIN)) < feastol); %% equality bounds on Vm
peq = find(abs(gen(:, PMAX) - gen(:, PMIN)) < feastol); %% equality bounds on Pg
qeq = find(abs(gen(:, QMAX) - gen(:, QMIN)) < feastol); %% equality bounds on Qg
if ~isempty(zb) || ~isempty(veq) || ~isempty(peq) || ~isempty(qeq)
  om2 = om;
  if ~isempty(zb)
    nzb = length(zb);
    Azb = sparse(1:nzb, zb, 1, nzb, nz);
    om2 = add_constraints(om2, 'zlims', Azb, zl(zb), zu(zb), uvs);
  end
  relax = 1;
  if ~isempty(veq)
    nveq = length(veq);
    Vmin = bus(veq, VMIN);
    Vmax = bus(veq, VMAX);
    Aveq = sparse(1:nveq, veq, 1, nveq, vv.N.Vm);
    bveq = (Vmin + Vmax) / 2;
    om2 = add_constraints(om2, 'veq', Aveq, bveq, bveq, {'Vm'});
    bus(veq, VMIN) = Vmin - relax;
    bus(veq, VMAX) = Vmax + relax;
  end
  if ~isempty(peq)
    npeq = length(peq);
    Pmin = gen(peq, PMIN);
    Pmax = gen(peq, PMAX);
    Apeq = sparse(1:npeq, peq, 1, npeq, vv.N.Pg);
    bpeq = (Pmin + Pmax) / 2 / baseMVA;
    om2 = add_constraints(om2, 'peq', Apeq, bpeq, bpeq, {'Pg'});
    gen(peq, PMIN) = Pmin - relax;
    gen(peq, PMAX) = Pmax + relax;
  end
  if ~isempty(qeq)
    nqeq = length(qeq);
    Qmin = gen(qeq, QMIN);
    Qmax = gen(qeq, QMAX);
    Aqeq = sparse(1:nqeq, qeq, 1, nqeq, vv.N.Qg);
    bqeq = (Qmin + Qmax) / 2 / baseMVA;
    om2 = add_constraints(om2, 'qeq', Aqeq, bqeq, bqeq, {'Qg'});
    gen(qeq, QMIN) = Qmin - relax;
    gen(qeq, QMAX) = Qmax + relax;
  end
end

if isempty(om2)
  [A, l, u] = linear_constraints(om);
else
  [vv, ll] = get_idx(om2);
  [A, l, u] = linear_constraints(om2);
end

%% WORKAROUND: MEX file appears to have V before theta in the ordering of vars
if ~isempty(A)
    tmp = A(:, vv.i1.Va:vv.iN.Va);
    A(:, vv.i1.Va:vv.iN.Va) = A(:, vv.i1.Vm:vv.iN.Vm);
    A(:, vv.i1.Vm:vv.iN.Vm) = tmp;
end
if ~isempty(N)
    tmp = N(:, vv.i1.Va:vv.iN.Va);
    N(:, vv.i1.Va:vv.iN.Va) = N(:, vv.i1.Vm:vv.iN.Vm);
    N(:, vv.i1.Vm:vv.iN.Vm) = tmp;
end

%% replace Inf with numerical proxies
l(l == -Inf) = -1e10;
u(u ==  Inf) =  1e10;
bus(bus(:) == -Inf) = -1e10;
bus(bus(:) ==  Inf) =  1e10;
gen(gen(:) == -Inf) = -1e10;
gen(gen(:) ==  Inf) =  1e10;
branch(branch(:) == -Inf) = -1e10;
branch(branch(:) ==  Inf) =  1e10;

%% so, can we do anything good about lambda initialization?
if all(bus(:, LAM_P) == 0)
  bus(:, LAM_P) = (10)*ones(nb, 1);
end

%% call the solver
if strcmp(alg, 'PDIPM')
  if mpopt.pdipm.step_control
	[buso, geno, brancho, f, info, g, jac, xr, pimul] = scpdipmopf(baseMVA, ...
			bus, gen, gencost, branch, col, A, l, u, mpoption(mpopt, []), ...
			N, fparm, H, Cw, z0, zl, zu);
  else
	[buso, geno, brancho, f, info, g, jac, xr, pimul] = pdipmopf(baseMVA, ...
			bus, gen, gencost, branch, col, A, l, u, mpoption(mpopt, []), ...
			N, fparm, H, Cw, z0, zl, zu);
  end
elseif strcmp(alg, 'TRALM')
  [buso, geno, brancho, f, info, g, jac, xr, pimul] = tralmopf(baseMVA, ...
          bus, gen, gencost, branch, col, A, l, u, mpoption(mpopt, []), ...
          N, fparm, H, Cw, z0, zl, zu);
else
  error('tspopf_solver.m: ''%s'' is not a valid value for MPOPT.opf.ac.solver', alg);
end
success = (info == 1);

%% Filter multipliers for non-binding constraints
% buso(find((buso(:, VMAX) - buso(:, VM)) > feastol*100), MU_VMAX) = 0;
% buso(find((buso(:, VM)   - buso(:, VMIN)) > feastol*100), MU_VMIN) = 0;
% geno(find((geno(:, PMAX) - geno(:, PG))/baseMVA > feastol*100), MU_PMAX) = 0;
% geno(find((geno(:, QMAX) - geno(:, QG))/baseMVA > feastol*100), MU_QMAX) = 0;
% geno(find((geno(:, PG)   - geno(:, PMIN))/baseMVA > feastol*100), MU_PMIN) = 0;
% geno(find((geno(:, QG)   - geno(:, QMIN))/baseMVA > feastol*100), MU_QMIN) = 0;
% Sf = sqrt(brancho(:, PF).^2 + brancho(:, QF).^2);
% St = sqrt(brancho(:, PT).^2 + brancho(:, QT).^2);
% brancho(find((brancho(:, RATE_A) - Sf)/baseMVA > feastol*100), MU_SF) = 0;
% brancho(find((brancho(:, RATE_A) - St)/baseMVA > feastol*100), MU_ST) = 0;

%% WORKAROUND: MEX files do not copy these columns to output matrices
brancho(:, [ANGMIN, ANGMAX]) = branch(:, [ANGMIN, ANGMAX]);
geno(:, PC1:APF) = gen(:, PC1:APF);

%% package up results
nlnN = getN(om, 'nln');
linN = getN(om, 'lin');

%% WORKAROUND: Reorder pimul to account for adding bounds on user vars to A, l, u
%%             and moving equality var bounds to general linear constraints
%%             Copy multipliers back to output data matrices.
if ~isempty(om2)
  del = [];
  np = length(pimul);
  if ~isempty(zb)
    vb = (np-nz+1:np);                      %% idx of relevant variable bounds
    k = nlnN + (ll.i1.zlims:ll.iN.zlims);   %% idx of corresp linear constraints
    pimul(vb(zb)) = pimul(k);
    del = [del k];
  end
  if ~isempty(veq)
    buso(veq, VMIN) = Vmin;                 %% restore original limits
    buso(veq, VMAX) = Vmax;
    vb = np - nxyz + (vv.i1.Vm:vv.iN.Vm);   %% idx of relevant variable bounds
    k = nlnN + (ll.i1.veq:ll.iN.veq);       %% idx of corresp linear constraints
    pimul(vb(veq)) = pimul(k);
    kl = find(pimul(k) >= 0);
    ku = find(pimul(k) <  0);
    buso(veq(kl), MU_VMIN) =  pimul(k(kl));
    buso(veq(ku), MU_VMAX) = -pimul(k(ku));
    del = [del k];
  end
  if ~isempty(peq)
    geno(peq, PMIN) = Pmin;                 %% restore original limits
    geno(peq, PMAX) = Pmax;
    vb = np - nxyz + (vv.i1.Pg:vv.iN.Pg);   %% idx of relevant variable bounds
    k = nlnN + (ll.i1.peq:ll.iN.peq);       %% idx of corresp linear constraints
    pimul(vb(peq)) = pimul(k);
    kl = find(pimul(k) >= 0);
    ku = find(pimul(k) <  0);
    geno(peq(kl), MU_PMIN) =  pimul(k(kl)) / baseMVA;
    geno(peq(ku), MU_PMAX) = -pimul(k(ku)) / baseMVA;
    del = [del k];
  end
  if ~isempty(qeq)
    geno(qeq, QMIN) = Qmin;                 %% restore original limits
    geno(qeq, QMAX) = Qmax;
    vb = np - nxyz + (vv.i1.Qg:vv.iN.Qg);   %% idx of relevant variable bounds
    k = nlnN + (ll.i1.qeq:ll.iN.qeq);       %% idx of corresp linear constraints
    pimul(vb(qeq)) = pimul(k);
    kl = find(pimul(k) >= 0);
    ku = find(pimul(k) <  0);
    geno(qeq(kl), MU_QMIN) =  pimul(k(kl)) / baseMVA;
    geno(qeq(ku), MU_QMAX) = -pimul(k(ku)) / baseMVA;
    del = [del k];
  end
  pimul(del) = [];                          %% delete all these "extras"
end

%% WORKAROUND: Remove dummy row for piecewise linear costs if there are none
any_pwl = any(gencost(:, MODEL) == PW_LINEAR);
if ~any_pwl
  pimul(nlnN+linN+1) = [];
end

%% ... continue to package up results
lam_nl  = pimul(1:nlnN);
lam_ln  = pimul((1:linN) + nlnN);
lam_var = pimul((1:nxyz) + nlnN+linN+any_pwl);
mu = struct( ...
  'var', struct('l', zeros(nxyz,1), 'u', zeros(nxyz,1)), ...
  'nln', struct('l', zeros(nlnN,1), 'u', zeros(nlnN,1)), ...
  'lin', struct('l', zeros(linN,1), 'u', zeros(linN,1)) );
idx1 = find(lam_var > 0);
idx2 = find(lam_var < 0);
idx3 = find(lam_nl  > 0);
idx4 = find(lam_nl  < 0);
idx5 = find(lam_ln  > 0);
idx6 = find(lam_ln  < 0);
mu.var.l(idx1) =  lam_var(idx1);
mu.var.u(idx2) = -lam_var(idx2);
mu.nln.l(idx3) =  lam_nl(idx3);
mu.nln.u(idx4) = -lam_nl(idx4);
mu.lin.l(idx5) =  lam_ln(idx5);
mu.lin.u(idx6) = -lam_ln(idx6);

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(buso, brancho, geno, om, xr(1:nxyz), mu, f);

raw = struct('xr', xr, 'pimul', pimul, 'info', info, 'g', g, 'dg', jac);

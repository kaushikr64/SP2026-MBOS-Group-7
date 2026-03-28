clc; clearvars; close all; restoredefaultpath;

addpath("Library\");
addpath("Export_Fig\");
load("EarthMoon_SystemConstants.mat");

%% Choose the Specific Orbit Family
% Here you can load in different orbit families, from pre-computed
% libraries (see directory). This section will plot the strictly periodic
% orbit family in the rotating frame. Additionally, a second plot is
% generated showing the Brouke stability values to aid in determining which
% specific orbits admit nerby center subspaces.

tmp = load("EM_CR3BP_FamData\L1_NorthernHaloFamily_Data.mat");
fn = fieldnames(tmp);
OrbFmly = tmp.(fn{1});

% Trajectory Plots
f1 = figure();
plot_Moon(100,EMsys.Xmoon*EMsys.Ls);
hold on;
grid on;
plot3(EMsys.Lp(1,1:2)*EMsys.Ls,EMsys.Lp(2,1:2)*EMsys.Ls,EMsys.Lp(3,1:2)*...
    EMsys.Ls,'LineStyle','none','Marker','^','MarkerSize',5,'LineWidth',2,'color',[0,0,0]);
for i = 1:10:OrbFmly.nOrb
    plot3(OrbFmly.X_hst(1,:,i)*EMsys.Ls,OrbFmly.X_hst(2,:,i)*EMsys.Ls,...
        OrbFmly.X_hst(3,:,i)*EMsys.Ls,'LineWidth',1,'color',...
        [0.5,0.5,0.9]);
end
axis equal;
xlabel('X [km]','Interpreter','latex','FontSize',18);
ylabel('Y [km]','Interpreter','latex','FontSize',18);
zlabel('Z [km]','Interpreter','latex','FontSize',18);
title('Earth-Moon Periodic Family','Interpreter','latex','FontSize',18);
legend({'Moon','$$L_{1,2}$$'},'Interpreter','latex','FontSize',14,'location','best');
f1.Position = [10,60,1500,1200];

% Broucke Stability Plot
f2 = figure();
plot(OrbFmly.BrkeStab(1,:),OrbFmly.BrkeStab(2,:),'LineStyle','--','Marker',...
    'x','color',[0,0,0]);
text(OrbFmly.BrkeStab(1,:),OrbFmly.BrkeStab(2,:),string(1:OrbFmly.nOrb));
grid on;
hold on;
fplot(@(a) -2.*a-2,'-','LineWidth',1.5); % Tangent
fplot(@(a) 2.*a-2,'--','LineWidth',1.5); % Period doubling
fplot(@(a) a+1,'-.','LineWidth',1.5); % Period Tripling
fplot(2,':','LineWidth',1.5); % Period Quadrupling
fplot(@(a) (a.^2)./4+2,'LineWidth',1.5); % Secondary Hopf
% axis equal;
xlabel('$$\alpha$$','Interpreter','latex','FontSize',18);
ylabel('$$\beta$$','Interpreter','latex','FontSize',18);
legend('Periodic Orbit Family','Tangent','Period Doubling','Period Tripling',...
    'Period Quadrupling','Secondary Hopf','Interpreter','latex','FontSize',14);
title('Family Brouke Stability Diagram','Interpreter','latex','FontSize',18);
f2.Position = [10,60,1500,1200];


%% Identify Initial Orbit for QPO Generation ------------------------------
% From the family just loaded, pick the specific orbit from which you would
% like to generate and continue quasi-periodic orbits from (QPOs).
% Additionally, specificy the poincare map (plane) via a point and normal
% direction (in 6D space). 
Per_Orb.ID = 82;
section.X_sec = [EMsys.Lp(:,1)+[0;0;0];zeros(3,1)]; % L1 Point
section.n_sec = [0;0;1;zeros(3,1)]; % Pointing along +Z
% -------------------------------------------------------------------------

% Compute the section values (zero if point lies exactly on section)
Per_Orb.Section.vals = PoincareSec(OrbFmly.X_hst(:,:,Per_Orb.ID), section.X_sec, section.n_sec);

% Order orbit states based on how close to poincare section they are
[Per_Orb.Section.vals_sort, Per_Orb.Section.IDs_sort] = sort(abs(Per_Orb.Section.vals));
Per_Orb.Section.X_sort = OrbFmly.X_hst(:,Per_Orb.Section.IDs_sort,Per_Orb.ID);

%% Pick the point along the initial orbit to rephase to -------------------
% Here we reshoot for an initial state on the Poincare map that is the
% periodic orbit (fixed point).Be careful with symmetric orbtis
% to pick the point that has the velocity direction you want. If this is
% picked incorrectly, it will fail to converge because the return is at the
% half...ish period.
Per_Orb.RePhased.ID0 = 1;
% -------------------------------------------------------------------------

% Grab the state cooresponding to the point of the orbit chosen. 
Per_Orb.RePhased.X0 = Per_Orb.Section.X_sort(:,Per_Orb.RePhased.ID0);

% Shoot for the exact state (fixed point) on the poincare map
fun = @(X) CR3BP_PerOrb_SF(X, EMsys.mu, section);
OutFuncTolerance = 1e-12; % Max difference in state

opts_newt.tol = 1e-12;
opts_newt.max_step = 1e-2;
opts_newt.max_iter = 50;
Per_Orb.RePhased.X0 = simple_newton(fun,Per_Orb.RePhased.X0,opts_newt);
Per_Orb.JC = CR3BP_JC(Per_Orb.RePhased.X0,EMsys.mu);

%% Integrate Periodic Orbit for QPO Intial Conditions
prop_opts.RelTol = 1e-12;
prop_opts.AbsTol = 1e-12;
prop_opts.Direction = 1;

% Propagate to map (should return exactly) and return the Poincare map
% partial, the full STM, and the orbit period.
[~,Per_Orb.RePhased.DP,Per_Orb.RePhased.STM,Per_Orb.RePhased.T] = CR3BP_Prop2Poincare(Per_Orb.RePhased.X0,...
    EMsys.mu,section,[0,20],prop_opts);

% Solve for the eigenvalues and eigenvectors of the Poincare map partial
% wrt to the state.
[Per_Orb.RePhased.vec,Per_Orb.RePhased.val] = SortEigsPhi_2(Per_Orb.RePhased.DP);


%% Generate Initial guess of QPO Fourier Coefficients ---------------------
% Now that we have the eigenvalues and cooresponding eigenvectors, we can
% choose which specific term (one unit circle), to set to the QPO
% coefficients. 
QPO.vecval_ID = [3]; % Eigenvector/value term to use
QPO.vecval_blendw = 0.5; % Weighting between the two if multiple
u1 = [1;1;1;0;0;0]; % Area projection vector 1
u2 = [0;0;0;1;1;1]; % Area projection vector 2
Nc = 10; % Number of coefficients to use within the Fourier expansion
Np = 1+2*Nc+2; % Number of points to propagate
rho = 1e-2; % Distance to step of fixed point (FP) periodic orbit
% -------------------------------------------------------------------------

[QPO.init.fv_guess, QPO.init.A] = OrbFP2QPO_vals(Per_Orb.RePhased.X0,rho,...
    Per_Orb.RePhased.vec(:,QPO.vecval_ID),Per_Orb.RePhased.val(QPO.vecval_ID),...
    Nc,u1,u2,QPO.vecval_blendw); 

% Solve for the first QPO ------------------------------------------------

opts.mu = EMsys.mu;
opts.Np = Np;
opts.Nc = Nc;
opts.section = section;
opts.con_area.u1 = u1;
opts.con_area.u2 = u2;
opts.con_area.A = QPO.init.A;
opts.con_JC.JC = Per_Orb.JC;

opts_newt = struct;
opts_newt.max_step = 1e0;
opts_newt.max_iter = 100;
opts_newt.tol = [1e-8*ones(6*Np,1);1e-5;1e-2]; % [position, JC, rel_area]
fun = @(X) CR3BP_QPO_SF(X, opts);
QPO.init.fv_true = simple_newton(fun,QPO.init.fv_guess,opts_newt);
[F_sol,dFdfv_sol] = fun(QPO.init.fv_true);


%% Generate QPO family
N_QPOs = 40;
A_scale = 1.5;
opts_newt.max_iter = 10;

QPO.fmly.fv_hst(:,1) = QPO.init.fv_true;
A_pri = QPO.init.A;
dFdfv_pri = dFdfv_sol;

for i = 2:N_QPOs
    dfv = sol2nextfv(dFdfv_pri,A_pri,A_scale);
    fv_i = QPO.fmly.fv_hst(:,i-1) + dfv;
    
    A_pri = A_scale*A_pri;
    opts.con_area.A = A_pri;
    fun = @(X) CR3BP_QPO_SF(X, opts);
    QPO.fmly.fv_hst(:,i) = simple_newton(fun,fv_i,opts_newt);
    [~,dFdfv_pri] = fun(QPO.fmly.fv_hst(:,i));
    fprintf("Orbit i = "+i+" complete\n");

end

%% Plot specific QPO
i_plt = 10;
f1 = plotQPOOrbit_3D(QPO.fmly.fv_hst(1:end-1,i_plt),@(X0) CR3BP_Prop2Poincare(X0,EMsys.mu,section,[0,20]),...
    1+1*Np+2,Nc,EMsys);
plot_Moon(100,EMsys.Xmoon*EMsys.Ls);
plot3(EMsys.Lp(1,1:2)*EMsys.Ls, EMsys.Lp(2,1:2)*EMsys.Ls, EMsys.Lp(3,1:2)*EMsys.Ls, ...
    'LineStyle','none','LineWidth',2,'Marker','^','MarkerSize',5,'Color','k');
text(EMsys.Lp(1,1:2)*EMsys.Ls, EMsys.Lp(2,1:2)*EMsys.Ls, EMsys.Lp(3,1:2)*EMsys.Ls, ...
    ["$$L_1$$","$$L_2$$"],'FontSize',16,'Interpreter','latex')
xlabel("X",'interpreter','latex','FontSize',18);
ylabel("Y",'interpreter','latex','FontSize',18);
zlabel("Z",'interpreter','latex','FontSize',18);
title("Representative QPO from Family",'interpreter','latex','fontsize',20);
f1.Position = [10,60,1500,1000];

i_sec = 10;
f2 = figure();
subplot(2,1,1);
plotQPOFamilyPoinSection(QPO.fmly.fv_hst(1:end-1,1:i_sec),0,Nc,EMsys,section);
xlabel("X",'interpreter','latex','FontSize',18);
ylabel("Y",'interpreter','latex','FontSize',18);
zlabel("Z",'interpreter','latex','FontSize',18);
title("QPO Family - Poincare Section",'interpreter','latex','fontsize',20);
subplot(2,1,2);
plotQPO_JC(QPO.fmly.fv_hst(1:end-1,1:i_sec),Nc,EMsys);
xlabel('$$\theta$$','interpreter','latex','FontSize',18);
ylabel('Jacobi Constant','interpreter','latex','FontSize',18);
title("Jacobi Constant Evolution Along Discrete Fourier Series",'interpreter','latex','fontsize',20);
f2.Position = [10,60,800,1000];


%% Gradient Checks
% info = checkJacobian(fun,QPO.fmly.fv_hst(:,5));


function info = checkJacobian(fun, x0, opts)

if nargin < 3
    opts = struct();
end
if ~isfield(opts,'eps'), opts.eps = 1e-7; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

eps_fd = opts.eps;

x0 = x0(:);
n = length(x0);

% Evaluate function + analytic Jacobian
[F0, J_analytic] = fun(x0);
m = length(F0);

J_fd = zeros(m,n);

% Finite difference (central)
for j = 1:n
    dx = zeros(n,1);
    dx(j) = eps_fd;

    F_plus  = fun_only(fun, x0 + dx);
    F_minus = fun_only(fun, x0 - dx);

    J_fd(:,j) = (F_plus - F_minus) / (2*eps_fd);
end

% Error metrics
abs_err = abs(J_fd - J_analytic);
rel_err = abs_err ./ max(abs(J_fd), 1e-12);

[max_abs_err, idx_abs] = max(abs_err(:));
[row_abs, col_abs] = ind2sub(size(abs_err), idx_abs);

[max_rel_err, idx_rel] = max(rel_err(:));
[row_rel, col_rel] = ind2sub(size(rel_err), idx_rel);

% Output struct
info.max_abs_err = max_abs_err;
info.max_rel_err = max_rel_err;
info.worst_abs_idx = [row_abs, col_abs];
info.worst_rel_idx = [row_rel, col_rel];

info.J_fd = J_fd;
info.J_analytic = J_analytic;

% Verbose output
if opts.verbose
    fprintf('--- Jacobian Check ---\n');
    fprintf('Max abs error: %.3e at (%d,%d)\n', max_abs_err, row_abs, col_abs);
    fprintf('Max rel error: %.3e at (%d,%d)\n', max_rel_err, row_rel, col_rel);

    fprintf('\nWorst ABS entry:\n');
    fprintf('Analytic: %.6e\n', J_analytic(row_abs,col_abs));
    fprintf('FD      : %.6e\n', J_fd(row_abs,col_abs));

    fprintf('\nWorst REL entry:\n');
    fprintf('Analytic: %.6e\n', J_analytic(row_rel,col_rel));
    fprintf('FD      : %.6e\n', J_fd(row_rel,col_rel));
end

end

function F = fun_only(fun, x)
    F = fun(x);
    if iscell(F)
        F = F{1};
    end
end


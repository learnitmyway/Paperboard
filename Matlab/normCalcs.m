% Hand calcs for simplified umat - 1D normal stress

clear;
clc;

% Symbolic variables
syms S eps eps_p kap real;

% Material parameters
E_const = 1000;
G_const = 500;
nu = 0.0;
Sy0 = 1;
H = 100;

% Elasticity matrix
C_mat(1,1) = E_const;
C_mat(2,2) = G_const;
C_mat(3,3) = G_const;

% Point A
eps_pA = 0.0;
kapA = eps_pA;
Sy = Sy0 + H * kapA;
f = norm(S) - Sy == 0;
SA = solve(f, S > 0.0);
eqn = C_mat(1,1) * (eps - eps_pA) == SA;
epsA = solve(eqn, eps > 0.0);

% Point B
epsB = 5e-3;
kap = eps_p;
Sy = Sy0 + H * kap;
S = C_mat(1,1) * (epsB - eps_p);
f = norm(S) - Sy == 0;
eps_pB = solve(f, eps_p);
eps_pB = eps_pB(1);
kapB = eps_pB;
SB = C_mat(1,1) * (epsB - eps_pB);

% Redfine symbolic variables
syms S eps eps_p kap real;

% Point C
eps_pC = eps_pB;
kapC = kapB;
Sy = Sy0 + H * kapC;
f = norm(S) - Sy == 0;
SC = solve(f, S<SB);
eqn = C_mat(1,1) * (eps - eps_pC) == SC;
epsC = solve(eqn, eps);

% Point D
epsD = -5e-3;
kap = 2 * eps_pB - eps_p;
Sy = Sy0 + H * kap;
S = C_mat(1,1) * (epsD - eps_p);
f = norm(S) - Sy == 0;
eps_pD = solve(f, eps_p < eps_pB);
eps_pD = eps_pD(2);
kapD = 2 * eps_pB - eps_pD;
SD = C_mat(1,1) * (epsD - eps_pD);

% Final results
eps33_an = [epsA,epsB,epsC,epsD]
S33_an = [SA,SB,SC,SD]

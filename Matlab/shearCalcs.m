% Hand calcs for simplified umat - 1D shear

clear;
clc;

% Symbolic variables
syms S gam gam_p kap real;

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
gam_pA = 0.0;
kapA = gam_pA;
Sy = Sy0 + H * kapA;
f = norm(S) - Sy == 0;
SA = solve(f, S > 0.0);
eqn = C_mat(3,3) * (gam - gam_pA) == SA;
gamA = solve(eqn, gam > 0.0);

% Point B
gamB = 2 * 5e-3;
kap = gam_p;
Sy = Sy0 + H * kap;
S = C_mat(3,3) * (gamB - gam_p);
f = norm(S) - Sy == 0;
gam_pB = solve(f, gam_p);
gam_pB = gam_pB(1);
kapB = gam_pB;
SB = C_mat(3,3) * (gamB - gam_pB);

% Redfine symbolic variables
syms S gam gam_p kap real;

% Point C
gam_pC = gam_pB;
kapC = kapB;
Sy = Sy0 + H * kapC;
f = norm(S) - Sy == 0;
SC = solve(f, S < SB);
eqn = C_mat(3,3) * (gam - gam_pC) == SC;
gamC = solve(eqn, gam);

% Point D
gamD = 2 * -5e-3;
kap = 2 * gam_pB - gam_p;
Sy = Sy0 + H * kap;
S = C_mat(3,3) * (gamD - gam_p);
f = norm(S) - Sy == 0;
gam_pD = solve(f, gam_p < gam_pB);
gam_pD = gam_pD(2);
kapD = 2 * gam_pB - gam_pD;
SD = C_mat(3,3) * (gamD - gam_pD);

gam_an = [gamA,gamB,gamC,gamD]
S_an = [SA,SB,SC,SD]

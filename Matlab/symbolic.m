% Gradients for simplified umat

clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables
syms H kap Sy Sy0 gam
S = sym('S', [3 1]);

% Symbolic variables
S1 = S(1);
S2 = S(2);
S3 = S(3);

% Yield function
f_pre = norm(S) - Sy;

% Yield function in scalar form
f = (S1^2 + S2^2 + S3^2)^(1/2) - Sy;

% flow vector = gradient of yield function with respect to S
f_S = [diff(f, S(1)); diff(f, S(2)); diff(f, S(3))];

% Magnitude of yield function
mag_f_S_pre = norm(f_S);

% Magnitude of yield function in scalar form
mag_f_S = (S1^2/(S1^2 + S2^2 + S3^2) + (S2)^2/(S1^2 + S2^2 + S3^2) + S3^2/(S1^2 + S2^2 + S3^2))^(1/2);
 
r = f_S / mag_f_S;

% gradient of flow vector with respect to S
r_S = ...
    [diff(r(1), S(1)), diff(r(1), S(2)), diff(r(1), S(3));...
    diff(r(2), S(1)), diff(r(2), S(2)), diff(r(2), S(3));...
    diff(r(3), S(1)), diff(r(3), S(2)), diff(r(3), S(3))];

% gradient of f with respect to S
f_S2 = ...
    [diff(f_S(1), S(1)), diff(f_S(1), S(2)), diff(f_S(1), S(3));...
    diff(f_S(2), S(1)), diff(f_S(2), S(2)), diff(f_S(2), S(3));...
    diff(f_S(3), S(1)), diff(f_S(3), S(2)), diff(f_S(3), S(3))];

% Hardening function - linear isotropic
Sy = Sy0 + H * kap;

% Derivative of Sy with respect to accumulated plastic strain
Sy_kap = diff(Sy, kap);

% Derivative of yield function with respect to Sy
f_Sy = diff(f, 'Sy');

% Gradient of normalised flow vector with respect to kap
r_Sy = ...
    [diff(r(1), 'Sy'); ...
    diff(r(2), 'Sy'); ...
    diff(r(3), 'Sy')];

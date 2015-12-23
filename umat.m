% umat for simplified example

function [stress,statev,ddsdde] = umat(stress,dstran,stran,statev,ddsdde,time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intialise elastic and plastic trial strains

eps_e = statev(1:6) + dstran;
eps_p = statev(7:12);

% Initialise accumulated plastic trial strain
kap = statev(13);

% Tolerance for yield function
tol = 10^-6;

% Max number of iterations for Newton-Rapshon
Newton = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material parameters
E_const = 1000;
G_const = 500;
nu = 0.0;
Sy0 = 1;
H = 100;

% Algorithmic modulus
ddsdde = zeros(6,6);
ddsdde(3,3) = E_const;
ddsdde(5,5) = G_const;
ddsdde(6,6) = G_const;

% Elasticity matrix
C_mat = zeros(3,3);
C_mat(1,1) = ddsdde(3,3);
C_mat(2,2) = ddsdde(5,5);
C_mat(3,3) = ddsdde(6,6);

% Stress vector (6x1)
% stress = stress + ddsdde * dstran;
stress = ddsdde * eps_e;

% Trial stress vector (3x1)
sigma_trial(1,1) = stress(3);
sigma_trial(2,1) = stress(5);
sigma_trial(3,1) = stress(6);

% Hardening function - linear isotropic hardening
Sy = Sy0 + H * kap;

% Yield function
f = norm(sigma_trial) - Sy;

% Plastic corrector step (Return-Mapping Algorithm)  

if f > tol
    
    % Initialisation (k = 0)
    Dlam = 0;
    sigma = sigma_trial;
    
    % Newton-Raphson
    
    for k = 1: Newton
           
        % Check for convergence

        if f < tol && norm([a;b]) < tol
            break
        end % if f < tol

        % Symbolic variables 
        S1 = sigma(1);
        S2 = sigma(2);
        S3 = sigma(3);

        % Substitute variables
        v1 = (S1^2 + S2^2 + S3^2);

        % Gradient of yield function with respect to sigma
        f_S(1,1) = S1/v1^(1/2);
        f_S(2,1) = S2/v1^(1/2);
        f_S(3,1) = S3/v1^(1/2);

        % Normalised flow vector 
        r = f_S / norm(f_S);

        % Gradient of flow vector with respect to sigma  
        r_S(1,1) = 1/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(1/2)) - S1^2/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2)) + (S1*((2*S1^3)/v1^2 - (2*S1)/v1 + (2*S1*S2^2)/v1^2 + (2*S1*S3^2)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)); 
        r_S(1,2) = (S1*((2*S2^3)/v1^2 - (2*S2)/v1 + (2*S1^2*S2)/v1^2 + (2*S2*S3^2)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)) - (S1*S2)/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2)); 
        r_S(1,3) = (S1*((2*S3^3)/v1^2 - (2*S3)/v1 + (2*S1^2*S3)/v1^2 + (2*S2^2*S3)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)) - (S1*S3)/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2));
        r_S(2,1) = (S2*((2*S1^3)/v1^2 - (2*S1)/v1 + (2*S1*S2^2)/v1^2 + (2*S1*S3^2)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)) - (S1*S2)/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2)); 
        r_S(2,2) = 1/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(1/2)) - S2^2/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2)) + (S2*((2*S2^3)/v1^2 - (2*S2)/v1 + (2*S1^2*S2)/v1^2 + (2*S2*S3^2)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)); 
        r_S(2,3) = (S2*((2*S3^3)/v1^2 - (2*S3)/v1 + (2*S1^2*S3)/v1^2 + (2*S2^2*S3)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)) - (S2*S3)/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2));
        r_S(3,1) = (S3*((2*S1^3)/v1^2 - (2*S1)/v1 + (2*S1*S2^2)/v1^2 + (2*S1*S3^2)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)) - (S1*S3)/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2)); 
        r_S(3,2) = (S3*((2*S2^3)/v1^2 - (2*S2)/v1 + (2*S1^2*S2)/v1^2 + (2*S2*S3^2)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2)) - (S2*S3)/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2)); 
        r_S(3,3) = 1/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(1/2)) - S3^2/((S1^2/v1 + S2^2/v1 + S3^2/v1)^(1/2)*v1^(3/2)) + (S3*((2*S3^3)/v1^2 - (2*S3)/v1 + (2*S1^2*S3)/v1^2 + (2*S2^2*S3)/v1^2))/(2*(S1^2/v1 + S2^2/v1 + S3^2/v1)^(3/2)*v1^(1/2));

        % gradient of f with respect to sigma
        f_S2(1,1) = 1/(v1)^(1/2) - S1^2/(v1)^(3/2);
        f_S2(1,2) = -(S1*S2)/(v1)^(3/2);
        f_S2(1,3) = -(S1*S3)/(v1)^(3/2);
        f_S2(2,1) = -(S1*S2)/(v1)^(3/2);
        f_S2(2,2) = 1/(v1)^(1/2) - S2^2/(v1)^(3/2);
        f_S2(2,3) = -(S2*S3)/(v1)^(3/2);
        f_S2(3,1) = -(S1*S3)/(v1)^(3/2);
        f_S2(3,2) = -(S2*S3)/(v1)^(3/2);
        f_S2(3,3) = 1/(v1)^(1/2) - S3^2/(v1)^(3/2);
        
        % Derivatives of hardening function
        Sy_kap = H;

        % Derivative of yield function with respect to Sy
        f_Sy = -1;

        % Derivative of yield function with respect to kap
        f_kap = f_Sy * Sy_kap;

        % Gradient of normalised flow vector with respect to Sy       
        r_Sy(1,1) = 0;
        r_Sy(2,1) = 0;
        r_Sy(3,1) = 0;

        % Gradient of normalised flow vector with respect to kap 
        r_kap = r_Sy * Sy_kap;
                
        % inverse of C
        invC = C_mat^-1;

        % increment in plasticity parameter
        
        a = invC * (sigma - sigma_trial) + Dlam * r;

        b = -kap + statev(13) + Dlam;
        
        invA = ...
            [invC + Dlam*r_S, Dlam*r_kap; ...
            zeros(1,3), -1];
        A = invA^-1;

        dlam = (f - [f_S', f_kap] * A * [a; b]) / ...
            ([f_S', f_kap] * A * [r; 1]);

        % Obtain increments in stress and internal variables

        DSkap = - A * [a; b] - dlam * A * [r; 1];
        DS = DSkap(1:3);
        Dkap = DSkap(4);

        % Update plastic strain and internal variables

        eps_p(3,1) = eps_p(3,1) - invC(1,:) * DS;
        eps_p(5,1) = eps_p(5,1) - invC(2,:) * DS;
        eps_p(6,1) = eps_p(6,1) - invC(3,:) * DS;

        eps_e(3,1) = eps_e(3,1) + invC(1,:) * DS;
        eps_e(5,1) = eps_e(5,1) + invC(2,:) * DS;
        eps_e(6,1) = eps_e(6,1) + invC(3,:) * DS;

        kap = kap + Dkap;

        Dlam = Dlam + dlam;

        sigma = sigma + DS;

        Sy = Sy0 + H * kap;

        f = norm(sigma) - Sy;

    end % for k=1 : Newton
    
    % If it does not converge
    
    if k > Newton - 1 && f > tol
        fprintf('Warning: did not converge\n');
    end
        
        % Store stress and strains in state variable array if plastic
        stress(3) = sigma(1);
        stress(5) = sigma(2);
        stress(6) = sigma(3);
        statev(1:6) = eps_e;
        statev(7:12) = eps_p;
        statev(13) = kap;
        
        % Update Algorithmic modulus if plastic
      
        C_mod = A(1:3,1:3);
        C_mod2 = invC + Dlam * f_S2;
        
        C_alg = C_mod - ((C_mod * r) * (f_S' * C_mod)) / ((f_S' * C_mod * r) - f_kap);
        C_alg2 = C_mod2 - ((C_mod2 * f_S) * (f_S' * C_mod2)) / (f_S' * C_mod2 * f_S + f_kap * norm(f_S));
        
        ddsdde(3,3) = C_alg(1,1);
        ddsdde(3,5) = C_alg(1,2);
        ddsdde(3,6) = C_alg(1,3);
        ddsdde(5,3) = C_alg(2,1);
        ddsdde(5,5) = C_alg(2,2);
        ddsdde(5,6) = C_alg(2,3);
        ddsdde(6,3) = C_alg(3,1);
        ddsdde(6,5) = C_alg(3,2);
        ddsdde(6,6) = C_alg(3,3);
        
    % Store stress and strains in state variable array if elastic   
else 
    sigma = sigma_trial;
    stress(3) = sigma(1);
    stress(5) = sigma(2);
    stress(6) = sigma(3);
    statev(1:6) = eps_e;
    statev(7:12) = eps_p;
    statev(13) = kap;
    
end % if f > tol (plasticity)

end % umat

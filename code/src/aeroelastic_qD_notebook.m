%% DIVERGENCE DYNAMIC PRESSURE q_D - MATLAB NOTEBOOK
% This script reproduces the latex derivation in a MATLAB notebook style.
% It is organized as a live-script compatible file with sections (%%).
% You can run it either as a standard .m script or open it in MATLAB and
% save it as a .mlx Live Script if desired.
%
% Main outputs:
%   1) symbolic expressions of structural and aerodynamic matrices
%   2) symbolic divergence dynamic pressure q_D
%   3) symbolic control effectiveness E_C(q)
%   4) optional numerical evaluation and plot if parameter values are given

clear; clc;

%% 1) Symbolic variables
% Geometric / structural / aerodynamic parameters
syms b ba c e xB zA gamma real positive
syms E J G A real positive
syms CL_alpha CL_beta Cm_beta real
syms beta q qD real
syms phi_w phi_theta real

% Assumptions useful for simplification
assumeAlso(gamma > 0);
assumeAlso(gamma < sym(pi)/2);

%% 2) Ritz shape functions
% Root clamp:
%   w(0)=0, w'(0)=0, theta(0)=0
% Chosen Ritz ansatz:
%   w(y)     = Nw(y) * phi_w
%   theta(y) = Nt(y) * phi_theta

syms y real
Nw = (y/b)^2;
Nt = y/b;

w_y     = Nw * phi_w;
theta_y = Nt * phi_theta;

fprintf('Shape functions:\n');
disp('Nw(y) ='); disp(Nw);
disp('Nt(y) ='); disp(Nt);

%% 3) Internal virtual work and structural stiffnesses
% delta Wi = int delta w,yy * EJ * w,yy dy + int delta theta,y * GJ * theta,y dy

Nw_yy = diff(Nw, y, 2);
Nt_y  = diff(Nt, y, 1);

Kw     = simplify(int(Nw_yy * E*J * Nw_yy, y, 0, b));
Ktheta = simplify(int(Nt_y  * G*J * Nt_y , y, 0, b));

K_internal = [Kw, 0;
              0, Ktheta];

fprintf('\nInternal stiffness terms:\n');
disp('Kw ='); disp(Kw);
disp('Ktheta ='); disp(Ktheta);
disp('K_internal ='); disp(K_internal);

%% 4) Wing aerodynamic virtual work -> aerodynamic stiffness matrix K_A
% AC displacement: h_AC = w + e*theta
% Strip lift:      Lw'  = q*c*CL_alpha*theta
% We isolate the matrix K_A such that:
%   delta W_A^w = delta(phi)^T * q * K_A * phi

hAC_shape = [Nw, e*Nt];
Lw_shape  = c*CL_alpha*Nt;

% Matrix assembled by inspection from the derivation
KA = [ 0, simplify(int(Nw * c*CL_alpha*Nt, y, 0, b));
       0, simplify(int(e*Nt * c*CL_alpha*Nt, y, 0, b)) ];

KA1  = KA(1,2);
KA2  = KA(2,2);

fprintf('\nAerodynamic stiffness matrix:\n');
disp('K_A ='); disp(KA);

%% 5) Aileron aerodynamic virtual work -> forcing vector F_a
%   delta W_A^a = delta(phi)^T * q * beta * F_a
% Aileron span: [b-ba, b]

Fa1 = simplify((c/b) * CL_beta * (b^3 - (b-ba)^3) / (3*b));
Fa2 = simplify((c/b) * (e*CL_beta + c*Cm_beta) * (b^2 - (b-ba)^2) / 2);
Fa  = [Fa1; Fa2];

fprintf('\nAileron forcing vector:\n');
disp('F_a ='); disp(Fa);

%% 6) Strut equivalent stiffness
% Geometry:
%   ell = z_A / sin(gamma)
%   u   = h_B sin(gamma)
%   F_Sz = -(EA/z_A) sin^3(gamma) * h_B
% Define:
%   eta_B = z_A / (b tan(gamma))

KSTR = simplify((E*A/zA) * sin(gamma)^3);
etaB = simplify(zA / (b*tan(gamma)));

Kstrut_matrix = simplify(KSTR * etaB^2 * [etaB^2,      -etaB*xB;
                                          -etaB*xB,    xB^2   ]);

fprintf('\nStrut equivalent quantities:\n');
disp('K_STR ='); disp(KSTR);
disp('eta_B ='); disp(etaB);
disp('Strut contribution matrix ='); disp(Kstrut_matrix);

%% 7) Total structural matrix K_S
% From PVW:
%   (K_S - q K_A) phi = 0
% with
%   K_S = K_internal + Kstrut_matrix

KS = simplify(K_internal + Kstrut_matrix);

KS1  = simplify(KS(1,1));
KS12 = simplify(KS(1,2));
KS2  = simplify(KS(2,2));

fprintf('\nTotal structural matrix K_S =\n');
disp(KS);

fprintf('\nExplicit structural coefficients:\n');
disp('K_S1 =');  disp(KS1);
disp('K_S12 ='); disp(KS12);
disp('K_S2 =');  disp(KS2);

%% 8) Closed-form expressions exactly as in the report
KS1_closed  = simplify(4*E*J/b^3 + E*A*zA^3*cos(gamma)^4/(b^4*sin(gamma)));
KS12_closed = simplify(-E*A*zA^2*cos(gamma)^3*xB/b^3);
KS2_closed  = simplify(G*J/b + E*A*sin(gamma)*cos(gamma)^2*zA*xB^2/b^2);

fprintf('\nCheck with the compact closed forms from the derivation:\n');
disp('simplify(K_S1 - K_S1_closed) =');  disp(simplify(KS1 - KS1_closed));
disp('simplify(K_S12 - K_S12_closed) ='); disp(simplify(KS12 - KS12_closed));
disp('simplify(K_S2 - K_S2_closed) =');   disp(simplify(KS2 - KS2_closed));

%% 9) Divergence condition and q_D
% Divergence is obtained from:
%   det(K_S - q_D K_A) = 0

KTOT_qD = simplify(KS - qD*KA);
det_KTOT_qD = simplify(det(KTOT_qD));

fprintf('\nDeterminant for divergence condition det(K_S - q_D K_A):\n');
disp(det_KTOT_qD);

% Solve symbolically for q_D
qD_sol = simplify(solve(det_KTOT_qD == 0, qD, 'Real', true));

fprintf('\nSymbolic solution(s) for q_D:\n');
disp(qD_sol);

% Closed-form scalar expression from the derivation
qD_closed = simplify((KS1*KS2 - KS12^2) / (KS1*KA2 - KS12*KA1));

fprintf('\nClosed-form q_D used in the report:\n');
disp(qD_closed);

fprintf('\nCheck symbolic consistency simplify(qD_sol - qD_closed):\n');
if numel(qD_sol) == 1
    disp(simplify(qD_sol - qD_closed));
else
    disp('Multiple symbolic roots returned by solve; inspect qD_sol manually.');
end

%% 10) Forced static problem with aileron
% For q < q_D:
%   (K_S - q K_A) phi = q beta F_a
% Define K_TOT(q)

KTOT = simplify(KS - q*KA);
det_KTOT = simplify(det(KTOT));
KTOT_inv = simplify(inv(KTOT));

phi = simplify(q * beta * KTOT_inv * Fa);
phi_w_expr     = simplify(phi(1));
phi_theta_expr = simplify(phi(2));

fprintf('\nStatic response under aileron forcing:\n');
disp('phi_w ='); disp(phi_w_expr);
disp('phi_theta ='); disp(phi_theta_expr);

%% 11) Manual formula for phi_theta from adjugate matrix
% For a 2x2 matrix:
%   K_TOT^{-1} = 1/det(K_TOT) * [K22 -K12; -K21 K11]
% Therefore:
%   phi_theta = q beta / det(K_TOT) * (K11*F_a2 - K21*F_a1)

KT11 = simplify(KTOT(1,1));
KT21 = simplify(KTOT(2,1));
KT22 = simplify(KTOT(2,2)); %#ok<NASGU>

phi_theta_manual = simplify(q*beta/det_KTOT * (KT11*Fa2 - KT21*Fa1));

fprintf('\nCheck phi_theta formula consistency:\n');
disp(simplify(phi_theta_expr - phi_theta_manual));

%% 12) Lift decomposition and control effectiveness E_C
% Total lift:
%   L = L_R + L_E
% where
%   L_R = q c ba CL_beta beta
%   L_E = q c (b CL_alpha / 2) phi_theta

LR = simplify(q * c * ba * CL_beta * beta);
LE = simplify(q * c * (b*CL_alpha/2) * phi_theta_expr);
L  = simplify(LR + LE);

EC = simplify(L / LR);
EC_closed = simplify(1 + q * (b/ba) * (CL_alpha/CL_beta) * ...
    (KT11*Fa2 - KT21*Fa1) / (2*det_KTOT));

fprintf('\nLift terms:\n');
disp('L_R ='); disp(LR);
disp('L_E ='); disp(LE);
disp('L   ='); disp(L);

fprintf('\nControl effectiveness E_C:\n');
disp(EC);

fprintf('\nCheck E_C against compact formula:\n');
disp(simplify(EC - EC_closed));

%% 13) Pretty-printed summary
fprintf('\n================ SUMMARY ================\n');
fprintf('K_w     = %s\n', char(Kw));
fprintf('K_theta = %s\n', char(Ktheta));
fprintf('q_D     = %s\n', char(qD_closed));
fprintf('E_C(q)  = %s\n', char(EC_closed));

%% 14) Optional numerical evaluation
% Fill values below to perform a numerical study.
% Keep USE_NUMERIC = false if you only want symbolic derivations.

USE_NUMERIC = false;

if USE_NUMERIC
    % ---------------- USER VALUES ----------------
    b_val        = 10.0;
    ba_val       = 2.0;
    c_val        = 1.5;
    e_val        = 0.15;
    xB_val       = 0.20;
    zA_val       = 1.20;
    gamma_val    = deg2rad(35);

    E_val        = 70e9;
    G_val        = 27e9;
    A_val        = 3.0e-4;
    J_val        = 1.5e-4;

    CL_alpha_val = 5.7;
    CL_beta_val  = 3.2;
    Cm_beta_val  = -0.45;
    beta_val     = deg2rad(5);
    % ---------------------------------------------

    pars = [b, ba, c, e, xB, zA, gamma, E, G, A, J, CL_alpha, CL_beta, Cm_beta, beta];
    vals = [b_val, ba_val, c_val, e_val, xB_val, zA_val, gamma_val, E_val, G_val, A_val, J_val, CL_alpha_val, CL_beta_val, Cm_beta_val, beta_val];

    qD_num = double(subs(qD_closed, pars, vals));

    fprintf('\nNumerical divergence pressure q_D = %.6g\n', qD_num);

    % Evaluate E_C(q) below divergence
    q_vec = linspace(0, 0.95*qD_num, 300);
    EC_fun = matlabFunction(subs(EC_closed, pars, vals), 'Vars', q);
    EC_vec = EC_fun(q_vec);

    figure;
    plot(q_vec, EC_vec, 'LineWidth', 1.5);
    xlabel('q');
    ylabel('E_C');
    title('Control effectiveness versus dynamic pressure');
    grid on;

    % Optional display of matrices
    KS_num = double(subs(KS, pars(1:end-1), vals(1:end-1)));
    KA_num = double(subs(KA, [b,c,e,CL_alpha], [b_val,c_val,e_val,CL_alpha_val]));
    Fa_num = double(subs(Fa, pars, vals));

    disp('Numerical K_S ='); disp(KS_num);
    disp('Numerical K_A ='); disp(KA_num);
    disp('Numerical F_a ='); disp(Fa_num);
end

%% 15) Notes
% 1) The latex expression writes K_{A_1} in one place and K_{A_{12}} in another.
%    In the actual 2x2 matrix structure used here, the correct off-diagonal term is:
%       K_A(1,2) = K_{A1} = cb/4 * C_{L/alpha}
% 2) The denominator in q_D is therefore:
%       K_S1*K_A2 - K_S12*K_A1
% 3) If needed, this .m script can be converted manually into a Live Script (.mlx)
%    directly from MATLAB: open the file, then Save As -> Live Script.

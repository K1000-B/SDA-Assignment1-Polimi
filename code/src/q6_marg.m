%% ========================================================================
%  ITERATIVE REDESIGN FOR DIVERGENCE PRESSURE
%  ------------------------------------------------------------------------
%  Goal:
%  Find a light structural redesign such that
%       qD_new >= 1.4 * qD_baseline
%
%  Allowed design variables:
%       EJ, GJ, EA, gamma
%
%  Fixed parameters:
%       x_B, z_A, aerodynamic data
%
%  Method:
%  1) Start from the baseline design
%  2) Search locally on:
%         - 1 parameter at a time
%         - 2 parameters at a time
%         - 3 parameters at a time
%         - all 4 parameters together
%  3) Keep the lightest feasible candidate
%  4) Repeat until the mass index no longer improves
%
%  This code uses the same simple 1-mode model as in the rest of the report.
%  If you already have your own q_D function, you only need to replace
%  the local function "compute_qD" at the end.
%% ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1) DATA FROM THE STATEMENT
%% ------------------------------------------------------------------------
p.EJ0        = 1e7;        % [N.m^2]
p.GJ0        = 1e6;        % [N.m^2]
p.EA0        = 1e8;        % [N]
p.CL_alpha   = 2*pi;       % [-]
p.b          = 14;         % [m]
p.b_a        = 2.5;        % [m]
p.c          = 1;          % [m]
p.e          = p.c/4;      % [m]
p.gamma0_deg = 15;         % [deg]

p.x_B = -0.1;              % [m]  
p.z_A = 2;                 % [m] 



p.gamma0 = deg2rad(p.gamma0_deg);

%% ------------------------------------------------------------------------
%  2) BASELINE DESIGN
%% ------------------------------------------------------------------------
X0 = [p.EJ0, p.GJ0, p.EA0, p.gamma0];   % [EJ, GJ, EA, gamma]

qD0 = compute_qD(X0, p);

qD_target = 1.4 * qD0;

fprintf('\n============================================================\n');
fprintf('BASELINE DESIGN\n');
fprintf('============================================================\n');
fprintf('qD_baseline   = %.6g Pa\n', qD0);
fprintf('qD_target     = %.6g Pa   (= 1.4 * qD_baseline)\n', qD_target);
fprintf('x_B           = %.6g m\n', p.x_B);
fprintf('z_A           = %.6g m\n', p.z_A);
fprintf('gamma_0       = %.3f deg\n', rad2deg(p.gamma0));
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
%  3) MASS INDEX
%% ------------------------------------------------------------------------

function M = mass_index(X, p)
    EJ    = X(1);
    GJ    = X(2);
    EA    = X(3);
    gamma = X(4);

    M = (EJ / p.EJ0 - 1) ...
      + (GJ / p.GJ0 - 1) ...
      + (EA / p.EA0) * (sin(p.gamma0) / sin(gamma)) - 1;
end


%% ------------------------------------------------------------------------
%  4) SEARCH SETTINGS
%% ------------------------------------------------------------------------
% Lower bounds: we do not go below the baseline values
lb = [p.EJ0, p.GJ0, p.EA0, p.gamma0];

% Upper bounds
ub = [5.0*p.EJ0, ...
      5.0*p.GJ0, ...
      5.0*p.EA0, ...
      deg2rad(60)];

% Initial local steps
step.EJ    = 0.1 * p.EJ0;      
step.GJ    = 0.1 * p.GJ0;      
step.EA    = 0.1 * p.EA0;      
step.gamma = deg2rad(1);      

% Minimum steps for stopping
step_min.EJ    = 0.001 * p.EJ0; 
step_min.GJ    = 0.001 * p.GJ0; 
step_min.EA    = 0.001 * p.EA0; 
step_min.gamma = deg2rad(0.001); 

% Step update factors
grow_factor   = 1.25;   % used if no feasible design is found yet
shrink_factor = 0.1;   % used during refinement

% Small tolerance to decide if the mass really improved
mass_tol = 1e-14;

% Safety cap
max_outer_iter = 60;

%% ------------------------------------------------------------------------
%  5) ALL SEARCH FAMILIES
%  Each row says which variables are allowed to move:
%  [EJ  GJ  EA  gamma]
%% ------------------------------------------------------------------------
subset_list = [
    1 0 0 0
    0 1 0 0
    0 0 1 0
    0 0 0 1
    1 1 0 0
    1 0 1 0
    1 0 0 1
    0 1 1 0
    0 1 0 1
    0 0 1 1
    1 1 1 0
    1 1 0 1
    1 0 1 1
    0 1 1 1
    1 1 1 1
];

%% ------------------------------------------------------------------------
%  6) ITERATIVE SEARCH
%% ------------------------------------------------------------------------
Xcur = X0;

history_iter   = [];
history_mass   = [];
history_qD     = [];
history_EJ     = [];
history_GJ     = [];
history_EA     = [];
history_gamma  = [];
history_move   = strings(0);

for outer = 1:max_outer_iter

    qD_cur   = compute_qD(Xcur, p);
    mass_cur = mass_index(Xcur, p);

    is_current_feasible = isfinite(qD_cur) && (qD_cur >= qD_target);

    % ---------------------------------------------------------------------
    % Search mode:
    % - "expand" : baseline phase, we only try larger values
    % - "refine" : feasible phase, we allow up/down local moves
    % ---------------------------------------------------------------------
    if is_current_feasible
        search_mode = "refine";
    else
        search_mode = "expand";
    end

    best_candidate = [];
    best_mass      = inf;
    best_qD        = NaN;
    best_move_name = "";

    % If the current design is already feasible, it is our starting reference
    if is_current_feasible
        best_candidate = Xcur;
        best_mass      = mass_cur;
        best_qD        = qD_cur;
        best_move_name = "keep current design";
    end

    % ---------------------------------------------------------------------
    % Explore all search families: singletons, doublets, triplets, quadruplet
    % ---------------------------------------------------------------------
    for s = 1:size(subset_list,1)

        subset = logical(subset_list(s,:));
        candidates = build_candidates(Xcur, step, lb, ub, subset, search_mode);

        for k = 1:size(candidates,1)

            Xtest = candidates(k,:);
            qDtest = compute_qD(Xtest, p);

            % Only feasible designs are allowed
            if isfinite(qDtest) && (qDtest >= qD_target)

                mtest = mass_index(Xtest, p);

                if mtest < best_mass - mass_tol
                    best_candidate = Xtest;
                    best_mass      = mtest;
                    best_qD        = qDtest;
                    best_move_name = subset_name(subset);
                end
            end
        end
    end

    % ---------------------------------------------------------------------
    % Decision after all local searches
    % ---------------------------------------------------------------------
    if isempty(best_candidate)
        % No feasible design found yet -> enlarge the search window
        step = scale_steps(step, grow_factor);

        fprintf('Iter %02d | mode = %-6s | no feasible candidate found -> enlarge steps\n', ...
            outer, search_mode);

    else
        candidate_is_better = isempty(Xcur) || ...
                              (~is_current_feasible) || ...
                              (best_mass < mass_cur - mass_tol);

        if candidate_is_better
            Xcur = best_candidate;

            fprintf(['Iter %02d | mode = %-6s | best move = %-18s | ', ...
                     'qD/qD0 = %.4f | mass index = %.6f\n'], ...
                     outer, search_mode, best_move_name, best_qD/qD0, best_mass);

        else
            % Current design is feasible, but nothing lighter was found
            step = scale_steps(step, shrink_factor);

            fprintf(['Iter %02d | mode = %-6s | no lighter feasible design -> ', ...
                     'shrink steps\n'], outer, search_mode);
        end
    end

    % ---------------------------------------------------------------------
    % Save history after this outer iteration
    % ---------------------------------------------------------------------
    qD_cur   = compute_qD(Xcur, p);
    mass_cur = mass_index(Xcur, p);

    history_iter(end+1,1)  = outer;
    history_mass(end+1,1)  = mass_cur;
    history_qD(end+1,1)    = qD_cur;
    history_EJ(end+1,1)    = Xcur(1);
    history_GJ(end+1,1)    = Xcur(2);
    history_EA(end+1,1)    = Xcur(3);
    history_gamma(end+1,1) = Xcur(4);

    if isempty(best_move_name)
        history_move(end+1,1) = "none";
    else
        history_move(end+1,1) = best_move_name;
    end

    % ---------------------------------------------------------------------
    % Stopping criterion:
    % Stop if:
    % - current design is feasible
    % - and all steps are already very small
    % - and no lighter feasible design was found
    % ---------------------------------------------------------------------
    small_steps = ...
        (step.EJ    <= step_min.EJ)    && ...
        (step.GJ    <= step_min.GJ)    && ...
        (step.EA    <= step_min.EA)    && ...
        (step.gamma <= step_min.gamma);

    if is_current_feasible && small_steps && ~candidate_is_better
        fprintf('\nStopping criterion reached: mass converged.\n');
        break;
    end
end

%% ------------------------------------------------------------------------
%  7) FINAL RESULT
%% ------------------------------------------------------------------------
Xstar = Xcur;
qDstar = compute_qD(Xstar, p);
Mstar  = mass_index(Xstar, p);

fprintf('\n============================================================\n');
fprintf('FINAL DESIGN\n');
fprintf('============================================================\n');
fprintf('EJ*          = %.6g  N.m^2   (%.4f x baseline)\n', Xstar(1), Xstar(1)/p.EJ0);
fprintf('GJ*          = %.6g  N.m^2   (%.4f x baseline)\n', Xstar(2), Xstar(2)/p.GJ0);
fprintf('EA*          = %.6g  N       (%.4f x baseline)\n', Xstar(3), Xstar(3)/p.EA0);
fprintf('gamma*       = %.4f deg      (%.4f x baseline)\n', rad2deg(Xstar(4)), Xstar(4)/p.gamma0);
fprintf('qD*          = %.6g Pa\n', qDstar);
fprintf('qD*/qD0      = %.6f\n', qDstar/qD0);
fprintf('qD*/qDtarget = %.6f\n', qDstar/qD_target);
fprintf('Mass index   = %.6f\n', Mstar);
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
%  8) CLEAR SUMMARY TABLES
%% ------------------------------------------------------------------------
T_design = table( ...
    ["EJ"; "GJ"; "EA"; "gamma"], ...
    [p.EJ0; p.GJ0; p.EA0; rad2deg(p.gamma0)], ...
    [Xstar(1); Xstar(2); Xstar(3); rad2deg(Xstar(4))], ...
    [Xstar(1)/p.EJ0; Xstar(2)/p.GJ0; Xstar(3)/p.EA0; Xstar(4)/p.gamma0], ...
    'VariableNames', {'Parameter','Baseline','Final','Final_over_Baseline'});

T_metrics = table( ...
    qD0, qD_target, qDstar, qDstar/qD0, qDstar/qD_target, Mstar, ...
    'VariableNames', {'qD_baseline','qD_target','qD_final', ...
                      'qDfinal_over_qD0','qDfinal_over_qDtarget','MassIndex'});

disp('Design summary:');
disp(T_design);

disp('Global metrics:');
disp(T_metrics);

%% ------------------------------------------------------------------------
%  9) SIMPLE VISUAL OUTPUT FOR THE REPORT
%  Only one figure: convergence history
%% ------------------------------------------------------------------------
figure('Color','w','Position',[100 100 700 430]);

yyaxis left
plot(history_iter, history_mass, '-o', 'LineWidth', 1.6, 'MarkerSize', 5);
ylabel('Mass index');
grid on;
hold on;

yyaxis right
plot(history_iter, history_qD/qD0, '-s', 'LineWidth', 1.6, 'MarkerSize', 5);
yline(qD_target/qD0, '--', 'Target', 'LineWidth', 1.3);
ylabel('q_D / q_{D,0}');

xlabel('Outer iteration');
title('Iterative convergence');

txt = sprintf(['x_B = %.3f m\n', ...
               'z_A = %.3f m\n', ...
               'q_D^*/q_{D,0} = %.3f\n', ...
               'q_D^*/q_D^{target} = %.3f\n', ...
               'Mass index = %.3f'], ...
               p.x_B, p.z_A, qDstar/qD0, qDstar/qD_target, Mstar);

text(0.03, 0.95, txt, ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'Margin',8);

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function qD = compute_qD(X, p)
% Computes the divergence dynamic pressure for the 1-mode model.
%
% X = [EJ, GJ, EA, gamma]

    EJ    = X(1);
    GJ    = X(2);
    EA    = X(3);
    gamma = X(4);

    % Basic angle check
    if gamma <= 0 || gamma >= pi/2
        qD = NaN;
        return;
    end

        % Spanwise attachment position of the strut
    y_B = p.z_A / tan(gamma);
    eta_B = y_B / p.b;

    % ---------------------------------------------------------------------
    % Geometric admissibility condition:
    % the strut must remain at least 0.5 m before the start of the aileron
    %
    % Aileron starts at: y = b - b_a
    % So we impose:      y_B <= b - b_a - 0.5
    %
    % Equivalent nondimensional condition:
    % eta_B <= (b - b_a - 0.5) / b
    % ---------------------------------------------------------------------
    eta_B_max = (p.b - p.b_a - 0.5) / p.b;

    if y_B <= 0 || y_B >= p.b || eta_B > eta_B_max
        qD = NaN;
        return;
    end

    % 1-mode wing stiffnesses
    K_w     = 4 * EJ / p.b^3;
    K_theta = GJ / p.b;

    % Strut stiffness
    K_STR = (EA / p.z_A) * sin(gamma)^3;

    eta = y_B / p.b;

    % Structural matrix terms
    a = K_w     + K_STR * eta^4;
    b =         - K_STR * p.x_B * eta^3;
    c =         - K_STR * p.x_B * eta^3;
    d = K_theta + K_STR * p.x_B^2 * eta^2;

    % Aerodynamic matrix terms
    K_A1 = (p.c * p.b / 4) * p.CL_alpha;
    K_A2 = (p.e * p.c * p.b / 3) * p.CL_alpha;

    % Divergence condition:
    % det( K_struct - qD * K_aero ) = 0
    num = a*d - c*b;
    den = a*K_A2 - c*K_A1;

    if abs(den) < 1e-14
        qD = NaN;
        return;
    end

    qD = num / den;

    % Keep only positive, physical values
    if ~isfinite(qD) || qD <= 0
        qD = NaN;
    end
end



function candidates = build_candidates(Xcur, step, lb, ub, subset, mode)
% Builds a small local grid around the current design.
%
% mode = "expand" -> only forward moves [0, +1, +2]
% mode = "refine" -> local symmetric moves [-1, 0, +1]

    if mode == "expand"
        move_values = [0, 1, 2];
    else
        move_values = [-1, 0, 1];
    end

    levels = cell(1,4);
    for i = 1:4
        if subset(i)
            levels{i} = move_values;
        else
            levels{i} = 0;
        end
    end

    [L1,L2,L3,L4] = ndgrid(levels{1}, levels{2}, levels{3}, levels{4});
    moves = [L1(:), L2(:), L3(:), L4(:)];

    step_vec = [step.EJ, step.GJ, step.EA, step.gamma];
    candidates = zeros(size(moves,1),4);

    for k = 1:size(moves,1)
        X = Xcur + moves(k,:) .* step_vec;

        % Respect bounds
        X = max(X, lb);
        X = min(X, ub);

        candidates(k,:) = X;
    end

    % Remove duplicates that may appear because of bound clipping
    candidates = unique(candidates, 'rows');
end

function name = subset_name(subset)
% Returns a readable name for the active variables.

    labels = ["EJ", "GJ", "EA", "gamma"];
    name = strjoin(labels(subset), ' + ');
end

function step = scale_steps(step, factor)
% Scales all local steps by the same factor.

    step.EJ    = factor * step.EJ;
    step.GJ    = factor * step.GJ;
    step.EA    = factor * step.EA;
    step.gamma = factor * step.gamma;
end
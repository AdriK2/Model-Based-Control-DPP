format compact;

%% =======================
% Parameters
% =======================
Voc = 39.71 / 3;
Isc = 9.91;
Ns  = 60 / 3;

q = 1.602176634e-19;
k = 1.380649e-23;
T = 298.15;
A = 1;

Rs  = 0.0907;
Rsh = 276.4340;

alpha = q / (Ns * k * T * A);

Is = (Isc + Isc * Rs / Rsh - Voc / Rsh) / ...
     (exp(alpha * Voc) - exp(alpha * Isc * Rs));

%% =======================
% User Input
% =======================
Im = input('Enter Im (Current): ');

%% =======================
% Solver Settings
% =======================
tol = 1e-10;
max_iter = 500;

% Voltage scanning range
V_min = -1000;
V_max =  1000;
N_scan = 200000;

V_scan = linspace(V_min, V_max, N_scan);

%% =======================
% Define Residual Function F(Vm, Im)
% MPP condition: dP/dV = 0
%
% dP/dV = I + V*dI/dV = 0
% =======================
F_func = @(V) localResidual(Im, V, alpha, Is, Rs, Rsh);

%% =======================
% Step 1: Brutal Scan to Find Sign-Change Interval
% =======================
F_scan = arrayfun(F_func, V_scan);

idx_all = find(F_scan(1:end-1) .* F_scan(2:end) <= 0);

if isempty(idx_all)
    error(['No sign-change interval found. ' ...
           'Try increasing scan range or checking Im.']);
end

if length(idx_all) > 1
    warning(['More than one sign-change interval found. ' ...
             'Using the first one.']);
end

best_idx = idx_all(1);

Vm_lo = V_scan(best_idx);
Vm_hi = V_scan(best_idx + 1);

F_lo = F_func(Vm_lo);
F_hi = F_func(Vm_hi);

fprintf('\n--- Sign-change interval found ---\n');
fprintf('Vm_lo = %.12f, F_lo = %.12e\n', Vm_lo, F_lo);
fprintf('Vm_hi = %.12f, F_hi = %.12e\n', Vm_hi, F_hi);

%% =======================
% Step 2: Bisection Method
% =======================
fprintf('\n--- Start Bisection ---\n');

Vm_found = NaN;

for iter = 1:max_iter

    Vm_mid = 0.5 * (Vm_lo + Vm_hi);
    F_mid  = F_func(Vm_mid);

    fprintf('Iter %4d | Vm = %14.10f | F = %+14.6e\n', ...
        iter, Vm_mid, F_mid);

    % Convergence check
    if abs(F_mid) < tol || abs(Vm_hi - Vm_lo) < tol
        Vm_found = Vm_mid;

        fprintf('\n--- Target HIT after %d iterations ---\n', ...
            iter);
        break;
    end

    % Update bisection interval
    if sign(F_mid) == sign(F_lo)
        Vm_lo = Vm_mid;
        F_lo  = F_mid;
    else
        Vm_hi = Vm_mid;
        F_hi  = F_mid;
    end
end

if isnan(Vm_found)
    warning(['Did not fully converge within max_iter. ' ...
             'Using midpoint.']);

    Vm_found = 0.5 * (Vm_lo + Vm_hi);
end

%% =======================
% Final Output
% =======================
F_final = F_func(Vm_found);

fprintf('\n=== RESULT ===\n');
fprintf('Im input  = %.8f A\n', Im);
fprintf('Vm found  = %.12f V\n', Vm_found);
fprintf('F(Vm,Im)  = %.12e\n', F_final);
fprintf('|F|       = %.3e\n', abs(F_final));

%% =======================
% Local Function
% =======================
function F = localResidual(Im, Vm, alpha, Is, Rs, Rsh)

    arg = alpha * (Vm + Im * Rs);

    % Prevent numerical overflow
    arg = min(arg, 100);

    exp_term = exp(arg);

    % dI/dV from the single-diode PV model
    dIdV = (-Is * alpha * exp_term - 1/Rsh) / ...
           (1 + Rs/Rsh + alpha * Is * Rs * exp_term);

    % MPP condition:
    % dP/dV = I + V*dI/dV = 0
    F = Im + Vm * dIdV;

end
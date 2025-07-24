% Clear Workspace and Close Figures
clear; close all; clc;

%% Import Data
[data, ~, ~] = xlsread('ECE370-Lab 5.xlsx', 'Sheet1');

% Adjust column indices based on your data layout
Vbe1 = data(:, 1);   Vcc1 = data(:, 2);  Vce1 = data(:, 3);  Ic1 = data(:, 4);
Vbe2 = data(:, 5);   Vcc2 = data(:, 6);  Vce2 = data(:, 7);  Ic2 = data(:, 8);
Vbe3 = data(:, 9);   Vcc3 = data(:, 10); Vce3 = data(:, 11); Ic3 = data(:, 12);
Vbe4 = data(:, 13);  Vcc4 = data(:, 14); Vce4 = data(:, 15); Ic4 = data(:, 16);

% Assigned base currents (in μA)
Ib1 = 0; Ib2 = 10; Ib3 = 20; Ib4 = 30;

% Convert currents from μA to A
Ic1 = Ic1 * 1e-6; Ic2 = Ic2 * 1e-6; Ic3 = Ic3 * 1e-6; Ic4 = Ic4 * 1e-6;

% Convert VBE from mV to V if necessary
if max(Vbe1) > 2, Vbe1 = Vbe1 * 1e-3; end
if max(Vbe2) > 2, Vbe2 = Vbe2 * 1e-3; end
if max(Vbe3) > 2, Vbe3 = Vbe3 * 1e-3; end
if max(Vbe4) > 2, Vbe4 = Vbe4 * 1e-3; end

% Determine max supply voltage (assume all sets share Vcc)
Vcc = max([Vcc1; Vcc2; Vcc3; Vcc4]);

%% Plot VCE vs. IC for all IB (2.a)
figure;
hold on;
plot(Vce1, Ic1 * 1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib1) ' μA']);
plot(Vce2, Ic2 * 1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib2) ' μA']);
plot(Vce3, Ic3 * 1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib3) ' μA']);
plot(Vce4, Ic4 * 1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib4) ' μA']);
xlabel('V_{CE} (V)');
ylabel('I_C (\muA)');
title('Collector Characteristic Curves for Various I_B');
legend('show');
grid on;
hold off;

%% β, α, γ Calculations (2.c, 2.d)
% Choose Vce > 0.3V as active region
idx2 = Vce2 > 0.3; beta2 = mean(Ic2(idx2)) / (Ib2 * 1e-6);
idx3 = Vce3 > 0.3; beta3 = mean(Ic3(idx3)) / (Ib3 * 1e-6);
idx4 = Vce4 > 0.3; beta4 = mean(Ic4(idx4)) / (Ib4 * 1e-6);

% Calculate α = β/(β+1) and γ ≈ α
alpha2 = beta2 / (beta2 + 1);
alpha3 = beta3 / (beta3 + 1);
alpha4 = beta4 / (beta4 + 1);

disp('β, α, γ values:');
disp(['Ib = 10 μA, β = ' num2str(beta2) ', α = ' num2str(alpha2) ', γ ≈ ' num2str(alpha2)]);
disp(['Ib = 20 μA, β = ' num2str(beta3) ', α = ' num2str(alpha3) ', γ ≈ ' num2str(alpha3)]);
disp(['Ib = 30 μA, β = ' num2str(beta4) ', α = ' num2str(alpha4) ', γ ≈ ' num2str(alpha4)]);

%% Voltage Transfer Characteristic (2.e)
Rc = 1e3; % Resistor value in ohms
Is = 1e-15; % Initial guess for reverse saturation current
n = 1; % Ideality factor
k = 1.38064852e-23; T = 300; q = 1.60217662e-19;
Vt = k * T / q; % Thermal voltage ~ 25.85 mV

Vbe_calc = linspace(0.5, 0.8, 100);
Ic_calc = Is * exp(Vbe_calc / (n * Vt));
Vce_calc = Vcc - Rc * Ic_calc;

figure;
plot(Vbe_calc, Vce_calc, 'b-', 'LineWidth', 2);
xlabel('V_{BE} (V)');
ylabel('V_{CE} (V)');
title('Calculated VTC (V_{CE} vs. V_{BE})');
grid on;

%% Part i: Extract Is and n from Exponential Region
% Choose an exponential region (e.g., 0.55V < Vbe < 0.7V)
exp_idx = (Vbe2 > 0.55) & (Vbe2 < 0.7);
Vbe_exp = Vbe2(exp_idx);
Ic_exp = Ic2(exp_idx);

% Remove invalid points
valid_idx = Ic_exp > 0;
Vbe_exp = Vbe_exp(valid_idx);
Ic_exp = Ic_exp(valid_idx);

% Fit ln(Ic) vs. Vbe
if length(Vbe_exp) > 2 && all(Ic_exp > 0)
    [p, ~] = polyfit(Vbe_exp, log(Ic_exp), 1);
    slope = p(1);
    intercept = p(2);
    n_calculated = 1 / (slope * Vt);
    Is_calculated = exp(intercept);

    % Compare with measured Is at Vbe = 0
    Is_measured = Ic2(find(Vbe2 == 0, 1));
    disp(['Calculated Is: ' num2str(Is_calculated) ' A']);
    disp(['Measured Is (VBE=0): ' num2str(Is_measured) ' A']);
    disp(['Calculated Ideality Factor (n): ' num2str(n_calculated)]);
else
    disp('No valid exponential region for fit. Check data range.');
end

%% Part j: Compare Measured gm to Datasheet
% Calculate gm = Ic / (n * Vt) for a chosen point
Ic_sample = mean(Ic_exp);
if exist('n_calculated', 'var') && ~isnan(n_calculated) && Ic_sample > 0
    gm_measured = Ic_sample / (n_calculated * Vt);
    disp(['Measured Transconductance (g_m): ' num2str(gm_measured) ' S']);
else
    disp('Unable to calculate transconductance. Check data.');
end

% Compare with datasheet gm (replace with actual datasheet value)
gm_datasheet = 0.01; % Example
disp(['Datasheet Transconductance (g_m): ' num2str(gm_datasheet) ' S']);

%% Part k: Largest Signal Amplitude and Voltage Gain
% Identify linear region from VTC (e.g., 0.6V < Vbe < 0.7V)
linear_idx = (Vbe2 > 0.6) & (Vbe2 < 0.7);
Vbe_linear = Vbe2(linear_idx);
Vce_linear = Vce2(linear_idx);

% Calculate voltage gain A = ΔVout / ΔVin
if length(Vbe_linear) > 1
    dV = diff(Vce_linear) ./ diff(Vbe_linear);
    gain = mean(dV);
    disp(['Voltage Gain (A): ' num2str(gain)]);
else
    disp('No valid linear region for gain calculation. Check data.');
end

%% Revised Load Line Section
% Choose a larger resistor value to see the load line more clearly
Rc_load = 500e3;  % 500 kΩ for a visible load line in μA range
Vcc_load = 5;     % 5 V supply
Ic_line = linspace(0, 10e-6, 100);  % 0 to 10 μA range

% Load line equation: Vce = Vcc - Ic*Rc
Vce_line = Vcc_load - Ic_line * Rc_load;

% Plot the measured DCIV characteristic
figure;
plot(Vce2, Ic2 * 1e6, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Measured DCIV'); % Ic in μA
hold on;

% Plot the adjusted load line
plot(Vce_line, Ic_line * 1e6, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Load Line');

% Label the axes and add legend
xlabel('V_{CE} (V)');
ylabel('I_C (\muA)');
title('DCIV Characteristic with Load Line');
legend('show');
grid on;
hold off;

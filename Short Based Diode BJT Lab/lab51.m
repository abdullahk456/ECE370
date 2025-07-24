% Clear Workspace and Close Figures
clear; close all; clc;

% Read Data from Excel File
[data, txt, raw] = xlsread('ECE370-Lab 5.xlsx', 'Sheet1');

% Extract data for each dataset
% Adjust the column indices based on your Excel file
Vbe1 = data(:, 1);
Vcc1 = data(:, 2);
Vce1 = data(:, 3);
Ic1 = data(:, 4);

Vbe2 = data(:, 5);
Vcc2 = data(:, 6);
Vce2 = data(:, 7);
Ic2 = data(:, 8);

Vbe3 = data(:, 9);
Vcc3 = data(:, 10);
Vce3 = data(:, 11);
Ic3 = data(:, 12);

Vbe4 = data(:, 13);
Vcc4 = data(:, 14);
Vce4 = data(:, 15);
Ic4 = data(:, 16);

% Assign Ib values (in μA)
Ib1 = 0;    % Dataset 1
Ib2 = 10;   % Dataset 2
Ib3 = 20;   % Dataset 3
Ib4 = 30;   % Dataset 4

% Combine the data
Vce = [Vce1; Vce2; Vce3; Vce4];
Ic = [Ic1; Ic2; Ic3; Ic4];
Ib = [Ib1 * ones(size(Vce1)); Ib2 * ones(size(Vce2)); Ib3 * ones(size(Vce3)); Ib4 * ones(size(Vce4))];

% Convert currents to Amperes
Ib = Ib * 1e-6;   % Convert Ib from μA to A
Ic = Ic * 1e-6;   % Convert Ic from μA to A

% Define Vcc as the maximum supply voltage
Vcc = max([Vcc1; Vcc2; Vcc3; Vcc4]);

% Plot VCE against IC for all IB on a Single Plot
figure;
hold on;
Ib_values = unique(Ib) * 1e6;  % Convert to μA for labeling
colors = lines(length(Ib_values));  % Generate distinct colors for each Ib

for i = 1:length(Ib_values)
    idx = Ib == Ib_values(i) * 1e-6;  % Convert back to A for indexing
    plot(Vce(idx), Ic(idx) * 1e6, 'o-', 'Color', colors(i, :), 'DisplayName', ['I_B = ' num2str(Ib_values(i)) ' μA']);
end
xlabel('V_{CE} (V)');
ylabel('I_C (μA)');
title('Collector Characteristic Curves');
legend('show');
grid on;
hold off;

% Label the Cutoff, Active, and Saturation Regions
figure(1);
hold on;

% Cutoff Region
x_cutoff = [0 0.2];
y_cutoff = [0 Vcc];
patch([x_cutoff fliplr(x_cutoff)], [y_cutoff(1) y_cutoff(1) y_cutoff(2) y_cutoff(2)], 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
text(0.1, max(Ic) * 0.8e6, 'Cutoff Region', 'Color', 'red');

% Saturation Region
x_saturation = [0 0.3];
y_saturation = [0 Vcc];
patch([x_saturation fliplr(x_saturation)], [y_saturation(1) y_saturation(1) y_saturation(2) y_saturation(2)], 'green', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
text(0.15, max(Ic) * 0.6e6, 'Saturation Region', 'Color', 'green');

% Active Region
x_active = [0.3 max(Vce)];
y_active = [0 Vcc];
patch([x_active fliplr(x_active)], [y_active(1) y_active(1) y_active(2) y_active(2)], 'blue', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
text(mean(x_active), max(Ic) * 0.9e6, 'Active Region', 'Color', 'blue', 'HorizontalAlignment', 'center');

hold off;

% Calculate β values
beta_values = [];
Ib_unique = unique(Ib);
for i = 1:length(Ib_unique)
    idx = (Ib == Ib_unique(i)) & (Vce > 0.3);
    if sum(idx) > 0
        beta = mean(Ic(idx) ./ Ib_unique(i));
        beta_values = [beta_values; Ib_unique(i) * 1e6, beta];
    end
end

% Display β values
disp('β values for each Ib in the active region:');
disp('Ib (μA)     β');
disp(beta_values);

% Calculate α and γ for each IB
alpha_values = beta_values(:,2) ./ (beta_values(:,2) + 1);
gamma_values = alpha_values;  % Approximation

% Display α and γ values
disp('α and γ values for each Ib:');
disp('Ib (μA)     α          γ');
for i = 1:length(alpha_values)
    fprintf('%.3f       %.4f     %.4f\n', beta_values(i,1), alpha_values(i), gamma_values(i));
end

% Calculate the Voltage Transfer Characteristic (VTC) based on DCIV Characteristic
Rc = 1e3;  % Example value in Ohms
Is = 1e-15;  % Saturation current
n = 1;       % Ideality factor
Vt = 25.85e-3;  % Thermal voltage at room temperature

% Define the range of Vbe for plotting
Vbe = linspace(0.5, 0.8, 100);  % Adjusted range based on expected data
Ic_calculated = Is * exp(Vbe / (n * Vt));  % Calculate Ic using diode equation
Vce_calculated = Vcc - Rc * Ic_calculated;  % Calculate Vce using circuit equation

% Ensure Vbe and Vce_calculated have the same dimensions
if length(Vbe) ~= length(Vce_calculated)
    error('Dimension mismatch: Vbe and Vce_calculated must have the same length.');
end

% Plot the calculated VTC
figure;
plot(Vbe, Vce_calculated, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated VTC');
xlabel('V_{BE} (V)');
ylabel('V_{CE} (V)');
title('Calculated Voltage Transfer Characteristic (V_{CE} vs. V_{BE})');
grid on;
legend('show');


% e.ii Comment on how output voltage changes with a small change in input voltage in the linear region
% (Include comments in your report)

%% f. Plot the Measured Voltage Transfer Characteristic

% Assuming you have measured Vbe and corresponding Vce
% For example, use data columns for Vbe and Vce

% Extract Vbe and Vce from your data
Vbe = data(:,6);  % Adjust the column index based on your data
Vce_measured = data(:,4);  % Assuming Vce is in column 4

% Plot the measured VTC
figure;
plot(Vbe, Vce_measured, 'ro', 'DisplayName', 'Measured Data');
xlabel('V_{BE} (V)');
ylabel('V_{CE} (V)');
title('Measured Voltage Transfer Characteristic (V_{CE} vs. V_{BE})');
grid on;
legend('show');

% f.i Compare the calculated and measured voltage transfer characteristics
figure;
plot(Vbe, Vce_measured, 'ro', 'DisplayName', 'Measured Data');
hold on;
plot(Vbe, Vce_calculated, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated VTC');
xlabel('V_{BE} (V)');
ylabel('V_{CE} (V)');
title('Comparison of Measured and Calculated VTC');
grid on;
legend('show');
hold off;

%% g. Label the Regions of Operation on the Measured VTC

% Add annotations to the plot
figure(3);
hold on;

% Cutoff Region: Vbe below turn-on voltage (~0.6V)
x_cutoff = [min(Vbe) 0.6];
y_cutoff = [0 Vcc];
patch([x_cutoff fliplr(x_cutoff)], [Vcc Vcc 0 0], 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
text(0.55, Vcc*0.8, 'Cutoff Region', 'Color', 'red', 'HorizontalAlignment', 'right');

% Active Region: Vbe between 0.6V and saturation voltage
x_active = [0.6 0.7];
patch([x_active fliplr(x_active)], [Vcc Vcc 0 0], 'blue', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
text(0.65, Vcc*0.5, 'Active Region', 'Color', 'blue', 'HorizontalAlignment', 'center');

% Saturation Region: Vbe above 0.7V
x_saturation = [0.7 max(Vbe)];
patch([x_saturation fliplr(x_saturation)], [Vcc Vcc 0 0], 'green', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
text(0.75, Vcc*0.2, 'Saturation Region', 'Color', 'green');

hold off;

%% h. Plot the Measured Transfer Characteristic of VBE to IC

% Plot Vbe vs. Ic
figure;
semilogy(Vbe, Ic*1e3, 'ko-', 'LineWidth', 1.5);
xlabel('V_{BE} (V)');
ylabel('I_C (mA)');
title('Transfer Characteristic (I_C vs. V_{BE})');
grid on;

%% h.i Find the VBE Equivalent Series Resistance

% The equivalent series resistance (r_e) can be found from the slope of the Ic vs. Vbe curve
% In the exponential region (active region), the slope is:
% dIc/dVbe = Ic / (n * Vt)
% So r_e = n * Vt / Ic

% Let's calculate r_e at a specific Ic
Ic_point = 1e-3;  % 1 mA
re = n * Vt / Ic_point;
disp(['Equivalent Series Resistance (r_e) at I_C = ' num2str(Ic_point*1e3) ' mA: ' num2str(re) ' Ohms']);

%% h.iii Find the Turn-On Voltage

% The turn-on voltage is the Vbe at which Ic starts to increase significantly
% From the plot, estimate Vbe_turn_on
Vbe_turn_on = 0.6;  % Adjust based on where Ic starts increasing

disp(['Turn-On Voltage (V_{be}) is approximately: ' num2str(Vbe_turn_on) ' V']);

% h.iii.1 Is this the same as the datasheet value?
% Check the datasheet for the specified transistor and compare

%% h.iii Find the Reverse Saturation Current (Is)

% From the diode equation: Ic = Is * exp(Vbe / (n * Vt))
% Rearranged: Is = Ic / exp(Vbe / (n * Vt))

% Choose a point in the exponential region
Vbe_sample = 0.7;  % V
Ic_sample = interp1(Vbe, Ic, Vbe_sample);  % Interpolate Ic at Vbe_sample

Is_calculated = Ic_sample / exp(Vbe_sample / (n * Vt));
disp(['Calculated Reverse Saturation Current (Is): ' num2str(Is_calculated) ' A']);

%% h.iv Find the Ideality Factor (n)

% Using two points in the exponential region
Vbe1 = 0.65;
Vbe2 = 0.7;
Ic1 = interp1(Vbe, Ic, Vbe1);
Ic2 = interp1(Vbe, Ic, Vbe2);

% Using the diode equation:
% ln(Ic2 / Ic1) = (Vbe2 - Vbe1) / (n * Vt)
n_calculated = (Vbe2 - Vbe1) / (Vt * log(Ic2 / Ic1));

disp(['Calculated Ideality Factor (n): ' num2str(n_calculated)]);

%% h.v Compare the Plot to the DCIV Characteristic of a PN Junction

% The plot of Ic vs. Vbe on a semilog scale should be linear in the exponential region,
% similar to the I-V characteristic of a diode (PN junction).

% (Include comments and observations in your report)

%% h.vi Extract the Transconductance, gm, of the BJT

% Transconductance gm = dIc/dVbe = Ic / (n * Vt)
% Calculate gm at a specific Ic

gm = Ic_sample / (n_calculated * Vt);
disp(['Transconductance (g_m) at I_C = ' num2str(Ic_sample*1e3) ' mA: ' num2str(gm) ' S']);

% Alternatively, plot gm vs. Ic
gm_values = Ic ./ (n_calculated * Vt);

figure;
plot(Ic*1e3, gm_values, 'b-', 'LineWidth', 2);
xlabel('I_C (mA)');
ylabel('g_m (S)');
title('Transconductance (g_m) vs. Collector Current (I_C)');
grid on;

%% Notes:

% - Ensure that all units are consistent throughout the calculations.
% - Replace placeholder values with actual values from your data and circuit parameters.
% - For accurate calculations, precise values of Is, n, and Rc should be determined from your data.
% - Include comments, observations, and comparisons with datasheet values in your report.


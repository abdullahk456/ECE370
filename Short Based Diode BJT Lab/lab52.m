% Clear Workspace and Close Figures
clear; close all; clc;

%% Import Data
[data, txt, raw] = xlsread('ECE370-Lab 5.xlsx', 'Sheet1');

% Adjust these column indices based on your actual data layout
% Example assumption:
% Ib=0 μA:   Vbe1=col1, Vcc1=col2, Vce1=col3, Ic1=col4
% Ib=10 μA:  Vbe2=col5, Vcc2=col6, Vce2=col7, Ic2=col8
% Ib=20 μA:  Vbe3=col9, Vcc3=col10,Vce3=col11,Ic3=col12
% Ib=30 μA:  Vbe4=col13,Vcc4=col14,Vce4=col15,Ic4=col16

Vbe1 = data(:,1);   Vcc1 = data(:,2);  Vce1 = data(:,3);  Ic1 = data(:,4);
Vbe2 = data(:,5);   Vcc2 = data(:,6);  Vce2 = data(:,7);  Ic2 = data(:,8);
Vbe3 = data(:,9);   Vcc3 = data(:,10); Vce3 = data(:,11); Ic3 = data(:,12);
Vbe4 = data(:,13);  Vcc4 = data(:,14); Vce4 = data(:,15); Ic4 = data(:,16);

% Assigned base currents (in μA)
Ib1 = 0;   
Ib2 = 10;  
Ib3 = 20;  
Ib4 = 30;

% Convert Ic from μA to A
Ic1 = Ic1 * 1e-6; Ic2 = Ic2 * 1e-6; Ic3 = Ic3 * 1e-6; Ic4 = Ic4 * 1e-6;

% Check if VBE is in mV instead of V (if max > 2, likely mV)
% Adjust each dataset as needed
if max(Vbe1) > 2, Vbe1 = Vbe1*1e-3; end
if max(Vbe2) > 2, Vbe2 = Vbe2*1e-3; end
if max(Vbe3) > 2, Vbe3 = Vbe3*1e-3; end
if max(Vbe4) > 2, Vbe4 = Vbe4*1e-3; end

% Determine max supply voltage (assuming all sets have the same Vcc)
Vcc = max([Vcc1; Vcc2; Vcc3; Vcc4]);

%% Plot VCE vs. IC for all IB (2.a)
figure;
hold on;
plot(Vce1, Ic1*1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib1) ' μA']);
plot(Vce2, Ic2*1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib2) ' μA']);
plot(Vce3, Ic3*1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib3) ' μA']);
plot(Vce4, Ic4*1e6, 'o-', 'DisplayName', ['I_B = ' num2str(Ib4) ' μA']);
xlabel('V_{CE} (V)');
ylabel('I_C (μA)');
title('Collector Characteristic Curves for Various I_B');
legend('show');
grid on;
hold off;

%% Label regions (2.b) - approximate
figure(1);
hold on;
ymax = max([Ic1; Ic2; Ic3; Ic4])*1e6;
text(0.1, ymax*0.8, 'Cutoff Region (low VBE)', 'Color', 'red');
text(2, ymax*0.8, 'Active Region (0.6V < VBE < ~0.7V)', 'Color', 'blue');
text(1, ymax*0.3, 'Saturation Region (low VCE)', 'Color', 'green');
hold off;

%% Calculate β in the active region (2.c)
% Choose Vce > 0.3 V as active region
idx2 = Vce2 > 0.3; beta2 = mean(Ic2(idx2)) / (Ib2*1e-6);
idx3 = Vce3 > 0.3; beta3 = mean(Ic3(idx3)) / (Ib3*1e-6);
idx4 = Vce4 > 0.3; beta4 = mean(Ic4(idx4)) / (Ib4*1e-6);

disp('β values in the active region:');
disp(['Ib = 10 μA, β = ' num2str(beta2)]);
disp(['Ib = 20 μA, β = ' num2str(beta3)]);
disp(['Ib = 30 μA, β = ' num2str(beta4)]);

%% Calculate α and γ (2.d)
% α = β/(β+1)
alpha2 = beta2/(beta2+1);
alpha3 = beta3/(beta3+1);
alpha4 = beta4/(beta4+1);

disp('α and γ values:');
disp(['For Ib=10 μA: α = ' num2str(alpha2) ', γ ≈ ' num2str(alpha2)]);
disp(['For Ib=20 μA: α = ' num2str(alpha3) ', γ ≈ ' num2str(alpha3)]);
disp(['For Ib=30 μA: α = ' num2str(alpha4) ', γ ≈ ' num2str(alpha4)]);

%% Voltage Transfer Characteristic (2.e)
% Use a known Rc, Is, n, and Vt
Rc = 1e3;       % Adjust based on your circuit
Is = 1e-15;     % Initial guess
n = 1;          
k = 1.38064852e-23; T=300; q=1.60217662e-19;
Vt = k*T/q;     % Thermal voltage ~ 25.85mV at room temp

Vbe_calc = linspace(0.5,0.8,100);
Ic_calc = Is * exp(Vbe_calc/(n*Vt));
Vce_calc = Vcc - Rc*Ic_calc;

figure;
plot(Vbe_calc, Vce_calc, 'b-', 'LineWidth', 2);
xlabel('V_{BE} (V)');
ylabel('V_{CE} (V)');
title('Calculated VTC (V_{CE} vs. V_{BE})');
grid on;

%% Measured VTC (2.f)
% Choose one dataset as measured VTC, e.g., Ib=10 μA dataset
Vbe_measured = Vbe2;
Vce_measured = Vce2;
Ic_measured = Ic2;

% Clean data: remove negative or zero currents/voltages if any
valid_idx = isfinite(Vbe_measured) & isfinite(Vce_measured) & isfinite(Ic_measured) & Ic_measured>0;
Vbe_measured = Vbe_measured(valid_idx);
Vce_measured = Vce_measured(valid_idx);
Ic_measured = Ic_measured(valid_idx);

figure;
plot(Vbe_measured, Vce_measured, 'ro-', 'DisplayName','Measured Data');
xlabel('V_{BE}(V)');
ylabel('V_{CE}(V)');
title('Measured VTC (V_{CE} vs. V_{BE})');
grid on;
legend('show');

% Compare measured and calculated (2.f.i)
figure;
plot(Vbe_measured, Vce_measured, 'ro', 'DisplayName','Measured VTC');
hold on;
plot(Vbe_calc, Vce_calc, 'b-', 'LineWidth',2, 'DisplayName','Calculated VTC');
xlabel('V_{BE}(V)');
ylabel('V_{CE}(V)');
title('Comparison of Measured and Calculated VTC');
grid on;
legend('show');
hold off;

%% Label Regions on Measured VTC (2.g)
figure;
plot(Vbe_measured, Vce_measured, 'ro-', 'DisplayName', 'Measured VTC');
hold on;
text(0.55, mean(Vce_measured)*0.8, 'Cutoff Region (~VBE<0.6V)', 'Color','red');
text(0.65, mean(Vce_measured)*0.5, 'Active Region', 'Color','blue');
text(0.7, mean(Vce_measured)*0.2, 'Saturation Region (low VCE)', 'Color','green');
xlabel('V_{BE}(V)');
ylabel('V_{CE}(V)');
title('Measured VTC with Regions Annotated');
grid on;
legend('show');
hold off;

%% Plot the Transfer Characteristic I_C vs. V_BE on semilog (2.h)
figure;
semilogy(Vbe_measured, Ic_measured*1e6, 'ko-', 'LineWidth',1.5);
xlabel('V_{BE}(V)');
ylabel('I_C (μA)');
title('Transfer Characteristic (I_C vs. V_{BE}) on semilog');
grid on;

%% h.i Equivalent Series Resistance (r_e)
Ic_point = 1e-3; % 1mA
re = n*Vt/Ic_point;
disp(['Equivalent Series Resistance (r_e) at Ic=1mA: ' num2str(re) ' Ohms']);

%% h.iii Turn-On Voltage
% Approximate where Ic starts to rise significantly
% For a small-signal BJT, ~0.6-0.7V is typical
Vbe_turn_on = 0.6;
disp(['Turn-On Voltage ~' num2str(Vbe_turn_on) ' V']);

%% Extract Ideality Factor (n) and Is
% Select an exponential region: e.g., 0.58V < VBE < 0.65V
exp_idx = (Vbe_measured > 0.58) & (Vbe_measured < 0.65);
Vbe_exp = Vbe_measured(exp_idx);
Ic_exp = Ic_measured(exp_idx);

if length(Vbe_exp) > 2 && all(Ic_exp > 0)
    % Fit ln(Ic) = ln(Is) + (Vbe/(nVt))
    [p,S,mu] = polyfit(Vbe_exp, log(Ic_exp), 1);
    slope = p(1);
    intercept = p(2);
    
    n_extracted = 1/(slope*Vt);
    Is_extracted = exp(intercept);
    
    disp(['Extracted Ideality Factor n: ' num2str(n_extracted)]);
    disp(['Extracted Reverse Saturation Current Is: ' num2str(Is_extracted) ' A']);
else
    warning('No valid exponential region found. Check Vbe range or data quality.');
    n_extracted = NaN;
    Is_extracted = NaN;
end

%% h.vi Transconductance gm
% gm ~ Ic/(n*Vt) at a chosen Ic in the exponential region
if ~isnan(n_extracted) && ~isempty(Ic_exp)
    Ic_sample = mean(Ic_exp);
    gm = Ic_sample/(n_extracted*Vt);
    disp(['Transconductance (g_m) at chosen Ic: ' num2str(gm) ' S']);
else
    disp('Could not compute g_m due to invalid n or empty exponential region.');
end

%% (2.k) Gain Approximation
% Check if we have a range for linear operation in VBE:
linear_idx = (Vbe_measured > 0.6) & (Vbe_measured < 0.7);
if any(linear_idx)
    % Compute derivative dVce/dVbe:
    dV = diff(Vce_measured)./diff(Vbe_measured);
    % For simplicity, pick the corresponding subset of dV:
    % Need to match indices for diff:
    Vbe_mid = Vbe_measured(1:end-1)+diff(Vbe_measured)/2;
    linear_idx_mid = (Vbe_mid > 0.6) & (Vbe_mid < 0.7);
    
    if any(linear_idx_mid)
        gain_approx = mean(dV(linear_idx_mid));
        disp(['Approximate Voltage Gain A: ' num2str(gain_approx)]);
    else
        disp('No suitable midpoint data in 0.6-0.7V for gain calculation.');
    end
else
    disp('No data in 0.6-0.7V range for gain calculation.');
end

disp('Analysis complete. Check plots and console output.');

%EGH451 Assessment 1
%Author: Ethan Yun Sang CHAN
%Student ID: n12622401

clear all; close all; clc

%% Module parameters (datasheet STC)
Vmpp = 37.6;      % V at MPP
Voc  = 45.3;      % V open circuit
Impp = 5.32;      % A at MPP
Isc  = 5.84;      % A short circuit
ki   = 0.0373;    % Isc temp coefficient (A/K)
kv   = -0.2810;   % Voc temp coefficient (V/K)

Tstc = 25 + 273.15;   % STC temperature (K)
Gstc = 1000;          % STC irradiance (W/m²)

% Array configuration
ns = 16;   % number of modules in series per string
np = 2;    % number of parallel strings

% Other operating conditions (module level)
G_600 = 600;           % 600 W/m²
G_200 = 200;           % 200 W/m²
G_1000 = 1000;         %1000 W/m²
T_5 = 5 + 273.15;   
T_45 = 45 + 273.15;  

%% ----- 1. Exact four‑parameter extraction at STC -----
function [Iph, Vt, Rs, Io] = FourParaSTC(Impp, Vmpp, Isc, Voc)
    % Exact extraction for the four‑parameter model (no shunt)
    Iph = Isc;                       % standard assumption
    delta = Isc - Impp;
    if delta <= 0
        error('Impp must be less than Isc');
    end
    % Vt from derivative condition at MPP
    Vt = (2*Vmpp - Voc) / (log(delta/Isc) + Impp/delta);
    % Rs from the MPP equation
    Rs = (Vt * log(delta/Isc) + Voc - Vmpp) / Impp;
    % Io from Voc condition
    Io = Isc / exp(Voc/Vt);
end

[Iph, Vt, Rs, Io] = FourParaSTC(Impp, Vmpp, Isc, Voc)

%% ----- 2. Function to correct module parameters for T and G -----
function [Iscgt, Vocgt, Voct, Vtt, Iot] = FourParaCorr(Isc, Voc, ki, kv, T, Tstc, Vt, G, Gstc)
    % Returns corrected module parameters for given T and G
    Isct = Isc + ki * (T - Tstc);
    Voct = Voc + kv * (T - Tstc);
    Vtt = Vt * (T / Tstc);               % thermal voltage scales with T
    Iot = Isct / exp(Voct / Vtt);        % saturation current
    Iscgt = Isct * (G / Gstc);           % photo current scales with G
    Vocgt = log(Iscgt / Iot) * Vtt;      % open‑circuit voltage
end

%% ----- 3. Compute module I‑V curves for three irradiances (T = 25°C) -----
points = 500;   % number of points for smooth curves

% STC (1000 W/m², 25°C)
V_stc = linspace(0, Voc, points);
I_stc = zeros(1, points);
for i = 1:points
    V = V_stc(i);
    I = Iph;   % initial guess
    for iter = 1:10
        I = Iph - Io * (exp((V + I*Rs)/Vt) - 1);
    end
    I_stc(i) = I;
end

% 600 W/m², 25°C
[Iscgt_600G, Vocgt_600G, Voct_600G, Vtt_600G, Iot_600G] = ...
    FourParaCorr(Isc, Voc, ki, kv, Tstc, Tstc, Vt, G_600, Gstc);
V_600G = linspace(0, Vocgt_600G, points);
I_600G = zeros(1, points);
for i = 1:points
    V = V_600G(i);
    I = Iscgt_600G;
    for iter = 1:10
        I = Iscgt_600G - Iot_600G * (exp((V + I*Rs)/Vtt_600G) - 1);
    end
    I_600G(i) = I;
end

% 200 W/m², 25°C
[Iscgt_200G, Vocgt_200G, Voct_200G, Vtt_200G, Iot_200G] = ...
    FourParaCorr(Isc, Voc, ki, kv, Tstc, Tstc, Vt, G_200, Gstc);
V_200G = linspace(0, Vocgt_200G, points);
I_200G = zeros(1, points);
for i = 1:points
    V = V_200G(i);
    I = Iscgt_200G;
    for iter = 1:10
        I = Iscgt_200G - Iot_200G * (exp((V + I*Rs)/Vtt_200G) - 1);
    end
    I_200G(i) = I;
end

%% ----- 4. Compute module MPP (before scaling) -----
% STC module
P_mod_stc = V_stc .* I_stc;
[Pmax_mod_stc, idx_mod_stc] = max(P_mod_stc);
Vmpp_mod_stc = V_stc(idx_mod_stc);
Impp_mod_stc = I_stc(idx_mod_stc);

% 600 W/m² module
P_mod_600G = V_600G .* I_600G;
[Pmax_mod_600G, idx_mod_600G] = max(P_mod_600G);
Vmpp_mod_600G = V_600G(idx_mod_600G);
Impp_mod_600G = I_600G(idx_mod_600G);

% 200 W/m² module
P_mod_200G = V_200G .* I_200G;
[Pmax_mod_200G, idx_mod_200G] = max(P_mod_200G);
Vmpp_mod_200G = V_200G(idx_mod_200G);
Impp_mod_200G = I_200G(idx_mod_200G);

%% ----- 5. Compute Array I‑V curves for three Conditions (STC), (T = 5°C, 1000 W/m2) & (T = 45°C, 200 W/m2) -----

% STC (1000 W/m², 25°C)
V_stc_array = linspace(0, Voc, points);
I_stc_array = zeros(1, points);
for i = 1:points
    V = V_stc_array(i);
    I = Iph;   % initial guess
    for iter = 1:10
        I = Iph - Io * (exp((V + I*Rs)/Vt) - 1);
    end
    I_stc_array(i) = I;
end

% 1000 W/m², 5°C
[Iscgt_1000G_5C, Vocgt_1000G_5C, Voct_1000G_5C, Vtt_1000G_5C, Iot_1000G_5C] = ...
    FourParaCorr(Isc, Voc, ki, kv, T_5, Tstc, Vt, G_1000, Gstc);
V_1000G_5C = linspace(0, Vocgt_1000G_5C, points);
I_1000G_5C = zeros(1, points);
for i = 1:points
    V = V_1000G_5C(i);
    I = Iscgt_1000G_5C;
    for iter = 1:10
        I = Iscgt_1000G_5C - Iot_1000G_5C * (exp((V + I*Rs)/Vtt_1000G_5C) - 1);
    end
    I_1000G_5C(i) = I;
end

% 200 W/m², 45°C
[Iscgt_200G_45C, Vocgt_200G_45C, Voct_200G_45C, Vtt_200G_45C, Iot_200G_45C] = ...
    FourParaCorr(Isc, Voc, ki, kv, T_45, Tstc, Vt, G_200, Gstc);
V_200G_45C = linspace(0, Vocgt_200G_45C, points);
I_200G_45C = zeros(1, points);
for i = 1:points
    V = V_200G_45C(i);
    I = Iscgt_200G_45C;
    for iter = 1:10
        I = Iscgt_200G_45C - Iot_200G_45C * (exp((V + I*Rs)/Vtt_200G_45C) - 1);
    end
    I_200G_45C(i) = I;
end

%% ----- 6. Compute Array MPP for three Conditions -----
% STC module
P_1000G_5C = V_1000G_5C .* I_1000G_5C;
[Pmax_1000G_5C, idx_1000G_5C] = max(P_1000G_5C);
Vmpp_1000G_5C = V_1000G_5C(idx_1000G_5C);
Impp_1000G_5C = I_1000G_5C(idx_1000G_5C);

% 200 W/m², 45°C PV Panel
P_200G_45C = V_200G_45C .* I_200G_45C;
[Pmax_200G_45C, idx_200G_45C] = max(P_200G_45C);
Vmpp_200G_45C = V_200G_45C(idx_200G_45C);
Impp_200G_45C = I_200G_45C(idx_200G_45C);

%% ----- 7. Scale to array level (ns series, np parallel) -----
V_array_stc  = V_stc * ns;
I_array_stc  = I_stc * np;
P_array_stc  = V_array_stc .* I_array_stc;

V_array_600G = V_600G * ns;
I_array_600G = I_600G * np;
P_array_600G = V_array_600G .* I_array_600G;

V_array_200G = V_200G * ns;
I_array_200G = I_200G * np;
P_array_200G = V_array_200G .* I_array_200G;

V_array_1000G_5C = V_1000G_5C * ns;
I_array_1000G_5C = I_1000G_5C * np;
P_array_1000G_5C = V_array_1000G_5C .* I_array_1000G_5C;

V_array_200G_45C = V_200G_45C * ns;
I_array_200G_45C = I_200G_45C * np;
P_array_200G_45C = V_array_200G_45C .* I_array_200G_45C;

%% ----- 8. Find array MPP -----
[Pmax_arr_stc, idx_arr_stc] = max(P_array_stc);
Vmpp_arr_stc = V_array_stc(idx_arr_stc);
Impp_arr_stc = I_array_stc(idx_arr_stc);
Pmax_arr_stc = P_array_stc(idx_arr_stc);

[Pmax_arr_600G, idx_arr_600G] = max(P_array_600G);
Vmpp_arr_600G = V_array_600G(idx_arr_600G);
Impp_arr_600G = I_array_600G(idx_arr_600G);
Pmax_arr_600G = P_array_600G(idx_arr_600G);

[Pmax_arr_200G, idx_arr_200G] = max(P_array_200G);
Vmpp_arr_200G = V_array_200G(idx_arr_200G);
Impp_arr_200G = I_array_200G(idx_arr_200G);
Pmax_arr_200G = P_array_200G(idx_arr_200G);

[Pmax_arr_1000G_5C, idx_arr_1000G_5C] = max(P_array_1000G_5C);
Vmpp_arr_1000G_5C = V_array_1000G_5C(idx_arr_1000G_5C);
Impp_arr_1000G_5C = I_array_1000G_5C(idx_arr_1000G_5C);
Pmax_arr_1000G_5C = P_array_1000G_5C(idx_arr_1000G_5C);

[Pmax_arr_200G_45C, idx_arr_200G_45C] = max(P_array_200G_45C);
Vmpp_arr_200G_45C = V_array_200G_45C(idx_arr_200G_45C);
Impp_arr_200G_45C = I_array_200G_45C(idx_arr_200G_45C);
Pmax_arr_200G_45C = P_array_200G_45C(idx_arr_200G_45C);

%% ----- 9. Combined Module‑level I‑V & P‑V plot (top/bottom subplots) -----
figure('Name','Ameresco Solar 200J-V PV Module Characteristics','NumberTitle','off','Color','white','Position',[100,100,1000,900]);

% Top subplot: I‑V
subplot(2,1,1);
hold on; box on; grid on;
plot(V_stc, I_stc, 'r', 'LineWidth',1.5, 'DisplayName',sprintf('STC (1000 W/m², 25°C)'));
plot(V_600G, I_600G, 'b', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 25°C',G_600));
plot(V_200G, I_200G, 'm', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 25°C',G_200));
xlabel('Voltage (V)','FontSize',14,'FontWeight','bold');
ylabel('Current (A)','FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',12);
title('Module I‑V Characteristics','FontSize',14,'FontWeight','bold');
xlim([0, 46]); ylim([0, 8]);

% Annotations for I‑V subplot
xLim = xlim; yLim = ylim; xRange = diff(xLim); yRange = diff(yLim);
% STC (red)
yline(Impp_mod_stc, '--r', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_mod_stc + 0.02*xRange, Impp_mod_stc + 0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.2f A', Vmpp_mod_stc, Impp_mod_stc), ...
    'Color','r','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','r','Margin',5,'HandleVisibility','off');
% 600 W/m² (blue)
yline(Impp_mod_600G, '--b', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_mod_600G + 0.03*xRange, Impp_mod_600G + 0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.2f A', Vmpp_mod_600G, Impp_mod_600G), ...
    'Color','b','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','b','Margin',5,'HandleVisibility','off');
% 200 W/m² (magenta)
yline(Impp_mod_200G, '--m', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_mod_200G + 0.03*xRange, Impp_mod_200G + 0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.2f A', Vmpp_mod_200G, Impp_mod_200G), ...
    'Color','m','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','m','Margin',5,'HandleVisibility','off');
set(gca,'FontSize',12);
hold off;

% Bottom subplot: P‑V
subplot(2,1,2);
hold on; box on; grid on;
plot(V_stc, P_mod_stc, 'r', 'LineWidth',1.5, 'DisplayName',sprintf('STC (1000 W/m², 25°C)'));
plot(V_600G, P_mod_600G, 'b', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 25°C',G_600));
plot(V_200G, P_mod_200G, 'm', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 25°C',G_200));
xlabel('Voltage (V)','FontSize',14,'FontWeight','bold');
ylabel('Power (W)','FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',12);
title('Module P‑V Characteristics','FontSize',14,'FontWeight','bold');
xlim([0, 46]); ylim([0, 290]);

% Annotations for P‑V subplot
xLim = xlim; yLim = ylim; xRange = diff(xLim); yRange = diff(yLim);
% STC (red)
yline(Pmax_mod_stc, '--r', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_mod_stc + 0.03*xRange, Pmax_mod_stc + 0.08*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_mod_stc, Pmax_mod_stc), ...
    'Color','r','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','r','Margin',5,'HandleVisibility','off');
% 600 W/m² (blue)
yline(Pmax_mod_600G, '--b', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_mod_600G + 0.03*xRange, Pmax_mod_600G + 0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_mod_600G, Pmax_mod_600G), ...
    'Color','b','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','b','Margin',5,'HandleVisibility','off');
% 200 W/m² (magenta)
yline(Pmax_mod_200G, '--m', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_mod_200G + 0.03*xRange, Pmax_mod_200G + 0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_mod_200G, Pmax_mod_200G), ...
    'Color','m','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','m','Margin',5,'HandleVisibility','off');
set(gca,'FontSize',12);
hold off;

%% ----- 10. Combined Array‑level I‑V & P‑V plot (top/bottom subplots) -----
figure('Name','Ameresco Solar 200J-V PV Array Characteristics','NumberTitle','off','Color','white','Position',[100,100,1000,900]);

% Top subplot: Array I‑V
subplot(2,1,1);
hold on; box on; grid on;
plot(V_array_stc, I_array_stc, 'r', 'LineWidth',1.5, 'DisplayName',sprintf('STC (1000 W/m², 25°C)'));
plot(V_array_1000G_5C, I_array_1000G_5C, 'b', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 5°C',G_1000));
plot(V_array_200G_45C, I_array_200G_45C, 'm', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 45°C',G_200));
xlabel('Voltage (V)','FontSize',14,'FontWeight','bold');
ylabel('Current (A)','FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',12);
title(sprintf('Array I‑V Characteristics (%dS × %dP, %.1f kWp)', ns, np, Pmax_arr_stc/1000), ...
      'FontSize',14,'FontWeight','bold');
xlim([0, 900]); ylim([0, 16]);

% Annotations for array I‑V
xLim = xlim; yLim = ylim; xRange = diff(xLim); yRange = diff(yLim);
% STC (red)
yline(Impp_arr_stc, '--r', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_arr_stc+0.01*xRange, Impp_arr_stc+0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.1f A', Vmpp_arr_stc, Impp_arr_stc), ...
    'Color','r','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','r','Margin',5,'HandleVisibility','off');
% 1000 W/m², 5°C (blue)
yline(Impp_arr_1000G_5C, '--b', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_arr_1000G_5C+0.10*xRange, Impp_arr_1000G_5C-0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.1f A', Vmpp_arr_1000G_5C, Impp_arr_1000G_5C), ...
    'Color','b','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','b','Margin',5,'HandleVisibility','off');
% 200 W/m², 45°C (magenta)
yline(Impp_arr_200G_45C, '--m', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_arr_200G_45C+0.03*xRange, Impp_arr_200G_45C+0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.1f A', Vmpp_arr_200G_45C, Impp_arr_200G_45C), ...
    'Color','m','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','m','Margin',5,'HandleVisibility','off');
set(gca,'FontSize',12);
hold off;

% Bottom subplot: Array P‑V
subplot(2,1,2);
hold on; box on; grid on;
plot(V_array_stc, P_array_stc, 'r', 'LineWidth',1.5, 'DisplayName',sprintf('STC (1000 W/m², 25°C)'));
plot(V_array_1000G_5C, P_array_1000G_5C, 'b', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 5°C',G_1000));
plot(V_array_200G_45C, P_array_200G_45C, 'm', 'LineWidth',1.5, 'DisplayName',sprintf('%.0f W/m², 45°C',G_200));
xlabel('Voltage (V)','FontSize',14,'FontWeight','bold');
ylabel('Power (W)','FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',12);
title(sprintf('Array P‑V Characteristics (%dS × %dP, %.1f kWp)', ns, np, Pmax_arr_stc/1000), ...
      'FontSize',14,'FontWeight','bold');
xlim([0, 900]); ylim([0, 10000]);

% Annotations for array P‑V
xLim = xlim; yLim = ylim; xRange = diff(xLim); yRange = diff(yLim);
% STC (red)
yline(Pmax_arr_stc, '--r', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_arr_stc-0.12*xRange, Pmax_arr_stc+0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_arr_stc, Pmax_arr_stc), ...
    'Color','r','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','r','Margin',5,'HandleVisibility','off');
% 1000 W/m², 5°C (blue)
yline(Pmax_arr_1000G_5C, '--b', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_arr_1000G_5C+0.03*xRange, Pmax_arr_1000G_5C+0.08*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_arr_1000G_5C, Pmax_arr_1000G_5C), ...
    'Color','b','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','b','Margin',5,'HandleVisibility','off');
% 200 W/m², 45°C (magenta)
yline(Pmax_arr_200G_45C, '--m', 'LineWidth',1.5, 'HandleVisibility','off');
text(Vmpp_arr_200G_45C+0.03*xRange, Pmax_arr_200G_45C+0.09*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_arr_200G_45C, Pmax_arr_200G_45C), ...
    'Color','m','BackgroundColor','white','FontSize',12,'FontWeight','bold',...
    'EdgeColor','m','Margin',5,'HandleVisibility','off');
set(gca,'FontSize',12);
hold off;
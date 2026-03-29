clear all
close all 
clc

Vmpp = 37.6;
Voc = 45.3;
Impp = 5.32;
Isc = 5.84;
ki=0.000373;
kv=-0.00281;
Tstc=25+273.15;
T_45=45+273.15;         %Change this to another temperature
Gstc=1000;
G_600=600;          %Change this to another irradiance
G_200=200;          %Change this to another irradiance
np=2;           % number of parallel strings
ns=16;          % number of modules in series per string

[Iph,Vt,Rs,Io]=FourParaSTC(Impp,Vmpp,Isc,Voc);

function [Iph,Vt,Rs,Io]=FourParaSTC(Impp,Vmpp,Isc,Voc)
% This function determines the 4 paramters for the
% four parameter model from Impp, Vmpp, Isc, Voc

Iph = Isc; 
Vt = ((2*Vmpp-Voc)*(Isc - Impp))/(Impp + (Isc -Impp)*log((Isc - Impp))/Isc);
Rs = (Vt*log((Isc - Impp)/Isc) + Voc - Vmpp)/Impp;
Io = Isc/(exp(Voc/Vt));

end

%Invoke Corrected parameter function for 600W/m2 plot, T remains STC
[Iscgt_600G,Vocgt_600G,Voct_600G,Vtt_600G,Iot_600G]=FourParaCorr(Isc,Voc,ki,kv,Tstc,Tstc,Vt,G_600,Gstc);

function [Iscgt,Vocgt,Voct,Vtt,Iot]=FourParaCorr(Isc,Voc,ki,kv,T,Tstc,Vt,G,Gstc)
% This function returns the corrected parameters based on
% custom temperature and irradiance from FourParaSTC

%Isc,Voc temperature dependency
Isct=Isc+ki*(T-Tstc);
Voct=Voc+kv*(T-Tstc);

%I0, Vt temperature dependency
Vtt=Vt*(T/Tstc);
Iot=Isct/(exp(Voct/Vtt));

%Isc, Voc irradiance AND temperature dependency
Iscgt=Isct*(G/Gstc);
Vocgt=log(Iscgt/Iot)*Vtt;
end

%Initialize 1st plot for PV and IV characteristics under STC conditions
points = 1000;
V_stc = linspace(0,Voc,points);
I_stc = zeros(1,points);


%PV and IV characteristics under STC conditions
I_stc = (Isc - Io*(exp((V_stc + I_stc.*Rs)./Vt)));

%I_array = I_stc*np;


% Iarray=I_600G*np;
% Vstring=V*ns;


%Initialize 2nd plot for PV and IV characteristics under 600W/m2
%irradiance, T remains STC
V_600G = linspace(0,Vocgt_600G,points);
I_600G = zeros(1,points);

%PV and IV characteristics under 600W/m2 Irradiance and Tstc
I_600G=Iscgt_600G-Iot_600G*exp((V_600G+I_600G.*Rs)./Vtt_600G);


%Invoke Corrected parameter function for 200W/m2 plot, T remains STC
[Iscgt_200G,Vocgt_200G,Voct_200G,Vtt_200G,Iot_200G]=FourParaCorr(Isc,Voc,ki,kv,Tstc,Tstc,Vt,G_200,Gstc);

%Initialize 3rd plot for PV and IV characteristics under 200W/m2
%irradiance, T remains STC
V_200G = linspace(0,Vocgt_200G,points);
I_200G = zeros(1,points);

%PV and IV characteristics under 200W/m2 Irradiance and Tstc
I_200G=Iscgt_200G-Iot_200G*exp((V_200G+I_200G.*Rs)./Vtt_200G);




% After I_stc and V_stc
P_stc = V_stc .* I_stc;
[Pmax_stc, idx_stc] = max(P_stc);
Vmpp_stc = V_stc(idx_stc);
Impp_stc = I_stc(idx_stc);

% After I_600G and V_600G
P_600G = V_600G .* I_600G;
[Pmax_600G, idx_600G] = max(P_600G);
Vmpp_600G = V_600G(idx_600G);
Impp_600G = I_600G(idx_600G);

% After I_200G and V_200G
P_200G = V_200G .* I_200G;
[Pmax_200G, idx_200G] = max(P_200G);
Vmpp_200G = V_200G(idx_200G);
Impp_200G = I_200G(idx_200G);




%figure('Name','I-V Characteristics','NumberTitle','off')
figure('Name','I-V Characteristics for Ameresco Solar 200W (24V) MODULE 200J-V','NumberTitle','off','Color', 'white', 'Position', [100, 100, 1000, 750])
hold on
box on
grid on

% STC curve
plot(V_stc,I_stc,"Linewidth",1.5,"Color",'r', ...
    'DisplayName',sprintf('STC (1000 W/m², 25°C)'))


% 2nd plot for PV and IV characteristics under 600W/m2, T remains STC
plot(V_600G,I_600G,"Linewidth",1.5,"Color",'b', ...
    'DisplayName',sprintf('%.0f W/m², %.0f°C',G_600,Tstc-273.15))

% 3rd plot for PV and IV characteristics under 200W/m2, T remains STC
plot(V_200G,I_200G,"Linewidth",1.5,"Color",'m', ...
    'DisplayName',sprintf('%.0f W/m², %.0f°C',G_200,Tstc-273.15))

xlabel("Voltage (V)","FontSize",16,'FontWeight','bold')
ylabel("Current (A)", 'Rotation', 0, 'FontSize', 16, 'FontWeight', 'bold')
legend('Location','northeast','FontSize',14,'FontWeight','bold','FontName','Arial')
title('I-V Characteristics for Ameresco Solar 200W (24V) MODULE 200J-V','FontSize',14,'FontWeight','bold')


% Get current axis limits for offset calculations
xLim = xlim;
yLim = ylim;
xRange = diff(xLim);
yRange = diff(yLim);

% STC (red)
%xline(Vmpp_stc, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(Impp_stc, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(Vmpp_stc + 0.02*xRange, Impp_stc + 0.02*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.2f A', Vmpp_stc, Impp_stc), ...
    'Color', 'r', 'BackgroundColor', 'white', 'FontSize', 16, ...
    'EdgeColor', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% 600 W/m² (blue)
%xline(Vmpp_600G, '--b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(Impp_600G, '--b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(Vmpp_600G + 0.02*xRange, Impp_600G + 0.02*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.2f A', Vmpp_600G, Impp_600G), ...
    'Color', 'b', 'BackgroundColor', 'white', 'FontSize', 16, ...
    'EdgeColor', 'b', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% 200 W/m² (magenta)
%xline(Vmpp_200G, '--m', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(Impp_200G, '--m', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(Vmpp_200G + 0.02*xRange, Impp_200G + 0.02*yRange, ...
    sprintf('V_{mpp}=%.1f V\nI_{mpp}=%.2f A', Vmpp_200G, Impp_200G), ...
    'Color', 'm', 'BackgroundColor', 'white', 'FontSize', 16, ...
    'EdgeColor', 'm', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Set axis limits
%xlim([0, 0.75]);  % Extend x-axis to show forecast to 0.7
ylim([0, 7]);   % Adjust y-axis to accommodate forecast
% Make axis tick numbers larger
set(gca, 'FontSize', 16)

hold off


figure('Name','P-V Characteristics for Ameresco Solar 200W (24V) MODULE 200J-V','NumberTitle','off','Color', 'white', 'Position', [100, 100, 1000, 750])
hold on
box on
grid on

% STC curve
plot(V_stc,V_stc.*I_stc,"Linewidth",1.5,"Color",'r', ...
    'DisplayName',sprintf('STC (1000 W/m², 25°C)'))

% 2nd plot for PV and IV characteristics under 600W/m2, T remains STC
plot(V_600G,V_600G.*I_600G,"Linewidth",1.5,"Color",'b', ...
    'DisplayName',sprintf('%.0f W/m², %.0f°C',G_600,Tstc-273.15))

% 3rd plot for PV and IV characteristics under 200W/m2, T remains STC
plot(V_200G,V_200G.*I_200G,"Linewidth",1.5,"Color",'m', ...
    'DisplayName',sprintf('%.0f W/m², %.0f°C',G_200,Tstc-273.15))

xlabel("Voltage (V)","FontSize",16,'FontWeight','bold')
ylabel("Power (W)", 'Rotation', 0, 'FontSize', 16, 'FontWeight', 'bold')
legend('Location','northeast','FontSize',14,'FontWeight','bold','FontName','Arial')
title('P-V Characteristics for Ameresco Solar 200W (24V) MODULE 200J-V','FontSize',14,'FontWeight','bold')

% --- In the P-V figure, after plotting the curves, add annotations ---

% Get current axis limits for offset calculations
xLim = xlim;
yLim = ylim;
xRange = diff(xLim);
yRange = diff(yLim);

% STC (red)
%xline(Vmpp_stc, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(Pmax_stc, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(Vmpp_stc + 0.02*xRange, Pmax_stc + 0.02*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_stc, Pmax_stc), ...
    'Color', 'r', 'BackgroundColor', 'white', 'FontSize', 12, ...
    'FontWeight', 'bold', 'EdgeColor', 'r', 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left', 'Margin', 5, 'HandleVisibility', 'off');

% 600 W/m² (green)
%xline(Vmpp_600G, '--b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(Pmax_600G, '--b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(Vmpp_600G + 0.02*xRange, Pmax_600G + 0.02*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_600G, Pmax_600G), ...
    'Color', 'b', 'BackgroundColor', 'white', 'FontSize', 12, ...
    'FontWeight', 'bold', 'EdgeColor', 'b', 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left', 'Margin', 5, 'HandleVisibility', 'off');

% 200 W/m² (magenta)
%xline(Vmpp_200G, '--m', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(Pmax_200G, '--m', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(Vmpp_200G + 0.02*xRange, Pmax_200G + 0.02*yRange, ...
    sprintf('V_{mpp}=%.1f V\nP_{mpp}=%.1f W', Vmpp_200G, Pmax_200G), ...
    'Color', 'm', 'BackgroundColor', 'white', 'FontSize', 12, ...
    'FontWeight', 'bold', 'EdgeColor', 'm', 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left', 'Margin', 5, 'HandleVisibility', 'off');

% Set axis limits
%xlim([0, 0.75]);  % Extend x-axis to show forecast to 0.7
ylim([0, 250]);   % Adjust y-axis to accommodate forecast
% Make axis tick numbers larger
set(gca, 'FontSize', 16)


hold off

%% Design Space (GA trainer)
clc; clear; close all

% Data
CLTO = 0.5;     % Lift (Take-off)
CDTO = 0.04;    % Drag (Take-off)
CDmin = 0.025;  % Drag
g = 32.174;     % Gravity (ft^2/s)
mu = 0.04;      % Ground Friction
eta = 0.8;      % Propeller efficiency
AR = 9;         % Wing Aspect Ratio
S = 200;        % Wing Area (ft^2)
Vcr = 150;      % Cruise Speed (KTAS)
Vcl = 25;       % Climb Speed (ft/s)
Vsc = 1.666;    % Service Ceiling Speed (ft/s)
V = 80;         % Airspeed at Sea-Level (KCAS)
VLOF = 65;      % Lift-off Speed (KCAS)
SG = 900;       % Ground Run (ft)
W = 2000;       % Gross Weight (lbf)
Hcr = 8000;     % Altittude @ Cruise (ft)
Hce = 20000;    % Altittude @ Ceiling (ft)
e = 1.78*(1-0.045*AR^0.68) - 0.64;  % Ostwald's Span Efficiency
k = 1/(pi*AR*e);                    % Lift Induced Drag Constant
rho1 = 0.002378*(1-0.0000068756*Hcr)^4.2561;    % Density @ 8000 ft (slugs/ft^3)
rho2 = 0.002378;                                % Density @ Sea-Level (slugs/ft^3)
rho3 = 0.002378*(1-0.0000068756*Hce)^4.2561;    % Density @ 20000 ft (slugs/ft^3)
% Dynamic Pressure LVCT (Level Constant Velocity Turn, lbf/ft^2)
q1 = 0.5*rho1*(Vcr*1.688)^2;
% Dynamic Pressure DRC (Desired Rate of Climb, lbf/ft^2)
q2 = 0.5*rho2*(V*1.688)^2; 
% Dynamic Pressure TOD (Desired Take-off Distance, lbf/ft^2)
q3 = 0.5*rho2*(VLOF*1.688/sqrt(2))^2;
% Dynamic Pressure DCAS (Desired Cruise Air Speed, lbf/ft^2)
q4 = 0.5*rho1*(Vcr*1.688)^2;
n = 2;          % Load Factor
WS = linspace(5,60,1000);   % Weight to Area Ratio (lbf/ft^2)
%% Equations
close all

% Thrust to Weight Ratio
TW1 = q1*((CDmin./WS)+k*WS*(n/q1)^2); % T/W @ Level Constant Velocity Turn (LCVT)
TW2 = Vcl/(V*1.688) + (q2./WS)*CDmin + (k/q2)*WS;    % T/W @ Desired Rate of Climb (DRC)
TW3 = ((VLOF*1.688)^2)/(2*g*SG) + (q3*CDTO)./WS + mu*(1-(q3*CLTO./WS));   % T/W @ Desired Take-off Distance (TOD)
TW4 = (q4*CDmin./WS) + (k/q4)*WS;    % T/W @ Desired Cruise Air Speed (DCAS)
TW5 = Vsc./(sqrt((2/rho3*WS)*(sqrt(k/3/CDmin)))) + 4*sqrt(k*CDmin/3);    % T/W @ Service Ceiling (SC)

% Thrust
T1 = TW1*W; % T @ LCVT (lbf)
T2 = TW2*W; % T @ DRC (lbf)
T3 = TW3*W; % T @ TOD (lbf)
T4 = TW4*W; % T @ DCAS (lbf)
T5 = TW5*W; % T @ SC (lbf)

% Power
P1 = T1*Vcr*1.688/(eta*550);  % P_BHP @ LCVT (BHP)
P2 = T2*V*1.688/(eta*550);    % P_BHP @ DRC (BHP)
P3 = T3*V*1.688/(eta*550);    % P_BHP @ TOD (BHP)
P4 = T4*Vcr*1.688/(eta*550);  % P_BHP @ DCAS (BHP)
P5 = T5*VLOF*1.688/(eta*550); % P_BHP @ SC (BHP)

% Power @ Sea Level
P1SL = P1/(1.132*(rho1/rho2) - 0.132); % P_BHP @ LCVT (BHP)
P2SL = P2/(1.132*(rho1/rho1) - 0.132); % P_BHP @ DRC (BHP)
P3SL = P3/(1.132*(rho2/rho2) - 0.132); % P_BHP @ TOD (BHP)
P4SL = P4/(1.132*(rho1/rho2) - 0.132); % P_BHP @ DCAS (BHP)
P5SL = P5/(1.132*(rho3/rho2) - 0.132); % P_BHP @ SC (BHP)

% Lift Coeficient --> Stall Velocity
CL1 = 2*WS/(rho2*(65*1.688)^2); % CL_max @ 65 KCAS
CL2 = 2*WS/(rho2*(60*1.688)^2); % CL_max @ 60 KCAS
CL3 = 2*WS/(rho2*(55*1.688)^2); % CL_max @ 55 KCAS
CL4 = 2*WS/(rho2*(50*1.688)^2); % CL_max @ 50 KCAS
CL5 = 2*WS/(rho2*(45*1.688)^2); % CL_max @ 45 KCAS
CL6 = 2*WS/(rho2*(40*1.688)^2); % CL_max @ 40 KCAS

figure;
set(gcf, 'Position', [100, 100, 700, 600]);
colororder({'k','k'})
yyaxis left
hold on
plot(WS, P1SL, 'r-', 'LineWidth', 1)
plot(WS, P2SL, 'b-', 'LineWidth', 1)
plot(WS, P3SL, 'LineWidth', 1, 'Color', '#D95319', 'LineStyle', '-')
plot(WS, P5SL, 'LineWidth', 1, 'Color', '#77AC30', 'LineStyle', '-')
plot(WS, P4SL, 'm-', 'LineWidth', 1)
hold off
ylabel('Power @ Sea-Level - $\mathrm{P}_\mathrm{BHP}$', 'Interpreter', 'latex')

yyaxis right
hold on
plot(WS, CL1, 'k--', 'LineWidth', 1)
plot(WS, CL2, 'k--', 'LineWidth', 1)
plot(WS, CL3, 'k--', 'LineWidth', 1)
plot(WS, CL4, 'k--', 'LineWidth', 1)
plot(WS, CL5, 'k--', 'LineWidth', 1)
plot(WS, CL6, 'k--', 'LineWidth', 1)
hold off
ylabel('$\mathrm{CL}_\mathrm{Max}$', 'Interpreter', 'latex')

xlim([5, 60])
ylim([0, 3])

% Add text labels for CLmax
text(WS(round(end/3)) + .5, CL1(round(end/3)), '65 KCAS', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 9, 'Rotation', 25)
text(WS(round(end/3)) + .5, CL2(round(end/3)), '60 KCAS', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 9, 'Rotation', 25)
text(WS(round(end/3)) + .5, CL3(round(end/3)), '55 KCAS', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 9, 'Rotation', 25)
text(WS(round(end/3)) + .5, CL4(round(end/3)), '50 KCAS', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 9, 'Rotation', 25)
text(WS(round(end/5)) + .5, CL5(round(end/5)), '45 KCAS', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 9, 'Rotation', 25)
text(WS(round(end/6)) + .5, CL6(round(end/6)), '40 KCAS', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 9, 'Rotation', 25)
legend('LCVT','DRC','TOD','DCAS','SC','Interpreter','latex','Location','northeast')
xlabel('Wing Loading - $W/S \left[\mathrm{lbf}/\mathrm{ft}^2 \right]$', 'Interpreter', 'latex')
title('Design Space for GA Trainer Aircraft', 'Interpreter', 'latex')
grid on

figure;
hold on
plot(WS,TW1,'r','LineWidth',1)
plot(WS,TW2,'b','LineWidth',1)
plot(WS,TW3,'LineWidth',1,'Color','#D95319')
plot(WS,TW4,'m','LineWidth',1)
plot(WS,TW5,'LineWidth',1,'Color','#77AC30')
hold off
xlim([5,60])
legend('LCVT','DRC','TOD','DCAS','SC','Interpreter','latex','Location','bestoutside')
xlabel('Wing Loading - $W/S \left[\mathrm{lbf}/\mathrm{ft}^2 \right]$','Interpreter','latex')
ylabel('Thrust Loading - $T/W$','Interpreter','latex')
title('Weight to Area Ratio Vs. Thrust to Weight Ratio','Interpreter','latex')
grid on

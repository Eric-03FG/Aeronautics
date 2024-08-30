%% Banking Constraint Diagram
clc; clear

% Data
S = 144.9;          % Wing Area (ft^2)
AR = 10.12;         % Wing Aspect Ratio
W = 3400;           % Gross Weight (lbf)
P = 310;            % Max Power @ Sea Level (HP)
CDmin = 0.02541;    % Zero Lift Coefficient
eta = 0.8;          % Propeller Efficiency
CLmax = 1.41;       % Max Lift Coefficient
nmax = 3.8;         % Max Load Factor
p = 0.002378;       % Density @ Sea Level (slugs/ft^3)
e = 1.78*(1-0.045*(AR^0.68))-0.64;  % Oswald's Span Efficiency
k = 1/(pi*AR*e);                    % Lift Induced Drag
Vs = sqrt(2*W/(p*S*CLmax));         % Stall Speed (ft/2)
Vs = Vs/1.688;                      % Stall Speed (KCAS)
ASmax = 180;                        % Max Airspeed (KCAS)
AS = linspace(10,ASmax,1000);       % Airspeed (KCAS)
T = (eta*550*P)/(ASmax);            % Thrust (lbf)
q = 0.5*p*(AS*1.688).^2;            % Dynamic Pressure (lbf/ft^2)

%% Equations
close all

% Max Banking Load Factor
n1 = (q*S/W).*sqrt((1/k)*(T./(q*S)-CDmin));
% Max Stall Load Banking
n2 = ((AS*1.688).^2)*p*CLmax*S/(2*W);
% Limit Load
n3 = nmax*ones(1,1000);
% Stall Speed (KCAS)
Vsg = Vs*ones(1,1000);
x = linspace(0,180,1000);
% Normal Operating Speed
VNO = 178*ones(1,1000);

% Encontrar el punto de intersección
idx_intersection = find(n1 <= n2, 1);

figure;
hold on

plot(AS(idx_intersection:end), n1(idx_intersection:end), 'b', 'LineWidth', 1)
plot(AS(1:idx_intersection), n2(1:idx_intersection), 'r', 'LineWidth', 1)
plot(AS(1:idx_intersection), n1(1:idx_intersection), 'b--', 'LineWidth', 1)
plot(AS(idx_intersection:end), n2(idx_intersection:end), 'r--', 'LineWidth', 1)

plot(x,n3,'k--','LineWidth',1)
plot(Vsg,x,'k--','LineWidth',1)
plot(VNO,x,'k--','LineWidth',1)

text(80, nmax+0.2, 'Limit Load, $3.8g$', 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment', 'left')
text(Vs-2, 2.5, '$1gV_S$', 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment', 'right')
text(178-2, 2.5, '$V_\mathrm{NO}$', 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment', 'right')

hold off
ylim([0 5])
xlim([0 180])
titleHandle = title('Banking Constraint Diagram', 'Interpreter', 'latex');
xlabelHandle = xlabel('Calibrated Airspeed, KAS', 'Interpreter', 'latex');
ylabelHandle = ylabel('Load Factor', 'Interpreter', 'latex');
legendHandle = legend('Max Banking Load Factor', 'Max Stall Load Banking', 'Interpreter', 'latex', 'location', 'southoutside');
grid on

set(titleHandle, 'FontSize', 12)
set(xlabelHandle, 'FontSize', 11)
set(ylabelHandle, 'FontSize', 11)
set(legendHandle, 'FontSize', 11)
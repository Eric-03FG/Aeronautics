%% Thrust Loading Analysis
% Se considera un ala semi-elíptica combinada con un perfil aerodinámico
% NACA4412 y un motor O.S. Max 0.46 AX II

clc; clear; close all;

% Importar datos de CL y Alpha a número de Reynolds de 200,000
filename = 'NACA4412_Re_200000.xlsx';
data = readmatrix(filename, 'Range', 'A2:B104');

% Graficar CL Vs. Alpha
figure;
plot(data(:,1), data(:,2), 'LineWidth', 1, 'Color', 'b');
grid on;
xlabel('$\alpha$ (Grados)', 'Interpreter', 'latex');
ylabel('$\mathrm{C}_\mathrm{L}$', 'Interpreter', 'latex');
title('Coeficiente de sustentaci\''on Vs. \''Angulo de ataque','Interpreter','latex')
legend('NACA4412, Re = 200,000','Interpreter','latex','Location','southeast')

% Parámetros definidos
P = 1.63 * 550;        % Potencia nominal del motor (ft*lbf/s)
Beta = 0.7;            % Porcentaje de la potencia a utilizar
g = 32.174;            % Gravedad (ft/s^2)
rho = 2.378e-3;        % Densidad del aire (slugs/ft^3)
n = 0.7;               % Eficiencia de la hélice
CL = data(46, 2);      % Coeficiente de sustentación @ Alpha = 3°
mb = 2.5/14.594;       % Masa de los componentes de la aeronave (slugs) --> 2.5kg

% Parámetros a modificar
% Definir rango de posibilidades
r1 = linspace(0.25, 0.5, 1000);     % Radio menor (ft)
r2 = linspace(1, 3, 1000);          % Radio mayor (ft)
[R1, R2] = meshgrid(r1, r2);

% Funciones a calcular
S = (5/8) * pi * R2 .* R1;              % Área de la elipse asimétrica (ft^2)
me = (4/14.594)*(S./max(S(:)));         % El mayor área del ala generará 4kg extra
mt = mb + me;                           % Masa total de la aeronave (slugs)
W = mt * g;                             % Peso de la aeronave (lbf)
Vs = sqrt((2 * W) ./ (CL * rho * S));   % Velocidad de stall (ft/s)
T = (n * P * Beta) ./ Vs;               % Empuje (lbf)
TW = T ./ W;                            % Thrust Loading (T/W)

% Encontrar el valor máximo de TW y los índices correspondientes
[max_TW, idx] = max(TW(:));

% Convertir el índice lineal a índices de la matriz
[row, col] = ind2sub(size(TW), idx);

% Obtener los valores de r1 y r2 correspondientes
r1_max = r1(col);
r2_max = r2(row);

% Mostrar los resultados
fprintf('El valor máximo de T/W es: %.4f\n', max_TW);
fprintf('El radio menor es: %.4f ft\n', r1_max);
fprintf('El radio mayor es: %.4f ft\n', r2_max);
fprintf('El área es: %.4f ft^2\n', S(idx));
fprintf('La cuerda raíz es: %.4f ft\n', r1_max*2);
fprintf('La envergadura es: %.4f ft\n', r2_max*2);

% Graficar las posibilidades
figure;
colormap hot;
surf(r1, r2, TW);
shading interp;
title('Thrust Loading - $T/W$', 'Interpreter', 'latex');
xlabel('$r_1$ (ft)', 'Interpreter', 'latex');
ylabel('$r_2$ (ft)', 'Interpreter', 'latex');
zlabel('$T/W$', 'Interpreter', 'latex');
colorbar;
hold on;
plot3(r1_max, r2_max, max_TW, '*b', 'MarkerSize', 10);
hold off;
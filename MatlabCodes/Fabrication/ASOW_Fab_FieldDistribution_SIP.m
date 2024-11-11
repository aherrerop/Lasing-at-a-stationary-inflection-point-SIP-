%% Field Amplitudes at every cell for a finite length structure
%
clear all; close all; clc
% Parameters SIP:
Freq_Fundamental = 386.6; % In THz
RealRefrIndex_Fundamental = 2.2585;
Freq_DownConverted = 193.3; % In THz
RealRefrIndex_DownConverted = 2.0508;
Radius = 10e-6; %m
Freq_Center = Freq_Fundamental; % Fundamental at 776 nm. In THz
RealRefrIndex = RealRefrIndex_Fundamental; % LiNbO3 at 776 nm
CouplLength = 0.490842584539929;
Alpha = 1.06595385416124;
Alpha_2 = 1.06225804311915;


Loss = 0;%4.7e-6; % Loss is positive!        
Gain = 0; % Gain is negative!          
ImagRefrIndex = Gain + Loss; % Imaginary part of the Effective refractive index. Negative if gain.
% Constants
EffRefrIndex = RealRefrIndex -1j*ImagRefrIndex;
SpeedLight = 3e8; % m/s

NumberOfUnitCells = 20; % Number of Unit Cells Used

% Let's find the values at each side:
PortValues = zeros(12,1);
[PortValues, TransferMatrix, T_aux] = ASOW_Fab_PortValues_OutputMatrices (NumberOfUnitCells, Freq_Center, CouplLength, Radius, Alpha, Alpha_2, EffRefrIndex);

E_0 = PortValues(1:6);
E_End = PortValues(7:12);

% Now let's find the values in the middle
% Forwards
E_f = zeros(6,NumberOfUnitCells+1);
E_f(:,1) = E_0;
% E_f(NumberOfUnitCells,:) = E_End;
for ii = 2:NumberOfUnitCells
    E_f(:,ii)= TransferMatrix*E_f(:,ii-1);
end
E_f(:,NumberOfUnitCells+1) = T_aux*E_f(:,ii);

% Backwards
E_b = zeros(6,NumberOfUnitCells+1);
E_b(:,end) = E_End;
E_b(:,end-1) = T_aux^-1*E_b(:,end);
for ii = 2:NumberOfUnitCells
    E_b(:,end-ii) = TransferMatrix^-1*E_b(:,end-ii+1);
end

UnitCell_vec = linspace(0,NumberOfUnitCells, NumberOfUnitCells+1); 
% Let's plot
w_line = 2;

figure(1)
hold on
plot(UnitCell_vec,abs(E_f(1,:)),'ok','linewidth',w_line)
plot(UnitCell_vec,abs(E_b(1,:)),'-b','linewidth',1)
hold off
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|E_1^+(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 15])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
%
figure(2)
hold on
plot(UnitCell_vec,abs(E_f(2,:)),'ok','linewidth',w_line)
plot(UnitCell_vec,abs(E_b(2,:)),'-b','linewidth',1)
hold off
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|E_1^-(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 15])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
%%
figure(3)
hold on
plot(UnitCell_vec, abs(E_f(2,:) + E_f(1,:)),'ok','linewidth',w_line)
plot(UnitCell_vec, abs(E_b(2,:) + E_b(1,:)),'-b','linewidth',1)
hold off
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|E_1(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 15])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');

%% Nonlinear polarization

% P = 2\varepsilon_0 d_eff E,
Abs_permittivity = 8.85e-12; % In F/m
d_33 = -25.2; % In pmV^-1
d_eff = (4/(pi^2))*d_33; % In pmV^-1

P_f(1,:) = 2*Abs_permittivity*d_eff.*E_f(1,:);
P_f(2,:) = 2*Abs_permittivity*d_eff.*E_f(2,:);
% Backward-calculated (from end to beginning) should be the same as
% forward-calculated (from beginning to end)
P_b(1,:) = 2*Abs_permittivity*d_eff.*E_b(1,:);
P_b(2,:) = 2*Abs_permittivity*d_eff.*E_b(2,:);

% Let's plot the polarization
w_line = 2;

figure(1)
hold on
plot(UnitCell_vec,abs(P_f(1,:)),'ok','linewidth',w_line)
plot(UnitCell_vec,abs(P_b(1,:)),'-b','linewidth',1)
hold off
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|P_1^+(n)|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 15])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
%
figure(2)
hold on
plot(UnitCell_vec,abs(P_f(2,:)),'ok','linewidth',w_line)
plot(UnitCell_vec,abs(P_b(2,:)),'-b','linewidth',1)
hold off
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|P_1^-(n)|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 15])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
%
figure(3)
hold on
plot(UnitCell_vec, abs(P_f(2,:) + P_f(1,:)),'ok','linewidth',w_line)
plot(UnitCell_vec, abs(P_b(2,:) + P_b(1,:)),'-b','linewidth',1)
hold off
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|P_1(n)|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 15])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');

%% Green's Function

% Now we have the quadratic polarization. Now we need the field amplitude
% at the downconverted frequency.

Green_Function = exp(i*Wavenumber*abs(z_t-z_s))/(4*pi*abs(z_t-z_s)^3);

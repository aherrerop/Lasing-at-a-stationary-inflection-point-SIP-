%% Calculate Transfer Function with Gain at a SIP
% Normalized Group Delay & Q.
clear all; close all; clc
% Parameters SIP, RBE R = 6:
Freq_Center = 193.54; % In THz
RealRefrIndex = 2.38446; % 450x220nm^2. Average between n_e = 2.427 & n_o = 2.34 (even & odd modes from Tarek's paper) 
CouplLength = 0.499443761395272;
Alpha = 1.02710447969426;
Alpha_2 = 1.01295874936608;
Radius = 6e-6; % In m

% Constants
Bending_Losses = 2.29e-5; % Imag refr index for bending losses
Waveguide_Prop_Losses = 6.88e-6; % Imag refr index for waveguide prop. losses
ImagRefrIndex = 0;% Bending_Losses + Waveguide_Prop_Losses;
EffRefrIndex_vec = RealRefrIndex - 1j*ImagRefrIndex;
SpeedLight = 3e8; % m/s

Freq_RBE_Top = 193.636517303461;

Freq_RBE_Bottom = 193.470964192839;

NumberOfUnitCells = [6 8 10]; % Number of Unit Cells Used
% Frequency Sweep
Freqmax = Freq_Center + 0.01;
Freqmin = Freq_Center - 0.01;
Freqsteps = 5000;
% Frequency Vector
Freq = linspace(Freqmin, Freqmax, Freqsteps);
Freq = Freq';
% Length of the finite length structure

% Vector with the field amplitudes at either side for all frequencies &
% lengths of the finite-length structure.
PortValues = zeros(12,Freqsteps,length(NumberOfUnitCells));
for ii=1:length(NumberOfUnitCells)
    for jj = 1: Freqsteps
        % This function generates the 12x12 system matrix with the finite
        % length structure transfer matrix and the boundary condition.
        [PortValues(:,jj,ii), TransferMatrices_Vec(:,:,jj)] = ASOW_Fab_PortValues (NumberOfUnitCells(ii), Freq(jj), CouplLength, Radius, Alpha, Alpha_2, EffRefrIndex_vec);
        det(TransferMatrices_Vec(:,:,jj));
    end
end

% Let's assign the righ units for the Freq_Center.
Freq_Center = 1e12*Freq_Center;
Freq = 1e12*linspace(Freqmin, Freqmax, Freqsteps);
Freq = Freq';

% Transfer Function & Q
TransferFunction = zeros(Freqsteps,length(NumberOfUnitCells));
PhaseTransferFunction = zeros(Freqsteps,length(NumberOfUnitCells));
% The group delay has one less element bc 
GroupDelay = zeros(Freqsteps-1,length(NumberOfUnitCells));
Q = zeros(Freqsteps-1,length(NumberOfUnitCells));
% Define the vector tau_0 (group delay going through a straight waveguide
% of the same length as each finite length structure).
GroupDelayBaseline = EffRefrIndex_vec*(NumberOfUnitCells)*(2*pi*Radius + 2*(Alpha + Alpha_2)*Radius) / (SpeedLight);

for ii = 1:length(NumberOfUnitCells)
    TransferFunction(:,ii) = PortValues(7,:,ii); % The incident field amplitude is the unity.
    TransferFunctiondB(:,ii) = 20*log10(abs(TransferFunction(:,ii))); % Take the magnitude in dB
    PhaseTransferFunction(:,ii) = angle(TransferFunction(:,ii)); % Take the phase of the transfer function
    % The group delay is the derivate of the phase of the transfer function
    % wrt the angular frequency (hence the 2pi dividing). We will add the
    % sign later.
    GroupDelay(:,ii) = diff(PhaseTransferFunction(:,ii))./((2*pi).*diff(Freq));
    [GroupDelay_min, index_min] = min(GroupDelay(:,ii)); % Because the group delay is negative now. The sign will be added later.
    GroupDelay_norm(:,ii) = GroupDelay(:,ii) / GroupDelayBaseline(ii);  % Normalize the group delay
    Freq_res(ii) = Freq(index_min); % We get the frequency of the highest peak of the group delay (lowest bc we have yet to change the sign)
    % We can use different definitions of Q. We use the unnormalized one.
    Q(:,ii) = pi*Freq_res(ii).*GroupDelay_norm(:,ii); % Q with normalized group delay
    Q_unnormalized(:,ii) = pi*Freq_res(ii).*GroupDelay(:,ii); % Q unnormalized
    ReflectionFunction(:,ii) = PortValues(2,:,ii); % The incident field amplitude is the unity.
    ReflectionFunctiondB(:,ii) = 20*log10(abs(ReflectionFunction(:,ii))); % Take the magnitude in dB
    PhaseReflectionFunction(:,ii) = angle(ReflectionFunction(:,ii)); % Take the phase of the Reflection function
end

Freq_plt = Freq./Freq_Center; % Frequency normalized to the SIP frequency
% Lossless check:
% Want to check that |S_11|^2 + |S_21|^2 = 1
ShouldBeOne = abs(TransferFunction).^2 + abs(ReflectionFunction).^2;
% %
% w_line = 2;
% figure(5)
% hold on
% plot(Freq_plt, ShouldBeOne,'-','linewidth',w_line)
% xline(1,'b','Linewidth',1.5)
% hold off
% grid on
% ylabel('$|T_f|^2 + |R_f|^2$','FontSize', 20,'Interpreter','latex');
% xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
% ylim([0 2])

%
%
close all;
figure(1)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, TransferFunctiondB(:,ii),'-','linewidth',w_line);
end
xline(1,'--b','Linewidth',1);
xline(Freq_RBE_Top/Freq_Center*1e12,'--r','Linewidth',1);
xline(Freq_RBE_Bottom/Freq_Center*1e12,'--r','Linewidth',1);
legend('$N$=6', '$N$=8','$N$=10','FontSize', 20,'Interpreter','latex');
hold off
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
ylabel('$\mid$ Transfer function $\mid$ (dB)','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
% ylim([-60 2])
xlim([Freq_plt(1) Freq_plt(end)])

%%
figure(2)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, PhaseTransferFunction(:,ii),'-','linewidth',w_line)
end
hold off
pbaspect([1.3 1 1]);
legend('$N$=6', '$N$=8','$N$=10','FontSize', 20,'Interpreter','latex');
set(gca,'FontSize',20,'FontName', 'Times New Roman');
% title('Transfer Function [dB] versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
ylabel('Angle Transfer functions','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
grid on
%%

figure(3)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, ReflectionFunctiondB(:,ii),'-','linewidth',w_line)
end
xline(1,'b','Linewidth',1.5)
xline(Freq_RBE_Top*1e12/Freq_Center,'r','Linewidth',1.5)
xline(Freq_RBE_Bottom*1e12/Freq_Center,'r','Linewidth',1.5)
hold off
grid on
legend('$N$=6', '$N$=8','$N$=10','FontSize', 20,'Interpreter','latex');
hold off
set(gca,'FontSize',20,'FontName', 'Times New Roman');
ylim([-10 1])
ylabel('$\mid$ Reflection function $\mid$ (dB)','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
%%
figure(4)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, PhaseReflectionFunction(:,ii),'-','linewidth',w_line)
end
hold off
pbaspect([1.3 1 1]);
legend('$N$=6', '$N$=8','$N$=10','FontSize', 20,'Interpreter','latex');
set(gca,'FontSize',20,'FontName', 'Times New Roman');
ylabel('Angle Reflection functions','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
grid on
%%
% Filter the discontinuities in the Angle Transfer Function,
% which causes discontinuities in the Group Delay:
index_1 = find(GroupDelay_norm(:,1)<0);
index_12 = find(GroupDelay_norm(:,1)<(-200));
index_1 = setdiff(index_1,index_12);
Freq_Q_plt_1 = Freq_plt(index_1);
GroupDelay_plt_1 = GroupDelay_norm(index_1,1);
Q_norm_plt_1 = Q_normalized(index_1,1);
Q_plt_1 = Q_unnormalized(index_1,1);
   
index_2 = find(GroupDelay_norm(:,2)<0);
index_21 = find(GroupDelay_norm(:,2)<(-200));
index_2 = setdiff(index_2,index_21);
Freq_Q_plt_2 = Freq_plt(index_2);
GroupDelay_plt_2 = GroupDelay_norm(index_2,2);
Q_norm_plt_2 = Q_normalized(index_2,2);
Q_plt_2 = Q_unnormalized(index_2,2);

index_3 = find(GroupDelay_norm(:,3)<0);
index_32 = find(GroupDelay_norm(:,3)<(-200));
index_3 = setdiff(index_3,index_32);
Freq_Q_plt_3 = Freq_plt(index_3);
GroupDelay_plt_3 = GroupDelay_norm(index_3,3);
Q_norm_plt_3 = Q_normalized(index_3,3);
Q_plt_3 = Q_unnormalized(index_3,3);

index_4 = find(GroupDelay(:,4)<0);
Freq_Q_plt_4 = Freq_plt(index_4);
GroupDelay_plt_4 = GroupDelay_norm(index_4,4);
Q_norm_plt_4 = Q_normalized(index_4,4);

index_5 = find(GroupDelay(:,5)<0);
Freq_Q_plt_5 = Freq_plt(index_5);
GroupDelay_plt_5 = GroupDelay_norm(index_5,5);
Q_norm_plt_5 = Q_normalized(index_5,5);

index_6 = find(GroupDelay(:,6)<0);
Freq_Q_plt_6 = Freq_plt(index_6);
GroupDelay_plt_6 = GroupDelay_norm(index_6,6);
Q_norm_plt_6 = Q_normalized(index_6,6);

% Q_1 = pi*Freq_res(1).*abs(GroupDelay_plt_1);
% Q_2 = pi*Freq_res(2).*abs(GroupDelay_plt_2);
% Q_3 = pi*Freq_res(3).*abs(GroupDelay_plt_3);
% Q_4 = pi*Freq_res(4).*abs(GroupDelay_plt_4);
% Q_5 = pi*Freq_res(5).*abs(GroupDelay_plt_5);
% Q_6 = pi*Freq_res(6).*abs(GroupDelay_plt_6);




%%
w_line = 2;
figure(5)
hold on
plot(Freq_Q_plt_1,-GroupDelay_plt_1,'-','linewidth',w_line)
plot(Freq_Q_plt_2,-GroupDelay_plt_2,'-','linewidth',w_line)
plot(Freq_Q_plt_3,-GroupDelay_plt_3,'-','linewidth',w_line)
% plot(Freq_Q_plt_4,-GroupDelay_plt_4,'-','linewidth',w_line)
% plot(Freq_Q_plt_5,-GroupDelay_plt_5,'-','linewidth',w_line)
% plot(Freq_Q_plt_6,-GroupDelay_plt_6,'-','linewidth',w_line)
hold off
set(gca,'FontSize',20,'FontName', 'Times New Roman');
legend('$N$=6', '$N$=8','$N$=10','FontSize', 20,'Interpreter','latex');
grid on
ylabel('$\tau_g / \tau_{0}$','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
xlim([1-2e-4 1+2e-4])
%%
figure(6)
hold on
plot(Freq_Q_plt_1,-Q_plt_1,'-','linewidth',w_line)
plot(Freq_Q_plt_2,-Q_plt_2,'-','linewidth',w_line)
plot(Freq_Q_plt_3,-Q_plt_3,'-','linewidth',w_line)
% plot(Freq_Q_plt_4,Q_4,'-','linewidth',w_line)
% plot(Freq_Q_plt_5,Q_5,'-','linewidth',w_line)
% plot(Freq_Q_plt_6,Q_6,'-','linewidth',w_line)
hold off
set(gca,'FontSize',20,'FontName', 'Times New Roman');
legend('$N$=6', '$N$=8','$N$=10','FontSize', 20,'Interpreter','latex');
grid on
ylabel('$Q$','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');


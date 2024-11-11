%% Plot dispersion diagram ASOW-SIP
%% 1

clear all; close all

% Parameters SIP:
Freq_Center = 193.54;
CouplCoeff_k1 = 0.49832327234602;
Radius = 10e-6;
Alpha = 1.15240339832846;
Alpha_2 = 0.980630321583591;

% Constants
Loss = 0; % Loss is positive!        
Gain = 0; % Gain is negative!          
ImagRefrIndex = Gain + Loss; % Imaginary part of the Effective refractive index. Negative if gain.
% Constants
EffRefrIndex = 2.362 -1j*ImagRefrIndex;
SpeedLight = 3e8; % m/s


% Frequency Sweep
Freqmax = Freq_Center + 0.040;
Freqmin = Freq_Center - 0.040;
Freqsteps = 5001;

Freq_low = linspace(Freqmin, Freq_Center, Freqsteps);
Freq_high = linspace(Freq_Center,Freqmax, Freqsteps);

Freq_low_plt = Freq_low./Freq_Center;
Freq_high_plt = Freq_high./Freq_Center;
TransferMatrices_Vec = zeros(6,6,Freqsteps);
Eigenvectors = zeros(6,Freqsteps);
for ii = 1:Freqsteps    
    [TransferMatrices_Vec(:,:,ii)] =  ASOW_Gain_TransferMatrix (Freq_low(ii), CouplCoeff_k1, Radius, Alpha, Alpha_2, EffRefrIndex);
end
[EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
for ii = 1:Freqsteps
    Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
    Sigma(ii) = ASOW_Gain_SIP_Finding_Hyperdistance(TransferMatrices_Vec(:,:,ii));
end

for ii = 1:6
    kdlow1 = Eigenvalues(1,:);
    kdlow2 = Eigenvalues(2,:);
    kdlow3 = Eigenvalues(3,:);
    kdlow4 = Eigenvalues(4,:);
    kdlow5 = Eigenvalues(5,:);
    kdlow6 = Eigenvalues(6,:);
end
%% 3

%clear all;

% Parameters SIP:
Freq_Center = 193.54;
CouplCoeff_k1 = 0.49832327234602;
Radius = 10e-6;
Alpha = 1.15240339832846;
Alpha_2 = 0.980630321583591;

% Frequency Sweep
Freqmax = Freq_Center + 0.040;
Freqmin = Freq_Center - 0.040;
Freqsteps = 5001;

Freq_low = linspace(Freqmin, Freq_Center, Freqsteps);
Freq_high = linspace(Freq_Center,Freqmax, Freqsteps);

Freq_low_plt = Freq_low./Freq_Center;
Freq_high_plt = Freq_high./Freq_Center;
TransferMatrices_Vec = zeros(6,6,Freqsteps);
Eigenvectors = zeros(6,Freqsteps);
for ii = 1:Freqsteps    
    [TransferMatrices_Vec(:,:,ii)] =  ASOW_Gain_TransferMatrix (Freq_high(ii), CouplCoeff_k1, Radius, Alpha, Alpha_2, EffRefrIndex);
end
[EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
for ii = 1:Freqsteps
    Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
    Sigma(ii) = ASOW_Gain_SIP_Finding_Hyperdistance(TransferMatrices_Vec(:,:,ii));
end

for ii = 1:6
    kdhigh1 = Eigenvalues(1,:);
    kdhigh2 = Eigenvalues(2,:);
    kdhigh3 = Eigenvalues(3,:);
    kdhigh4 = Eigenvalues(4,:);
    kdhigh5 = Eigenvalues(5,:);
    kdhigh6 = Eigenvalues(6,:);
end
%% 4
w_line=3;
% close all
BC=[255   , 135    ,0]/255;
GC='g';%[153   , 255    ,0]/255;

figure(1);
hold on;
plot(real(kdhigh3)/pi,Freq_high_plt,'-','linewidth',w_line,'color','g');
plot(real(kdhigh4)/pi,Freq_high_plt,'-','linewidth',w_line,'color','r');
plot(real(kdhigh5)/pi,Freq_high_plt,'--','linewidth',w_line,'color',BC);
plot(real(kdhigh6)/pi,Freq_high_plt,'--','linewidth',w_line,'color','b');


plot(real(kdlow1)/pi,Freq_low_plt,'-','linewidth',w_line,'color','b');
plot(real(kdlow2)/pi,Freq_low_plt,'-','linewidth',w_line,'color',BC);
plot(real(kdlow3)/pi,Freq_low_plt,'--','linewidth',w_line,'color','g');
plot(real(kdlow4)/pi,Freq_low_plt,'--','linewidth',w_line,'color','r');
plot(real(kdlow5)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
plot(real(kdlow6)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
plot(real(kdhigh1)/pi,Freq_high_plt,'-','linewidth',w_line,'color','k');
plot(real(kdhigh2)/pi,Freq_high_plt,'-','linewidth',w_line,'color','k');

set(gca,'FontSize',24,'FontName', 'Times New Roman');
% set(gca,'XTick',[], 'YTick', []);
grid on
xlabel('Re($kd / \pi$)','FontSize', 36,'Interpreter','latex');
ylabel('$\omega / \omega_s$','FontSize', 36,'Interpreter','latex');
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
% pbaspect([1.3 1 1]);
%

% figure(2)
% hold on;
% plot(imag(kdhigh3)/pi,Freq_high_plt,'-','linewidth',w_line,'color','b');
% plot(imag(kdhigh5)/pi,Freq_high_plt,'-','linewidth',w_line,'color','r');
% plot(imag(kdhigh4)/pi,Freq_high_plt,'--','linewidth',w_line,'color','g');
% plot(imag(kdhigh2)/pi,Freq_high_plt,'--','linewidth',w_line,'color',BC);
% 
% 
% plot(imag(kdlow4)/pi,Freq_low_plt,'-','linewidth',w_line,'color',BC);
% plot(imag(kdlow1)/pi,Freq_low_plt,'-','linewidth',w_line,'color','g');
% plot(imag(kdlow2)/pi,Freq_low_plt,'--','linewidth',w_line,'color','r');
% plot(imag(kdlow3)/pi,Freq_low_plt,'--','linewidth',w_line,'color','b');
% plot(imag(kdlow5)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
% plot(imag(kdlow6)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
% plot(imag(kdhigh6)/pi,Freq_high_plt,'-','linewidth',w_line,'color','k');
% plot(imag(kdhigh1)/pi,Freq_high_plt,'-','linewidth',w_line,'color','k');
% 
% xlabel('Im($kd / \pi$)','FontSize', 24,'Interpreter','latex');
% ylabel('$\omega / \omega_s$','FontSize', 24,'Interpreter','latex');
% axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
% grid on;
% % pbaspect([1.3 1 1]);
% 
% 
% % % 
% figure(3);
% hold on;
% plot(real(kdhigh3)/pi,imag(kdhigh3)/pi,'-','linewidth',w_line,'color',BC);
% plot(real(kdhigh4)/pi,imag(kdhigh4)/pi,'-','linewidth',w_line,'color','b');
% plot(real(kdhigh6)/pi,imag(kdhigh6)/pi,'-','linewidth',w_line,'color','r');
% plot(real(kdhigh5)/pi,imag(kdhigh5)/pi,'-','linewidth',w_line,'color','g');
% 
% plot(real(kdlow3)/pi,imag(kdlow3)/pi,'-','linewidth',w_line,'color',BC);
% plot(real(kdlow1)/pi,imag(kdlow1)/pi,'-','linewidth',w_line,'color','g');
% plot(real(kdlow2)/pi,imag(kdlow2)/pi,'-','linewidth',w_line,'color','r');
% plot(real(kdlow4)/pi,imag(kdlow4)/pi,'-','linewidth',w_line,'color','b');
% 
% 
% plot(real(kdhigh2)/pi,imag(kdhigh2)/pi,'-','linewidth',w_line,'color','k');
% plot(real(kdhigh1)/pi,imag(kdhigh1)/pi,'-','linewidth',w_line,'color','k');
% plot(real(kdlow5)/pi,imag(kdlow5)/pi,'-','linewidth',w_line,'color','k');
% plot(real(kdlow6)/pi,imag(kdlow6)/pi,'-','linewidth',w_line,'color','k');
% 
% grid on;
% set(gca,'FontSize',24,'FontName', 'Times New Roman');
% xlabel('Re($kd / \pi$)','FontSize', 24,'Interpreter','latex');
% ylabel('Im($kd / \pi$)','FontSize', 24,'Interpreter','latex');
% grid on;
% % pbaspect([1.3 1 1]);
% axis([-0.8,0.8, -0.4,0.4]);
% 
% figure(4);
% hold on;
% % plot(Freq_low_plt,Sigma,'b.');
% plot(Freq_high_plt,Sigma,'b.');
% hold off;
% ylabel('Coalescence parameter, $\sigma$','FontSize', 24,'Interpreter','latex');
% xlabel('$\omega / \omega_s$','FontSize', 24,'Interpreter','latex');
% set(gca,'FontSize',24,'FontName', 'Times New Roman');
% xlim([Freqmin/Freq_Center, Freqmax/Freq_Center]);
% grid on;
% 

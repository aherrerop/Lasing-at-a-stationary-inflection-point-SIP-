%% EPD_Finding_Mock_Code
clear all; close all; clc
%{
To find the SIP we work with three files: 
- ASOW_Fab_SIP_Finding_Script
- ASOW_Fab_SIP_Finding_Function
- ASOW_Fab_SIP_Finding_Hyperdistance

This file is the EPD_Finding_Mock_Code. In it we minimize the function
EPD_Finding_Function, which depends on the EPD_Finding_Hyperdistance, where
we calculate the Hyperdistance, sigma. It is a figure of merit which
quantifies the coalescence of the eigenvectors by calculating the angles
between them.

This code initializes the parameters within a set interval and finds the
parameters that will give an EPD by minimizing sigma. In
EPD_Finding_Function, you can impose conditions on the values of the
parameters or the position of the EPD in the BZ.
%}

% Start the loop: 
Sigma_EPD = [];
Eigenvalues_EPD = [];
X_EPD = [];
Ksd_pi_EPD = [];
Eigenvectors_EPD = [];
NumberOfIterations = 1e3;

% Eigenvalues_X = zeros(6,1);
% Set a Radius, a gap, a width
Radius = 6e-6; % In m
Gap = 150e-9; % In m
Waveguide_Width = 450e-9; % In m
Comparison = Waveguide_Width/(2*Radius)+Gap/(2*Radius)+1;

% Other parameters
Freq_Center = 193.54; % In THz
EffRefrIndex = 2.38446; % For a 450 x 220 nm^2 waveguide cross section
% EffRefrIndex = 2.362; 

kk = 1;
for ii = 1:NumberOfIterations
    % Here we repeat this code snipet for each tunable parameter:
    % Set a parameter_min & a parameter_max.
   	CouplLength_min = 2e-6; % In m
    CouplLength_max = 5e-6; % In m
    % Generate a parameter_seed between parameter_min & parameter_max.
    CouplLength_seed = CouplLength_min + rand()*(CouplLength_max - CouplLength_min);
    
    % Set a parameter_min & a parameter_max.
   	Angle_min = 0.8;
    Angle_max = 1.2;
    % Generate a parameter_seed between parameter_min & parameter_max.
    Alpha_seed = Angle_min + rand()*(Angle_max - Angle_min);
    Alpha_2_seed = acos(Comparison-cos(Alpha_seed));
           
    
    % Initial guess (randomly assigned between X0_min & X0_max).
    X0 = [CouplLength_seed, Alpha_seed, Alpha_2_seed];
    opts = optimset('MaxIter',1e8, 'Display','none');
    % Gives X at the SIP. Sigma is the eigenvector coalescence parameter (the more similar the
    % eigenvalues, the smaller it is (tends to zero)).
    [X,sigma, exitflag] = fminsearch(@ASOW_Fab_SIP_Finding_Function, X0, opts);
    if sigma > 1e-3
        continue
    end
    if exitflag>0  % If the optimization is successful
        CouplLength = X(1);
        Alpha = X(2);
        Alpha_2 = X(3);
        
        [TransferMatrix] = ASOW_Fab_TransferMatrix_DirectionalCoupler (Freq_Center, CouplLength, Radius, Alpha, Alpha_2, EffRefrIndex);
        [Eigenvectors, EigenValues_X] = eig(TransferMatrix);
        Eigenvalues = diag(EigenValues_X);
        [Eigenvalues_X,index_X] = sort(Eigenvalues);
        Ksd_X = -log(Eigenvalues_X)./(1j);
        Ksd_pi = Ksd_X/pi;
        Eigenvectors_X = sort(Eigenvectors(:,index_X));
        
    end
    Sigma_EPD(kk) =sigma ;
    Eigenvalues_EPD(:,kk) = Eigenvalues_X;
    X_EPD(:,kk) = X;
    Ksd_pi_EPD(:,kk) = Ksd_pi;
    Eigenvectors_EPD(:,:,kk) = Eigenvectors_X;
    kk = kk + 1;
end
%
% Now we find the EPD by choosing the minimum sigma, the coalescence
% parameter. Check that the eigenvalues are the expected ones.
[Sigma_min, index] = min(Sigma_EPD);
% Sigma_min;
% Eigenvalues_min = Eigenvalues_EPD(:,index);
% Ksd_pi_min = Ksd_pi_EPD(:,index);
% Eigenvectors_min = Eigenvectors_EPD(:,:,ii);
% X_min = X_EPD(:,index);
% Gamma_min = X_min(2) + X_min(3);
%
% Minimum coupling length
[CouplLength_min, index_length] = min(X_EPD(1,:));
% Sigma_EPD(index_length)
% X_EPD(:,index_length)

%
% Plotting the Dispersion Diagram of the code
% With these parameters, add a dispersion diagram plotting code to check it
% is an EPD.
    % Plot minimum length
% CouplLength = X_EPD(1,index_length);
% Alpha = X_EPD(2,index_length);
% Alpha_2 = X_EPD(3,index_length);
    % Plot minimum sigma
CouplLength = X_EPD(1,index);
Alpha = X_EPD(2,index);
Alpha_2 = X_EPD(3,index);

 % Frequency Sweep
Freqmax = Freq_Center + 0.04;
Freqmin = Freq_Center - 0.04;
Freqsteps = 5000;
Freq = linspace(Freqmin, Freqmax, Freqsteps);

TransferMatrices_Vec = zeros(6,6,Freqsteps);
Eigenvectors = zeros(6,Freqsteps);
NumberOfUnitCells = 1;
for ii = 1:Freqsteps    
    TransferMatrices_Vec(:,:,ii) =  ASOW_Fab_TransferMatrix_DirectionalCoupler (Freq(ii), CouplLength, Radius, Alpha, Alpha_2, EffRefrIndex);
end
[EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
for ii = 1:Freqsteps
    Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
    Sigma(ii) = ASOW_Fab_SIP_Finding_Hyperdistance(TransferMatrices_Vec(:,:,ii));
end

%
% plot dispersion
    
Eigenvalues_plt = Eigenvalues'/pi;
Eigenvalues_plt_real = real(Eigenvalues_plt);
Eigenvalues_plt_imag = imag(Eigenvalues_plt);
Freq_plt = (ones(6,1)*Freq./Freq_Center)';

%
figure(1)
subplot(2,2,1)
hold on
plot(Eigenvalues_plt_real,Freq_plt,'.')
% plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
% line_color = ['b' 'g' 'y' 'c' 'm' 'r'];
line_color = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]	; 	[0.9290, 0.6940, 0.1250] ; [0.4940, 0.1840, 0.5560] ; [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]];  
for k = 1:1
    plot(Eigenvalues_plt_real(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
end
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
xlabel('Re[kd/\pi]')
ylabel('Normalized Frequency')
grid on
% set(gca,'FontSize',24,'FontName', 'Times New Roman');

subplot(2,2,2)
hold on
plot(Eigenvalues_plt_imag,Freq_plt,'.')
% plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
for k = 1:1
    plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
end
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
xlabel('Im[kd/\pi]')
ylabel('Normalized Frequency')
grid on
% set(gca,'FontSize',24,'FontName', 'Times New Roman');

subplot(2,2,3)
hold on
plot(Eigenvalues_plt_real,Eigenvalues_plt_imag,'.')
for k = 1:6
    plot(Eigenvalues_plt_real(:,k),Eigenvalues_plt_real(:,k),'.','Color',line_color(k,:))
end
axis([-1,1, -1,1])
xlabel('Re[kd/\pi]')
ylabel('Im[kd/\pi]')
grid on

subplot(2,2,4)
plot(Freq_plt(:,k),Sigma,'.')
% axis([-1,1, -1,1])
ylabel('Normalized Frequency')
grid on
ylabel('Sigma')

% Sigma_min = min(Sigma)


%% Now let's plot all the dispersion diagrams
close all; clc
kk = 1;
for jj = 1:NumberOfIterations
    if Sigma_EPD(jj) > 1e-3
        continue;
    end
    CouplLength = X_EPD(1,jj);
    Alpha = X_EPD(2,jj);
    Alpha_2 = X_EPD(3,jj); % ...
%     Radius = 6e-6; %m
%     Freq_Center = 193.54; % In THz
%     EffRefrIndex = 2.38446; % Dimensionless
  
    % For the angles
%     CosSum = cos(Alpha) + cos(Alpha_2);
%     Comparison = Waveguide_Width/(2*Radius)+Gap/(2*Radius)+1;
    Gap(jj) = [(2*Radius*(cos(Alpha)+cos(Alpha_2)-1)-Waveguide_Width)];
%     if  Gap < 150e-9
%         continue;
%     end

    TransferMatrices_Vec = zeros(6,6,Freqsteps);
    Eigenvectors = zeros(6,Freqsteps);
    NumberOfUnitCells = 1;
    for ii = 1:Freqsteps
        TransferMatrices_Vec(:,:,ii) =  ASOW_Fab_TransferMatrix_DirectionalCoupler (Freq(ii), CouplLength, Radius, Alpha, Alpha_2, EffRefrIndex);
    end
    [EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
    for ii = 1:Freqsteps
        Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
        Sigma(ii) = ASOW_Fab_SIP_Finding_Hyperdistance(TransferMatrices_Vec(:,:,ii));
    end
   
    %
    % plot dispersion
    
    Eigenvalues_plt = Eigenvalues'/pi;
    Eigenvalues_plt_real = real(Eigenvalues_plt);
    Eigenvalues_plt_imag = imag(Eigenvalues_plt);
    Freq_plt = (ones(6,1)*Freq./Freq_Center)';
    
    %
    figure(jj)
    kk = kk + 1;
    subplot(2,2,1)
    hold on
    plot(Eigenvalues_plt_real,Freq_plt,'.')
    % plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
    % line_color = ['b' 'g' 'y' 'c' 'm' 'r'];
    line_color = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]	; 	[0.9290, 0.6940, 0.1250] ; [0.4940, 0.1840, 0.5560] ; [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]];
    for k = 1:1
        plot(Eigenvalues_plt_real(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
    end
    axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
    xlabel('Re[kd/\pi]')
    ylabel('Normalized Frequency')
    grid on
    % set(gca,'FontSize',24,'FontName', 'Times New Roman');
    
    subplot(2,2,2)
    hold on
    plot(Eigenvalues_plt_imag,Freq_plt,'.')
    % plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
    for k = 1:1
        plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
    end
    axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
    xlabel('Im[kd/\pi]')
    ylabel('Normalized Frequency')
    grid on
    % set(gca,'FontSize',24,'FontName', 'Times New Roman');
    
    
    subplot(2,2,3)
    hold on
    plot(Eigenvalues_plt_real,Eigenvalues_plt_imag,'.')
    for k = 1:6
        plot(Eigenvalues_plt_real(:,k),Eigenvalues_plt_real(:,k),'.','Color',line_color(k,:))
    end
    axis([-1,1, -1,1])
    xlabel('Re[kd/\pi]')
    ylabel('Im[kd/\pi]')
    grid on
    
    subplot(2,2,4)
    plot(Freq_plt(:,k),Sigma,'.')
    % axis([-1,1, -1,1])
    ylabel('Normalized Frequency')
    grid on
    ylabel('Sigma')

end
    
    
    
    
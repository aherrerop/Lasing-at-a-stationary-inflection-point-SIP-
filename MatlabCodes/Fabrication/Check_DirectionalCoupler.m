%% Plot
clear all; close all; clc

Freq = 193.54; % In THz
RefrIndex_odd = returnEffIndexOdd(Freq, Freq, 1); % Input in THz
RefrIndex_even = returnEffIndexEven(Freq, Freq, 1);
SpeedLight = 3e8; % In m/s
Wavenumber_even = 2*pi*RefrIndex_even*Freq*1e12/SpeedLight;
Wavenumber_odd = 2*pi*RefrIndex_odd*Freq*1e12/SpeedLight;

Length = 1e-6*linspace(0.2,15,1e3);

for ii = 1:length(Length)
    [T1(:,:,ii)] = ModifiedNates_Transfer_Matrix_Coupled_Distributed(Length(ii), Wavenumber_even, Wavenumber_odd);
    Effective_Coupling(ii) = 1/abs(T1(1,3,ii));
    Effective_Transmission(ii) = abs(T1(1,2,ii))*Effective_Coupling(ii);
    Det(ii) = det(T1(:,:,ii)); % Should be 1 for all Lengths
end
%
figure(1)
hold on
plot(Length*1e6,Effective_Coupling,'b','linewidth',1.5)
plot(Length*1e6,Effective_Transmission,'k','linewidth',1.5)
plot(Length*1e6,Effective_Coupling.^2 + Effective_Transmission.^2,'r','linewidth',1.5)
hold off
grid on
legend('\kappa','\tau','\kappa^2 + \tau^2')

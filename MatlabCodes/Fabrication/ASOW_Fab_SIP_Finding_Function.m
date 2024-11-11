% This mock code is called in EPD_Finding_Mock_Code, where its input (the parameters) is optimized
% to minimize the hyperdistance sigma, calculated in
% EPD_Finding_Hyperdistance.
% In the bottom part of the code is where we can impose conditions on the values of the
% parameters or the position of the EPD in the BZ.
function[sigma] = ASOW_Fab_SIP_Finding_Function(X)
 CouplLength = X(1);
 Alpha = X(2);
 Alpha_2 = X(3); 
 Radius = 6e-6; % In m
 Freq_Center = 193.54; % In THz
 EffRefrIndex = 2.38446; % Dimensionless
%  EffRefrIndex = 2.362;

 Gap = 150e-9; % In m
 Waveguide_Width = 450e-9; % In m

% Get the transfer matrix
[TransferMatrix] = ASOW_Fab_TransferMatrix_DirectionalCoupler (Freq_Center, CouplLength, Radius, Alpha, Alpha_2, EffRefrIndex);
% Calculate the hyperdistance, sigma
[sigma] = ASOW_Fab_SIP_Finding_Hyperdistance(TransferMatrix);

% Enforce SIP

[~,EigenValues] = eigenshuffle(TransferMatrix);
Eigenvalues = -log(EigenValues)./(1j*pi);
for ii = 1:6
    if abs(real(Eigenvalues(ii))) < 0.1 | abs(real(Eigenvalues(ii))) > 0.9 
        sigma = 1e3;
    end
end
% We can enforce values of the solutions to be physical, we make a condition giving a huger sigma, which will be avoided when we minimize:
% 
if  CouplLength > 5e-5
    sigma = 1e3;
end

% For the angles
CosSum = cos(Alpha) + cos(Alpha_2);
Comparison = Waveguide_Width/(2*Radius)+Gap/(2*Radius)+1;
Gap = (2*Radius*(cos(Alpha)+cos(Alpha_2)-1)-Waveguide_Width);
if  Gap < 140e-9 | Gap > 160e-9
    sigma = 1e3;
end

end

function [TransferMatrix] = ASOW_Fab_TransferMatrix_DirectionalCoupler (Freq, CouplLength, Radius, Alpha, Alpha_2, EffRefrIndex)
    % Returns the Transfer Function of Scheuer's SOW when given the
    % parameters
    SpeedLight = 3e8; % In m/s
    % Freq is in THz
    pA = -(pi^2*EffRefrIndex*Radius*1e12)/SpeedLight; % Phase of Quarter circle 
    pB = -(4*pi*EffRefrIndex*Radius*Alpha*1e12)/SpeedLight; % Phase of twice alpha arc
    pD = -(4*pi*EffRefrIndex*Radius*Alpha_2*1e12)/SpeedLight; % Phase of twice alpha_2 arc

    PhaseSectionA = Freq * pA;
    PhaseSectionB = Freq * pB;
    PhaseSectionD = Freq * pD;
    
    % Tp1: Transfer Matrix from Z0+ to Zc- (from begining of the unit cell to first coupling point)
    % Tp2: Transfer Matrix from Zc+ to Zc'- (from first to second coupling
    % point)

    Tp1(1,1)= exp(1j*PhaseSectionA);
    Tp1(2,2)= exp(-1j*PhaseSectionA);
    Tp1(3,3)= exp(1j*PhaseSectionB);
    Tp1(4,4)= exp(-1j*PhaseSectionB);
    Tp1(5,5)= exp(1j*PhaseSectionA);
    Tp1(6,6)= exp(-1j*PhaseSectionA);
    
    Tp2(1,1)= exp(1j*PhaseSectionA);
    Tp2(2,2)= exp(-1j*PhaseSectionA);
    Tp2(3,3)= exp(1j*PhaseSectionD);
    Tp2(4,4)= exp(-1j*PhaseSectionD);
    Tp2(5,5)= exp(1j*PhaseSectionA);
    Tp2(6,6)= exp(-1j*PhaseSectionA);
     
    % T1, T2 are coupling matrices. They are directional coupling matrices.
    % Step 1: Get the refractive indexes  
    RefrIndex_odd = returnEffIndexOdd(Freq, Freq, 1);
    RefrIndex_even = returnEffIndexEven(Freq, Freq, 1);
    % Step 2: Get the wavenumbers
    Wavenumber_even = RefrIndex_even*Freq*1e12/SpeedLight;
    Wavenumber_odd = RefrIndex_odd*Freq*1e12/SpeedLight;
    % Step 3: Get the 4x4 transfer matrix of the coupling section
    TransferMatrix_Coupling = ModifiedNates_Transfer_Matrix_Coupled_Distributed(CouplLength, Wavenumber_even, Wavenumber_odd);
    T1(1:4,1:4) = TransferMatrix_Coupling;
    T2(3:6,3:6) = TransferMatrix_Coupling;
    % Step 4: Add the "phase accumulation" at the other port
    T1(5,5) = 1;
    T1(6,6) = 1;
    T2(1,1) = 1;
    T2(2,2) = 1;
    
TransferMatrix = T2*Tp2*T1*Tp1;

end
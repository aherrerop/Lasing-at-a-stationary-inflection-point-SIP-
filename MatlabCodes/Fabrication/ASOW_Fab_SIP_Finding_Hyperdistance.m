% This mock code is called in EPD_Finding_Mock_Code, through the function
% EPD_Finding_Function, through which sigma is minimized.
% The hyperdistance, sigma, is the sqrt of the angle between the
% eigenvectors. Minimizing sigma means the eigenvectors coalesce. It has to
% be adapted for each EPD. Here we show the example for an SIP.
function [sigma] = ASOW_Fab_SIP_Finding_Hyperdistance(TransferMatrix)

[V,~]=eig(TransferMatrix);% Calculate eigenvector matrix

C = combnk(1:6,3) ;% Get combination of 3 out of 6 

for index=1:length(C)
    % Get the eigenvectors of this combination
    V1=V(:,C(index,1));
    V2=V(:,C(index,2));
    V3=V(:,C(index,3));
    % Calculate the angles between them
    C12=(1-abs(V2'*V1)^2)^0.5;
    C13=(1-abs(V3'*V1)^2)^0.5;
    C23=(1-abs(V3'*V2)^2)^0.5;
    % Find the mean (group) these angles.
    HD(index)=mean(abs([C12 C13 C23]));

end
% The lowest combination is the one that shows the SIP. 
sigma=min(HD);
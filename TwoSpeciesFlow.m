function y = TwoSpeciesFlow(x,N,Npct, TauD1,TauD2,TauF,Tpct,TauT )
% fit to equation 3 in Köhler RH, Schwille P, Webb WW, Hanson MR.   (2000)  Journal of Cell Science 113 (22): 3921–3930.
% Advection with triplet state
%
% Copyright (c) 2024 GalbraithLab 2024 - JA Galbraith, CG Galbraith
% All rights reserved.
% see License.txt file for details
%
% CONSTANTS:

StructureParam= 5;    % 1.2715 / .2543 Leica FCS system
% FIT PARAMETERS:
% TauD =   diffusion time constant
% TauF =  flow time constant
% N = Number of particles
% TauT = Triplet time constant
% Tpct = percent triplet stat

%% Diffusion -  species 1

TeeOverTau1 = x ./ TauD1;

A1= 1 ./ (1+ TeeOverTau1) ;
A=A1';
Aflow1=A;
B1 = sqrt(1 ./ (1 + (1/StructureParam^2) .* TeeOverTau1)) ;
B = B1';
% Flow1 exponent
C1 =  (x ./ TauF) .^ 2;
C1=C1';
C = -C1 .* Aflow1; % note: A = 1 / ( 1 + t/TauD) for two species there are two Dtau (1 and 2)

D1=(Npct) .* A .*B .* exp(C);

%% Diffusion - second species
TeeOverTau2 = x ./ TauD2;

A1= 1 ./ (1+ TeeOverTau2) ;
A=A1';
Aflow2=A;
B1 = sqrt(1 ./ (1 + (1/StructureParam^2) .* TeeOverTau2)) ;
B = B1';
% Flow2 exponent
C1 =  (x ./ TauF) .^ 2;
C1=C1';
C = -C1 .* Aflow2; % note: A = 1 / ( 1 + t/TauD) for two species there are two Tau

SecondTerm=1-Npct;

D2=SecondTerm .* A .*B .*exp(C);  % with flow - the exp term

%% sum of mobilities with flow for each part
Dtot= (1/N) .* (D1 + D2);
G = Dtot  ; % G(tau) correlation function without triplet correction

%% Triplet state correction
S = (1 - Tpct + Tpct .* exp(-x ./ TauT));

%%
S=S';
y = (S .* G);
y = y';

end

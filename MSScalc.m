% MSScalc
%
% Copyright (c) 2024 GalbraithLab 2024 - JA Galbraith, CG Galbraith
% All rights reserved.
% see License.txt file for details
%
%  Calculates the Moment Scaling Spectrum (MSS) - references:
%   Ferrari, R., Manfroi, A. J. & Young, W. R. (2001) Physica D 154, 111–137.
%   Ewers, H, Smith, A, Sbalzarini, IF, Lilie, H, Koumoutsakos, P, Helenius,A.  (2005)  PNAS  102 (42) 15110-15115.
%
%  Uses polyfitZero.m which is a function by Mark Mikofski
%  and can be obtained at: https://www.mathworks.com/matlabcentral/fileexchange/35401-polyfitzero
%
% Input data set is a mat file containing trajectory info generated from TrackIt program results
% The mat file contains two variables:  "TrackItData" - Trajectory ID Number, X, Y, Step #, ROI (not used here)
% and "file" - name of the dataset which is used for labeling output files
%
% The output consists of three graphs: histograms of the MSS slope and diffusion coefficient and a scatter plot of the diffusion coefficient versus the MSS slope.
% The scaling exponents for each moment and the slope of the MSS curve are saved as .mat and .csv files.  
%
%

% Select file
[NewFile,path] = uigetfile('mat');
load(strcat(path,NewFile))
TrackItBak=TrackItData;  
TotalNumberOfTrajectories=max(TrackItData(:,1));
TotalNumberOfJumps=size(TrackItData,1);

% Enter track parameters
prompt = {'Minimum Number of Steps','Maximum Number of Steps','Number of moments','Lag divider - MOSAIC 3, Saxton 4','Regression Coeff cutoff','Max Diffusion Coeff (um^2/sec)'};
dlg_title = 'Moment parameters';
num_lines=1;
def = {'15','125','7','3', '.5','50'}; 
answer = inputdlg(prompt,dlg_title,num_lines,def);
MinTrkLen=str2double(char(answer(1,1)));
MaxTrkLen=str2double(char(answer(2,1)));
NumberOfMoments=str2double(char(answer(3,1)));
LagSize= str2double(char(answer(4,1))); 
RegressionCutoff=str2double(char(answer(5,1))); 
MaxDiffusionCoeff=str2double(char(answer(6,1))); 

Timebase=0.008;  % for 8ms acquisition time 63 frames = 1/2 second, 125 frames = 1 sec

% Find trajectories that are outside of length constraints  MinTrkLen < TRAJECTORY LENGTH < MaxTrkLen
TrajectoryIndex = zeros(TotalNumberOfTrajectories, 2);
TooShortTrajectories = [];
TooLongTrajectories = [];

indicesCell = accumarray(TrackItData(:, 1), (1:size(TrackItData, 1))', [], @(x) {x});

for i = 1:TotalNumberOfTrajectories
    indices = indicesCell{i};
    if isempty(indices)
        continue;
    end
    TrajectoryIndex(i, 1) = indices(1);
    TrajectoryIndex(i, 2) = indices(end);
    
    trajLen = TrajectoryIndex(i, 2) - TrajectoryIndex(i, 1) + 1;
    
    if trajLen < MinTrkLen
        TooShortTrajectories(end + 1, 1) = i; 
    end
    
    if trajLen > MaxTrkLen - 1
        TooLongTrajectories(end + 1, 1) = i;
    end
end

AllTrajectoryLengths=(TrajectoryIndex(:,2) - TrajectoryIndex(:,1)+1) ; % Length of all trajectory found by TrackIt

NumTooShortTrajec=size(TooShortTrajectories,1);
NumTooLongTrajec=size(TooLongTrajectories,1);

% Remove out of range trajectories
OutOfRangeNdx = ismember(TrackItData(:,1), TooShortTrajectories) | ismember(TrackItData(:,1), TooLongTrajectories);
TrackItData(OutOfRangeNdx, :) = [];

% all trajectories shorter than MinTrkLen and longer than MaxTrkLen have been removed but... we now need to reindex the start and end positions of the tracks which are not sequentially numbered anymore
TotalNumberOfLengthFilteredJumps=size(TrackItData,1);

[TrackItIndex(:,1), TrackItIndex(:,2)]=unique(TrackItData(:,1),'first'); % first position of each trajectory: TrkNum (1) TrackIt track number, (2) starting array row
NumberOfGoodLengthTrajectories=size(TrackItIndex,1);
GoodLengthTrajectoriesTrackItID=unique(TrackItData(:,1));
MnD={NumberOfGoodLengthTrajectories}; % cell array with all moment curves for use in later calculations "Mean n Displacements"
MnD{1}=strcat('Min Trk Len:',num2str(MinTrkLen),32,'Max Trk Len:',num2str(MaxTrkLen),32,'LagSize:',num2str(LagSize),32,'Pct curve analyzed:','Regression cutoff',RegressionCutoff,'Max Diff Coeff',MaxDiffusionCoeff);

% Fit of log-log moment curves for each trace - overwritten each time
Ewers_MSD_LinearDiffusion=zeros(NumberOfGoodLengthTrajectories,4,NumberOfMoments);

h = waitbar(0,''); % Intialize progress bar
for i =1: NumberOfGoodLengthTrajectories

 PctProg=i / NumberOfGoodLengthTrajectories; % NumberOfTrajectories
     WaitBarLabel=strcat('Calculating moments for trajectory:',32,num2str(i),strcat(32,'/',32,num2str(NumberOfGoodLengthTrajectories))); % allows for spaces between strings
     if mod(i,1) == 0  % Update progress bar
        waitbar(PctProg,h,WaitBarLabel);
    end

Trkstart=TrackItIndex(i,2);
if i<NumberOfGoodLengthTrajectories
    Trkend=TrackItIndex(i+1,2)-1;
else
    Trkend=size(TrackItData,1);
end
TrkLength=Trkend-Trkstart+1;
Trajectory= TrackItData(Trkstart:Trkend,:);% TrackItData - x column 2, y column 3

numDataPts = size(Trajectory,1); %# number of data points
NumberOfDeltaT = floor(numDataPts/LagSize); % # for MSD, dt should be up to 1/4 of number of data points (Saxton), MOSAIC uses 1/3 and TrackIt uses 60%-90%
dataScaled=Trajectory;
dataScaled(:,4)=Timebase*(Trajectory(:,4)-Trajectory(1,4)); % start timesteps at 1

% # calculate msd for all deltaT's
Moments=zeros(NumberOfDeltaT,NumberOfMoments+2); % last two columns are lag (integer and time)

for j= 1:NumberOfMoments % calculate moments
 
    for dt = 1:NumberOfDeltaT  % different size steps [1]: 1-2, 2-3,  [2]: 1-3, 2-4, [3]: 1-4, 2-5,
        % calculate each step distance
        dx= dataScaled(1+dt:end,2) - dataScaled(1:end-dt,2);
        dy=dataScaled(1+dt:end,3) - dataScaled(1:end-dt,3);
        distance=sqrt(dx.^2+dy.^2);

        % raise each distance to power of moment
        DistanceRaisedToPower = distance .^j;
        Moments(dt,j)=mean(DistanceRaisedToPower); 
        Moments(dt,NumberOfMoments+1)=numDataPts-length(DistanceRaisedToPower);
        Moments(dt,NumberOfMoments+2)=Timebase*Moments(dt,NumberOfMoments+1);
    end
end % end of j th  moment calculations
MnD{i+1}=Moments;
end
MnD=MnD';
delete(h);

% Calculate diffusion coefficient (slope) from 2nd moment (MSD)
h = waitbar(0,''); % Intialize progress bar
for i=2:NumberOfGoodLengthTrajectories+1 % start at 2 since 1 is text info

    PctProg=i / NumberOfGoodLengthTrajectories; % NumberOfTrajectories
    WaitBarLabel=strcat({'Analyzing MSD curve:  '},num2str(i),strcat(32,'/',32,num2str(NumberOfGoodLengthTrajectories))); 
    if mod(i,1) == 0  % Update progress bar
        waitbar(PctProg,h,WaitBarLabel);
    end

    for k=1:NumberOfMoments
    
     % Ewers method of log-log plot - diffusion coefficient is exp(intercept) because of log-log
        [Linearfit, goodness]=polyfit(log(MnD{i,1}(:,9)),log(MnD{i, 1}(:,k)),1) ; % first order equation
        Ewers_MSD_LinearDiffusion(i-1,1,k)=0.25*exp(Linearfit(2)); % diffusion coefficient
        Ewers_MSD_LinearDiffusion(i-1,2,k)=Linearfit(1); % slope
        Ewers_MSD_LinearDiffusion(i-1,3,k)=Linearfit(2);  % intercept (could also back out from diffusion coeff)
        Ewers_MSD_LinearDiffusion(i-1,4,k)=goodness.rsquared;

    end  % end of moment calculation
end  % end of calculation for each trajectory

delete(h);  % delete progress bar

% apply input constraints
EwersDiffusionRegressionCoeff=Ewers_MSD_LinearDiffusion(:,4,2);
PoorRegression=find(EwersDiffusionRegressionCoeff<RegressionCutoff);
TooBigDIdx=find(Ewers_MSD_LinearDiffusion(:,1,2)>MaxDiffusionCoeff);

%Using slope of log-log plot
for n=1:NumberOfMoments
    MSSdataEwers(:,n)=Ewers_MSD_LinearDiffusion(:,2,n);
end
MSSdataEwers(:,NumberOfMoments+1)= Ewers_MSD_LinearDiffusion(1:NumberOfGoodLengthTrajectories,1,2);

%check for non monotonically increasing curves MSS curves
maxPossibleDownturns = NumberOfGoodLengthTrajectories * (NumberOfMoments - 1); % preallocate for speed and then trim
Edex = zeros(maxPossibleDownturns, 1);
counter = 0;

for p = 1:NumberOfGoodLengthTrajectories
    for m = 2:NumberOfMoments
        if MSSdataEwers(p, m) < MSSdataEwers(p, m - 1)
            counter = counter + 1;
            Edex(counter, 1) = p;
        end
    end
end

Edex = Edex(1:counter);

if ~isempty(Edex)
Edex(1)=[];
end

MonoDex=unique(Edex);

Edex=cat(1,Edex, PoorRegression, TooBigDIdx);
Edex=unique(Edex);  % remove duplicates
MSSdata=MSSdataEwers;
MSSdata(:,NumberOfMoments+1)= Ewers_MSD_LinearDiffusion(:,1,2);  % add diffusion coefficients to array
MSSdata(:,NumberOfMoments+2)=  GoodLengthTrajectoriesTrackItID;% add track number for retrival with future analysis
MSSdata(:,NumberOfMoments+3)=  AllTrajectoryLengths(GoodLengthTrajectoriesTrackItID);% add track length 
if ~isempty(Edex)
MSSdata(Edex,:)=[];
end

% MSSdata: columns 1-NumberofMoments, MSD Ewers diffusion coefficient, track number, track length
Xmss=1:NumberOfMoments;

% fit curve through origin ala Ewers using original MSSdataset 
for v=1:size(MSSdata,1)    
    ZeroFit=polyfitZero(Xmss', MSSdata(v,1:NumberOfMoments)',1); % polyfitZero is a function by Mark Mikofski 
    % https://www.mathworks.com/matlabcentral/fileexchange/35401-polyfitzero
    MSSslope(v,1)=ZeroFit(1,1);
end

% Smss: 1) slope, 2) diffusion coefficient, 3) track number, 4) track length
Smss=MSSslope(:,1); % Slope mss
Smss(:,2:4)=MSSdata(:,NumberOfMoments+1:end);

%%  Output ...

fig=figure;
histogram(Smss(:,1),20,'Normalization','percentage','BinWidth',0.05);
title('Ewers MSS slope');
xlabel('Ewers MSS slope');
ylabel('Percentage');
subtitle(file,"FontSize",10,'Interpreter','none') % filename
FigFile=strcat(file,'_MSSslopeHistogram');
print(fig,FigFile,"-dtiff")
close;

fig=figure;
scatter((Smss(:,2)),Smss(:,1));
title('Diffusion coefficient vs MSS slope');
subtitle(file,"FontSize",9,'Interpreter','none') % filename
xlabel('Diffusion coefficient um^2 / sec');
ylabel('MSS slope');
ylim([0 0.8])
FigFile=strcat(file,'_DiffCoeff_MSSslopeScatterPlot');
print(fig,FigFile,"-dtiff")
close;

fig=figure;
histogram(Smss(:,2),50,'Normalization','percentage');
title('Diffusion coefficient');
xlabel('Ewers 2nd moment (log-log) diffusion coefficient (um^2/sec)');
subtitle(file,"FontSize",9,'Interpreter','none') % filename
ylabel('Percentage');

FigFile=strcat(file,'_DiffusionHistogram');
print(fig,FigFile,"-dtiff")
close;

save(strcat(file,'_MomentScalingSpectrum.mat'),'Smss','MSSdata','-mat');
writematrix(MSSdata,strcat(file,'_MomentScalingSpectrumData.csv'));
writematrix(Smss,strcat(file,'_MSSslope.csv'));

 clear;
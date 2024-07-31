% TransportRatio
%
% Copyright (c) 2024 GalbraithLab 2024 - JA Galbraith, CG Galbraith
% All rights reserved.
% see License.txt file for details
%
%  Uses trapezoidal rule to determine area under each half of fluorescent intensity trace.  
%  Midpoint is considered peak of signal - i.e. the photoactivation spot.
%  Input data is stored in Excel file with Row1 the experimental name and line direction
%  Odd columns are perpendicular traces (perp)
%  Even columns are parallel traces (para)
%  "left" side of trace is cell leading edge
%
%  Reads either an excel or txt file
%  Output: ".mat" file with experimental names, traces, calculated ratios for traces / spline fit in both directions, and spline ratios and group statistics (mean and std)

clear;

[NewFile,path] = uigetfile('*');
[~,~,anExt]=fileparts(NewFile);
idx = find(strcmp(anExt,{'.xlsx','.txt'})) ;
if isempty(idx)
    uiwait(msgbox('Unsupported file type. Choose .xlsx or .txt file','Error','error'));
else
    FullPath=strcat(path,NewFile);
    FileExt=NewFile((max(strfind(NewFile,'.')))+1:end);
    FileRoot=NewFile(1:(max(strfind(NewFile,'.')))-1);
    switch FileExt
        case 'xlsx'
            [Traces,ExpNames]=xlsread(FullPath);
        case 'txt'
            TempTxt=importdata(FullPath);
            ExpNames=TempTxt.colheaders;
            Traces=TempTxt.data;
    end
end

NumExp=size(ExpNames,2);

% create array to store "raw" and "spline" results
% row 1 is raw, row 2 is spline / columns match with ExpNames
TransportRatioCalculated=zeros(2,NumExp);
TransportRatioParallel=zeros(2,NumExp/2);
TransportRatioPerpendicular = zeros(2, NumExp / 2);

% perp = line profile perpendicular to cell edge
% para = line profile parallel to cell edge
% grab last 4 characters to determine line direction
% ExpNames(r)

for r=1:NumExp

%  copy column vector to variable "X" for analysis
X=Traces(:,r);
X = X(all(~isnan(X),2),:);  %  strip out NaN in case traces are different lengths

% finds direction of line assuming perp or para as descriptor and no other uses of letter p in name
% temp=char(ExpNames{r});
% LineDirection=temp((strfind(temp,'p'):(strfind(temp,'p'))+3));

[Xmin,minIndex]=min(X); % find minimum value to use for baseline subtraction
Xshift=X-X(minIndex);  % subtract baseline
[Xmax,MaxIndex]=max(Xshift); % max value will be midpoint
Xnorm=Xshift/Xmax; % baseline subtracted and peak normalized

Startpt=1;
Endpt=size(Xnorm,1);

% Check to see if the peak is centered in the trace
RightWidth=Endpt-MaxIndex;
LeftWidth=MaxIndex-1;
if RightWidth ~= LeftWidth
    % Determine the starting or ending point based on which side is shorter
    if RightWidth < LeftWidth
        Startpt = MaxIndex - RightWidth;
        % disp([ExpNames{r}, ' Peak not centered - right side is shorter']);
    else
        Endpt = MaxIndex + LeftWidth;
        % disp([ExpNames{r}, ' Peak not centered - left side is shorter']);
    end
end

LeftArea=trapz(Xnorm(Startpt:MaxIndex-1));
RightArea=trapz(Xnorm(MaxIndex+1:Endpt));
Ratio=LeftArea/RightArea;

%% Use a spline fit to better approximate peak

SplineStep=.1;
a = Startpt:Endpt;
b = Xnorm(Startpt:Endpt); 
c = Startpt:SplineStep:Endpt;

Xspline=spline(a,b,c);

[d,e]=max(Xspline);
[Ymax,Yidex]=max(Xspline); % max value will be used as midpoint

[XsplineMin,XsplineMinIndex]=min(Xspline);
XsplineShift=Xspline-Xspline(XsplineMinIndex);
[XsplineMax,XsplineIndex]=max(XsplineShift); % max value will be midpoint
XsplineNorm=XsplineShift/XsplineMax; % baseline subtracted and peak normalized

XsplineStartpt=1;
XsplineEndpt=size(XsplineNorm,2);

% Check to see if the peak is not centered
XsplineRightWidth=XsplineEndpt-XsplineIndex;
XsplineLeftWidth=XsplineIndex-1;
if XsplineRightWidth ~= XsplineLeftWidth
    % Determine the starting or ending point based on which side is shorter
    if XsplineRightWidth < XsplineLeftWidth
        XsplineStartpt = XsplineIndex - XsplineRightWidth;
        % disp([ExpNames{r}, ' Peak not centered - right side is shorter']);
    else
        XsplineEndpt = XsplineIndex + XsplineLeftWidth;
        % disp([ExpNames{r}, ' Peak not centered - left side is shorter']);
    end
end

XsplineLeftArea=trapz(XsplineNorm(XsplineStartpt:XsplineIndex-1));
XsplineRightArea=trapz(XsplineNorm(XsplineIndex+1:XsplineEndpt));
XsplineRatio=XsplineLeftArea/XsplineRightArea;

TransportRatioCalculated(1,r)=Ratio;
TransportRatioCalculated(2,r)=XsplineRatio;

end

% Perpendicular Data Set
Ndx = 1:2:NumExp;
TransportRatioPerpendicular(1, :) = TransportRatioCalculated(1, Ndx);
TransportRatioPerpendicular(2, :) = TransportRatioCalculated(2, Ndx);
AdvectionRatioPerpendicularNames = ExpNames(1,Ndx);
% Stats
PerpendicularMean(1,1)=mean(TransportRatioPerpendicular(1,:));
PerpendicularMean(2,1)= mean(TransportRatioPerpendicular(2,:));
PerpendicularStd(1,1)=std(TransportRatioPerpendicular(1,:));
PerpendicularStd(2,1)= std(TransportRatioPerpendicular(2,:));

% Parallel Data Set
Ndx = 2:2:NumExp;
TransportRatioParallel(1, :) = TransportRatioCalculated(1, Ndx);
TransportRatioParallel(2, :) = TransportRatioCalculated(2, Ndx);
TransportRatioParallelNames = ExpNames(1,Ndx);
% Stats
ParallelMean(1,1)=mean(TransportRatioParallel(1,:));
ParallelMean(2,1)= mean(TransportRatioParallel(2,:));
ParallelStd(1,1)=std(TransportRatioParallel(1,:));
ParallelStd(2,1)= std(TransportRatioParallel(2,:));

% Format for boxplots
PerpendicularRawBoxPlot=TransportRatioPerpendicular(1,:)';
PerpendicularSplineBoxPlot=TransportRatioPerpendicular(2,:)';
ParallelRawBoxPlt=TransportRatioParallel(1,:)';
ParallelSplineBoxPlt=TransportRatioParallel(2,:)';

save(strcat(path,FileRoot,'_TransportRatio.mat'),'ExpNames','Traces','TransportRatioCalculated','PerpendicularRawBoxPlot','PerpendicularSplineBoxPlot','ParallelRawBoxPlt','ParallelSplineBoxPlt','ParallelMean','ParallelStd','PerpendicularMean','PerpendicularStd');

clear;
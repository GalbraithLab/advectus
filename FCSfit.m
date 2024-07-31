% FCSfit 
%
% Copyright (c) 2024 GalbraithLab 2024 - JA Galbraith, CG Galbraith
% All rights reserved.
% see License.txt file for details
%
% Fits raw FCS curve to diffusion and advection using Equation 3 Köhler et. al. (2000)  Journal of Cell Science 113 (22): 3921–3930.
% Dependencies: Uses the Curvefit toolbox and TwoSpeciesFlow.m to define fitting equation
%  ProTip:  [fList, pList] = matlab.codetools.requiredFilesAndProducts('FCSfit.m') will list this information
%
%  The calculated FCS values are stored in FCSdata array
%  The initial seed values are stored in FCSseeds array
%
%  Reads an Excel spreadsheet containing the Leica autocorrelation data where
%  Column A is time (ms), and the others are the data traces.  The first two rows
%  (formatted as text variables) are the experiment identifiers - row 1 trial/run number, and row 2 experimental date
%  The sheet name contains the experimental condition/treatment and is added to the filenames and graph titles.  
% The user can choose the desired tab to analyze if multiple sheets are present (for different conditions).  
% 
%The output consists of graphs (saved as tiffs) for each trace, with the model fit overlayed on the raw data and a plot of the residuals.
% A .mat file and separate CSV files contain the data, model fit, seeds, limits, residuals, and std of residual.
% The FCSdata array contains the fitted parameters: Trace number, regression coefficient, number of molecules,
% pct of species A, Deq species1, Deq species2, flow, triplet pct, and triplet time constant.
%
%
if (~exist ('Experiment','var'))

[FileName,path]=uigetfile({'*.xlsx'});
FullPath=strcat(path,FileName);
FileExt=FileName((max(strfind(FileName,'.')))+1:end);
sheets = sheetnames(FullPath);  % get names of all sheets in file
select=cell(1,size(sheets,1));
for i =1:size(sheets,1)
    select{1,i}=sheets(i);
end

if size(sheets,1)>1 % select from multiple sheets
    SheetList = select;
    [SheetIdx] = listdlg('ListString', SheetList,'SelectionMode', 'Single', 'PromptString', 'Select Sheet', 'Initialvalue', 1,'Name', 'EXPERIMENTS');
    ExperimentList=sheets(SheetIdx);
    NumExperiments=size(ExperimentList,1);
else % only one sheet
    ExperimentList=sheets;
    NumExperiments=1;
end

Experiment=ExperimentList(1);
PlotTitle=Experiment;
FCS_Table = readtable(FullPath,'Sheet',Experiment,'ReadVariableNames',true);
ColumnNames=FCS_Table.Properties.VariableNames;
RunNames=ColumnNames;
RunNames(1)=[];
XLSdata=table2array(FCS_Table);
XLSdata(1,:)=[]; %remove NaN from 1st row
XLStime=XLSdata(:,1); % create separate time array - will be X data
XLSdata(:,1)=[]; % remove first column with time
NumTraces=size(XLSdata,2); % number of trials
% Summary data for each trace:  TraceNumber, Deq, Flow, SystemCorr, TauSystem, Rsq
FCSdata=zeros(NumTraces,9);
FCSseeds=zeros(NumTraces,7);
FCSfitLimits=NaN(NumTraces,14);
FCSresid=NaN(size(XLStime,1),NumTraces);
FCSstdResid=NaN(size(XLStime,1),NumTraces);
FCSModelFit=zeros(size(XLStime,1),NumTraces);

end

% FCS microscope specific parameters
LateralRadius=0.2543;
Zo=1.2715;

prompt = {'D1 min','D1 max','D2 min','D2 max','Flow min','Flow max','Starting N avg points'};
dlg_title = 'Search bounds';
num_lines=1;
def = {'4.5','30', '0.01','2.5','0.0001','10','5'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
D1lower= str2double(char(answer(1,1)));
D1upper=str2double(char(answer(2,1))); 
D2lower= str2double(char(answer(3,1)));
D2upper=str2double(char(answer(4,1))); 
Flowlower= str2double(char(answer(5,1)));
Flowupper=str2double(char(answer(6,1)));
Navg=str2double(char(answer(7,1)));

TauD1upper=(1000*LateralRadius^2)/(4*D1upper);
TauD1lower=(1000*LateralRadius^2)/(4*D1lower);
TauD2upper=(1000*LateralRadius^2)/(4*D2upper);
TauD2lower=(1000*LateralRadius^2)/(4*D2lower);
TauFlowupper=2000*LateralRadius/Flowupper;
TauFlowlower=2000*LateralRadius/Flowlower;

  h = waitbar(0,''); % Intialize progress bar
for TraceNumber =1:  NumTraces 

            PctProg=TraceNumber/NumTraces;
            WaitBarLabel=strcat({'Calculating G(t) using 2 species with flow model   Run:  '},num2str(TraceNumber),strcat('/',num2str(NumTraces))); % allows for spaces between strings
            if mod(i,1) == 0  % Update progress bar
                waitbar(PctProg,h,WaitBarLabel);
            end
% Starting values
    FirstN=mean(XLSdata(1:Navg,TraceNumber));  % use average of first Navg points

    Nstart=(1/FirstN)*1;  
    
    Npctstart=1;                          
    TauD1start=1.6167;  % 10 um^2/sec
    TauD2start=161.671; % 0.1 um^2/sec
    TauFstart=508.6; %  1um/sec
    TpctStart=0.15; 
    TauTStart=0.05; 

   % Upper and lower bounds
   N_lower = (1/FirstN)-0.055;
   N_upper= (1/FirstN)+0.055;

    Npct_upper=1; % max 
    Npct_lower=.5; % 50-50 

    TauD1_lower = TauD1upper;
    TauD1_upper = TauD1lower;

    TauD2_lower =   TauD2upper; 
    TauD2_upper =   TauD2lower;       

    TauF_lower =TauFlowupper ; 
    TauF_upper = TauFlowlower; 

    Tpct_lower = 0.0;
    Tpct_upper = 0.35 ; 

    TauT_lower = 0;
    TauT_upper =0.5;

    ft = fittype('TwoSpeciesFlow(x,N1,N2,TauD1,TauD2,TauF,Tpct,TauT)');
   
   DiffMinDelta=.05;  
   DiffMaxDelta=.05;  
   FuncTol=1e-7;  

   MaxIterations=100000;
   MaxFuncEvals=20000;

    [cfit,goodness]= fit(XLStime,XLSdata(:,TraceNumber),ft,"StartPoint", [Nstart Npctstart TauD1start TauD2start TauFstart TpctStart TauTStart ], ...
        "Algorithm",'Trust-Region',"MaxIter",MaxIterations,"MaxFunEvals",MaxFuncEvals,"DiffMinChange",DiffMinDelta,"DiffMaxChange",DiffMaxDelta,"TolFun",FuncTol,"TolX",1e-6,  ...
        'Lower', [N_lower,Npct_lower,TauD1_lower,TauD2_lower,TauF_lower,Tpct_lower,TauT_lower], "Upper", [N_upper ,Npct_upper,TauD1_upper, TauD2_upper, TauF_upper , Tpct_upper, TauT_upper]);

    CV=coeffvalues(cfit); % coefficients are listed ALPHABETICALLY -  1) N, 2) Npct, 3) D1fit, 4) D2fit, 5) Flowfit, 6) TauSys, 7) Tpct
    Nfit=CV(1);
    Npctfit=CV(2);

    D1fit=CV(3);
    D1eq=(LateralRadius^2)/(4*D1fit);

    D2fit=CV(4);
    D2eq=(LateralRadius^2)/(4*D2fit);

    FlowFit=CV(5);
    Flow= 2* LateralRadius / FlowFit; % from Schwille - distance to cross beam = 2 * radius

    TauT=CV(6);  % Remember variables are returned alphabetically even though fed in different order 
    Tpct=CV(7);    % so "Ta" comes before "Tp" even though Tp was passed first

    RegCoeff=goodness.adjrsquare; % adjusted r^2

    yfit=TwoSpeciesFlow(XLStime,Nfit,Npctfit,D1fit,D2fit,FlowFit,Tpct,TauT); % generate fit curve with estimated parameters
    FCSresid(:,TraceNumber)=XLSdata(:,TraceNumber)-yfit;
    ResidStdev=std(FCSresid(:,TraceNumber));  % standard deviation of residuals
    FCSstdResid(:,TraceNumber)=FCSresid(:,TraceNumber)/ResidStdev;
    fig=figure('Visible','off');

    tiledlayout(8,3)
    nexttile([6,3])
    set(fig,'position',[100,100,1100,950])
    semilogx(XLStime,yfit)
    hold on;
    semilogx(XLStime,XLSdata(:,TraceNumber))
    yline(0,'--');
    hold off;
    CellRun=FCS_Table.Properties.VariableNames{TraceNumber+1};
    year=FCS_Table(1,TraceNumber+1);
    year=table2array(year);
    yr=num2str(year);
    CellRun=strcat(yr,':',CellRun);

    title(strcat('2 Species Deq with flow and Triplet state correction:',32,32,char(PlotTitle),',',32,'Trace:',num2str(TraceNumber),',',32,CellRun,32,32,'r^2=',num2str(RegCoeff)),'Interpreter','none');

    subtitletxt=strcat('D1eq =',32,num2str(D1eq*1000,5),32,'µm^2/sec',32,':',32,num2str(Npctfit*100,4),'%',32,'D2eq=',32,num2str(D2eq*1000,5),32,'µm^2/sec',32,32,'Flow=',...
        32,num2str(Flow*1000,5),32,'µm/sec',32,32,'Triplet %=',32,num2str(Tpct,4),32,'Triplet Tau=',32,num2str(TauT,3),32,'ms',32,32,'N=',32,num2str(Nfit+ 0,4),32,'molecules');
    subtitle(subtitletxt);
    grid on;
    nexttile([2,3])
    semilogx(XLStime,FCSresid(:,TraceNumber))

    ResLimit=1.25*max(abs(FCSresid(:,TraceNumber)));
    ylim([-ResLimit,ResLimit]);
    yline(0,'-');
    st1=yline(ResidStdev,'--');
    st2=yline(-ResidStdev,'--');
    legend(st1,'±1 stdev', 'FontSize',14);
    legend('boxoff');
    title('Raw residuals');
    ax=gca;
    ax.TitleHorizontalAlignment="right";

    grid on;

%%  
FigFile=strcat(PlotTitle,"_2SpeciesFlow",num2str(TraceNumber));
print(fig,FigFile,"-dtiff")

close(fig);
FCSModelFit(:,TraceNumber)=yfit;

FCSdata(TraceNumber,1)=TraceNumber;
FCSdata(TraceNumber,2)=RegCoeff;
FCSdata(TraceNumber,3)=Nfit;
FCSdata(TraceNumber,4)=Npctfit;

FCSdata(TraceNumber,5)=D1eq*1000; % convert time from milliseconds to sec
FCSdata(TraceNumber,6)=D2eq*1000;

FCSdata(TraceNumber,7)=Flow*1000;
FCSdata(TraceNumber,8)=Tpct;
FCSdata(TraceNumber,9)=TauT;

% Seed values for curve fit
FCSseeds(TraceNumber,1)= Nstart;
FCSseeds(TraceNumber,2)= Npctstart;
FCSseeds(TraceNumber,3)= TauD1start;
FCSseeds(TraceNumber,4)= TauD2start;
FCSseeds(TraceNumber,5)= TauFstart;
FCSseeds(TraceNumber,6)= TpctStart;
FCSseeds(TraceNumber,7)= TauTStart;

% Bounding limits
FCSfitLimits(TraceNumber,1)=TraceNumber;
FCSfitLimits(TraceNumber,2)=N_lower;
FCSfitLimits(TraceNumber,3)=N_upper;
FCSfitLimits(TraceNumber,4)=Npct_upper;
FCSfitLimits(TraceNumber,5)=Npct_lower;
FCSfitLimits(TraceNumber,6)=TauD1_lower;
FCSfitLimits(TraceNumber,7)=TauD1_upper;
FCSfitLimits(TraceNumber,8)=TauD2_lower;
FCSfitLimits(TraceNumber,9)=TauD2_upper;
FCSfitLimits(TraceNumber,10)=TauF_lower;
FCSfitLimits(TraceNumber,11)=TauF_upper;
FCSfitLimits(TraceNumber,12)=Tpct_lower;
FCSfitLimits(TraceNumber,13)=Tpct_upper;
FCSfitLimits(TraceNumber,14)=TauT_lower;
FCSfitLimits(TraceNumber,15)=TauT_upper;

end
delete(h);

% save computed parameters both as Matlab variables and CSV for import into Excel
FitFile=strcat(PlotTitle,"_2SpeciesFlow.mat");
save(FitFile, "FCSdata","FCSseeds","FCSfitLimits","FCSresid","FCSstdResid","FCSModelFit","XLStime","XLSdata");

writematrix(FCSdata,strcat(PlotTitle,'_FCS2SpeciesFlowData.csv'));
writematrix(FCSseeds,strcat(PlotTitle,'_FCS2SpeciesFlowSeeds.csv'));
writematrix(FCSfitLimits,strcat(PlotTitle,'FCS2SpeciesFlowfitLimits.csv'));
writematrix(FCSresid,strcat(PlotTitle,'FCS2SpeciesFlowResiduals.csv'));
writematrix(FCSstdResid,strcat(PlotTitle,'_2SpeciesFlowFCSstdResiduals.csv'));
writematrix(FCSModelFit,strcat(PlotTitle,'_2SpeciesFlowFCSmodelFit.csv'));

clear;
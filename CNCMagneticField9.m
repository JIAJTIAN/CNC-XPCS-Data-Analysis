%% Magnetic Field Under XPCS 9.1
% 1. Import Equilibrium test, Repetition test 
% 2. Imported Dose study
% The 8% data is problematic. Let's avoid using it

clc,clear
close all

%% Read Data
BetaArray = [0.138, 0.135, 0.133, 0.129, 0.128, 0.129, 0.129, 0.130];
AttArray = {'0.007','0.036','0.190','1.000'};
EpsArray = {'1.334','9.994','99.994'};

% Generate Dose Array
DoseIndex = 1;
for AttIndex = 1:length(AttArray)
    for EpsIndex = 1:length(EpsArray)
        Att = AttArray{AttIndex};
        Eps = EpsArray{EpsIndex};
        DoseArray(DoseIndex) = str2num(Att)*str2num(Eps);
        DoseIndex = DoseIndex +1;
    end
end
DoseArrayTotal = unique(DoseArray);
DoseLolim = 5;% Lower dose limit INDEX;
DoseUplim = 7; % Higher dose limit INDEX

% Threshold = 0.2;
% RefinedLength = 6;
EqDiff = 0.01;
RpDiff = 0.01;
RpNum = 2; % used in the equil and rpt check, this is the number of times the experiment was repeated 
PassRate = 0.95;
% SmoothnessThreshold = 0;
PlotIndex = 1;

% Folder names for Jiajun (comment out next 3 lines when not in use)
MainFolder = 'E:\文档\Beamline\Data\Results_Asamples\Results_Asamples_Renamed2';
ResultFolder = 'E:\文档\Beamline\Data\CNCUnderMagneticField\Results\';
SaveFolderName = 'CNCMagneticFild8\';

% Folder names for Andrew (comment out next 3 lines when not in use)
% MainFolder = 'C:\Users\anbae\Documents\Research\Isotrpoic_Result_CNC';
% ResultFolder = 'C:\Users\anbae\Documents\Research\CNCUnderMagneticFieldResults\';
% SaveFolderName = 'CNCMagneticField8\';
SaveFolderNameDose = 'Dose Study\'; % folder inside SaveFolderName for the dose study

SaveFolder = [ResultFolder,SaveFolderName];
SaveFolderDose = [ResultFolder,SaveFolderName, SaveFolderNameDose]; % save path for the dose study


CelluloseType = 'CNC';
Solvent = 'PG';
PlotCell{1} = [9, 12, 17, 22]; % folders numbers for concentration variation [9, 12, 17, 22, 32]
PlotCell{2} = [23, 24, 25]; % folder numbers for MagneticfieldExposureTime Variation
PlotCell{3} = [35, 3, 6]; % folder numbers CaCl2 variation

% CaseNum = 1;
% Index = find(contains(C,'bla'));;


% Creat folder to store results
if ~exist(SaveFolder, 'dir')
    mkdir(SaveFolder)
end

% Creat folder to store results for Dose study folder
if ~exist(SaveFolderDose, 'dir')
    mkdir(SaveFolderDose)
end


% Extract key words
NameList = dir(MainFolder);
ExtractIndex = 1;
for NameIndex = 3:length(NameList)
    SampleFolderList{ExtractIndex} = NameList(NameIndex).name; % Names of all samples
    KeywordsCell = strsplit( NameList(NameIndex).name,{CelluloseType,Solvent});
    ExpNumArray{ExtractIndex} = KeywordsCell{1}(1:end-1); % Experiment munber
    VolumeFractionArray{ExtractIndex} = KeywordsCell{2}(1:end-1); % Volum fraction value
    ExpConditionArray{ExtractIndex} = KeywordsCell{3}(2:end); % Experiment condition
    ExtractIndex = ExtractIndex + 1; % Index
end

    

    % Sub1. Loop for different Experiments
    for ExpNumIndex = 1:length(SampleFolderList)
        SubfolderName1 = SampleFolderList{ExpNumIndex};
        
        RunIndex = 0; % I think this is how many doses are going to be analyzed
        % Sub2. Loop for Att
        for AttIndex = 1:length(AttArray)
            Att = AttArray{AttIndex};
            SubfolderName2 = [SubfolderName1,'_Att',num2str(Att)];
            
            % Sub3. Loop for Eps
            for EspIndex = 1:length(EpsArray)
                Esp = EpsArray{EspIndex};
                
                
                % Dose filter, filtering out doses that do not have good data
                clear Dose
                Dose = str2num(Att).*str2num(Esp);
                
                if Dose >= DoseArrayTotal(DoseLolim) &&  Dose <= DoseArrayTotal(DoseUplim) %&& str2num(Att)<= AttUpBond
                   
                
                    RunIndex = RunIndex + 1;
                    SubfolderName3 = [SubfolderName2,'_Eps',num2str(Esp),' ms'];


                    SaveNameCell{RunIndex,ExpNumIndex} = [SubfolderName3,'.mat'];
                    MatSavePath = [MainFolder,'\',SubfolderName1,'\',SaveNameCell{RunIndex,ExpNumIndex}];

                    % Record dose number
                    DoseMatrix(RunIndex,ExpNumIndex) = Dose;

                    % Sub4. Open 2 uid folders
                    UidFolders = dir([MainFolder,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\*uid=*']);
                    for UidIndex = 1:length(UidFolders)
                        Uid{UidIndex} = UidFolders(UidIndex).name;
                        Csvfiles = dir([MainFolder,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\',UidFolders(UidIndex).name]);
                        CsvSize = length(Csvfiles);
                        for CsvIndex = 3:CsvSize
                            CsvName = Csvfiles(CsvIndex).name;
                            switch CsvName(end-5:end-4)
                                case 'g2'

                                    g2CSVDataCell{UidIndex} = readmatrix([MainFolder,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\',UidFolders(UidIndex).name,'\',CsvName]);

                                case '2b'
                                    g2bCSVDataCell{UidIndex} = readmatrix([MainFolder,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\',UidFolders(UidIndex).name,'\',CsvName]);

                                otherwise
                                    continue
                        end
                    end
                end
                else
                    continue
                end
                
                % Store data in mat files
                
                save(MatSavePath,'g2CSVDataCell','g2bCSVDataCell');
                
            end
        end
        
    end
    
    clear UidIndex g2CSVDataCell g2bCSVDataCell DArray G2LengendArray
    
    %% Data Processing and Graphing
    
   
    for CaseNum = 1:3 % determiens what variable (concentration, time, or CaCl2) is changed
        PlotArray = PlotCell{CaseNum}; % the specific folders for this case
        
        switch CaseNum
            case 1
                ExpName = 'Concentration Variation';
                Xarray = [1 2 4 6]; %[1 2 4 6 8]
                XUnit = ' (wt%)';
            case 2
                ExpName = 'MagneticfieldExposureTime Variation';
                Xarray = [0 3 28]; % 
                XUnit = 'Time(h)';
            case 3
                ExpName = 'CaCl2 Variation';
                Xarray = [0.5 1 2];
                XUnit = 'CaCl_2 Unit';   
        end % end of the switch statement
        
        for qNum = 1:5 % go through all of the q values (max of qNum is 5)
            Beta = BetaArray(qNum);
            for PlotIndex = 1:length(Xarray) % this determines what folders to go into for this case 
                ColorVec = 0.8*jet(length(Xarray)); % to be used later in plotting C vs tau
                
                ExpNumIndex = PlotArray(PlotIndex); 
                SavePath = [MainFolder,'\',SampleFolderList{ExpNumIndex}];
                
                % Segment size is the number of different atten and exposure conditions
                % Exp size is the number of different main folders
                [SegmentSize,ExpSize] = size(SaveNameCell);  
                DoseArray = DoseMatrix(:,PlotIndex);
                [DoseArray,DoseIndex] = sort(DoseArray); % Sort index by dose
                for SacnIndex = 1:RunIndex 
                    ScanNameCell{SacnIndex} = SaveNameCell{DoseIndex(SacnIndex),ExpNumIndex}; % creates new name cell with only certain cell
                end
                % Read mat file in folders
                for SegmentIndex = 1:SegmentSize % segment size is the number of different conditions
                    load([SavePath,'\',ScanNameCell{SegmentIndex}]);

                    for UidIndex = 1:length(g2CSVDataCell) % Uid is the experiment number
                        g2CSVData = g2CSVDataCell{UidIndex};
                        g2bCSVData = g2bCSVDataCell{UidIndex};
                        qValue = g2CSVData(1,qNum+2);

                        % calculate all segments data
                        g2Temp = g2CSVData(2:end,qNum+2);
                        g2bTemp = g2bCSVData(2:end,qNum+2);
                        TauTemp = g2CSVData(2:end,2);

                        % Store all segments data
                        TauCell{SegmentIndex,UidIndex} = TauTemp; % this is a 6x2 cell
                        G2Cell{SegmentIndex,UidIndex} = g2Temp;
                        G2bCell{SegmentIndex,UidIndex} = g2bTemp;

                    end
                end
            
                %% Equilibrum Check
                for SgIndex = 1:SegmentSize
                    for UidIndex = 1:RpNum
                        SgTau = TauCell{SgIndex,UidIndex};
                        SeG2 = G2Cell{SgIndex,UidIndex};
                        SgG2b = G2bCell{SgIndex,UidIndex};
                        EqResult = EquilChk(SgTau,SeG2,SgG2b,EqDiff,PassRate); % Pass 1, not pass 0
                        ChkMat(SgIndex,UidIndex) = EqResult;
                    end
                end
                %% Repetition Check
                RpIndex = 0;
                for SgIndex = 1:SegmentSize
                    ChkSgArray(SgIndex) = ChkMat(SgIndex,1)*ChkMat(SgIndex,2);
                    if ChkSgArray(SgIndex)
                        RpIndex = RpIndex + 1;
                        RpResult = RptChk(TauCell{SgIndex,2},G2Cell{SgIndex,1},G2bCell{SgIndex,2},RpDiff);
                        TauCellRp{RpIndex}= RpResult(:,1);
                        G2CellRp{RpIndex} = RpResult(:,2);
                        %DoseArray(RpIndex) = DoseRecord_Array(SgIndex);
                    else
                        continue
                    end
                end
%             %% Plot G2 vs. TAU (Dose Study)
%                 figure;
%                 ColorVec = 0.8*jet(3);
%                 
%                 for PlotIndex = 1:RpIndex
%                     G2LineArray(PlotIndex) = line(TauCellRp{PlotIndex},G2CellRp{PlotIndex},...
%                         'LineStyle','none','Marker','o','LineWidth',2,'Color',ColorVec(PlotIndex,:));
%                     G2LengendArray{PlotIndex} = ['Dose = ', num2str(DoseArray(PlotIndex))];
% 
%                 end
%                 G2Legend = legend(G2LineArray,G2LengendArray,'Interpreter', 'none');
%                 set(gca,'XScale','log','XTick',10.^(linspace(-3,2,6)))
%                 xlabel('\tau [s]')
%                 xlim([10^-3 10^1])
%                 ylabel('g_{2}')
%                 ylim([0.9 1.3])
%                 title('g2 vs \tau');
%                 box off
% 
%                 FigureTitle = ['Dose_Study_',SampleFolderList{ExpNumIndex},'_qNum=', char(string(qNum))];
%                 SaveTitle = [FigureTitle,'.png'];
%                 SavePath = strcat(SaveFolderDose,SaveTitle);
%                 saveas(gcf,SavePath);
                
                %% Data Assembly
                % Combines the uids
                
                AssembleData = SgPackaging(TauCellRp,G2CellRp);
                TAU = AssembleData(:,1);
                G2 = AssembleData(:,2);
                if CaseNum == 1 && PlotIndex == 4 % this gets rid of the 1.8989 dose for the 6% weight concentration because the data was not connected
                    TAU = TauCellRp{1, 1};
                    G2 = G2CellRp{1, 1};
                end
                
                
                % Get rid of zero tau
                ZeroIndex = TAU == 0; %creates a logical vector
                TAU(ZeroIndex) = [];
                G2(ZeroIndex) = [];
                
                % Define G2_Inf 
                BLDataY = G2(end-9:end);
                G2_Inf = mean(BLDataY); %takes the mean of the last 10 data points and defines that as G2 inf

                % Beta Correction
                C = (G2-G2_Inf)/(Beta - G2_Inf + 1); %square of G1
                
               
                [TAU,SortIndex] = sort(TAU); %sorts Tau from least to greatest
                C = C(SortIndex); %sorts g1 according to tau
               
                %% Data Fitting
                FitResult = Str_and_Sin_exp(TAU, C); % fits the data with streched and single exponential
                % Str_and_Sin_exp was renamed from something like StrExp because both streched and single exponential fittings are done
                
                CFitKWW = FitResult{1};
                CoefArrayKWW = FitResult{2};
                GammaKWW(PlotIndex) = CoefArrayKWW(1);
                littleGamma(PlotIndex) = CoefArrayKWW(2); % records little Gamma
        

                CFitSingle = FitResult{4};
                CoefArraySingle = FitResult{5};
                GammaSingle(PlotIndex) = CoefArraySingle(1);

                % Data to be used in the plots
                qArray(PlotIndex) = qValue;
                PlotTAU{PlotIndex} = TAU;
                PlotG2{PlotIndex} = C;
                PlotG2FitKWW{PlotIndex} = CFitKWW;
                PlotG2FitSingle{PlotIndex} = CFitSingle;
             end   
                
        
      %% Plot Fit C vs Tau  
      figure; 
        for PlotIdx = 1:length(Xarray) % data is fitted witht he streched exponential
            KWWLineArray(PlotIdx) = line(PlotTAU{PlotIdx}, PlotG2{PlotIdx}, 'LineStyle','none','Marker','o','LineWidth',2,'Color',ColorVec(PlotIdx,:));
            KWWFitLineArray(PlotIdx) = line(PlotTAU{PlotIdx}, PlotG2FitKWW{PlotIdx}, 'LineWidth',2,'Color',ColorVec(PlotIdx,:));
            KWWLegendArray{PlotIdx} = append(num2str(Xarray(PlotIdx)), XUnit); 
        end
        
        xlabel( '\tau [s]');
        ylabel('C');
        G2Legend = legend(KWWLineArray, KWWLegendArray,'Interpreter', 'none');%,'Location','northwest');
        set(gca,'XScale','log') % sets x scale to logarithimic
        
        FigureTitle = ['C vs tau ', ExpName, ' qNum = ', num2str(qNum)];
        title(FigureTitle);
        SaveTitle = ['C vs tau_', ExpName, '_qNum_', num2str(qNum), '.png'];
        SavePath = strcat(SaveFolder,SaveTitle);
        saveas(gcf,SavePath);
        clear KWWLegendArray % clear the legend array so legend stuff from previous case doesn't carry over to next case
        
%     G2Legend = legend(G2LineArray,G2LengendArray,'Interpreter', 'none');
%     set(gca,'XScale','log','XTick',10.^(linspace(-3,2,6)))
%     xlabel('\tau [s]')
%     xlim([10^-3 10^2])
%     ylabel('g_{2}')
%     ylim([0.9 1.2])
%     title('g2 vs \tau');
%     box off
%     
%     FigureTitle = ExpName;
%     
%     SaveTitle = [FigureTitle,'.png'];
%     SavePath = strcat(SaveFolder,SaveTitle);
%     saveas(gcf,SavePath);
%     
%     % Store D Array
%     DCell{CaseNum} = DArray;
    end % ends the q array
end % end of the for loop for the different cases





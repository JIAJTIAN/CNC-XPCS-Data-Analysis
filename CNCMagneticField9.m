%% Magnetic Field Under XPCS 9
% 1. Import Equilibrium test, Repetition test 
% 2. Imported Dose study


clc,clear
close all

%% Read Data
for CaseNum = 1:2
    AttArray = {'0.007','0.036','0.190','1.000'};
    EpsArray = {'1.334','9.994','99.994'};
    DoseArrayTota = [0.0093 0.0480 0.0700 0.2535 0.3598 0.7000 1.3340 1.8989 3.5998 9.9940 18.9989 99.9940];
    DoseNum = 6;
    Doselim = DoseArrayTota(DoseNum);
    Threshold = 0.2;
    RefinedLength = 6;
    SmoothnessThreshold = 0;
    Beta = 0.12;
    PlotIndex = 1;
    % Find main folder
    MainFolfer = 'E:\文档\Beamline\Data\Results_Asamples\Results_Asamples_Renamed2';
    ResultFolder = 'E:\文档\Beamline\Data\CNCUnderMagneticField\Results\';
    SaveFolderName = 'CNCMagneticFild8\';
    SaveFolder = [ResultFolder,SaveFolderName];
    CelluloseType = 'CNC';
    Solvent = 'PG';
    PlotCell{1} = [8 11 14 16 18 21 27 31];
    PlotCell{3} = [8 1 2 0];
    PlotCell{2} = [22 23 25 24];
    % CaseNum = 1;
    % Index = find(contains(C,'bla'));;
    PlotArray = PlotCell{CaseNum};
    
    switch CaseNum
        case 1
            ExpName = 'Concentratioin Variation';
            Xarray = [1 2 3 4 5 6 7 8];
            XUnit = 'Concentration(wt%)';
        case 3
            ExpName = 'CaCl2 Variation';
            Xarray = [1 2 0.5];
            XUnit = 'CaCl_2 Amount';
        case 2
            ExpName = 'MagneticfieldExposureTime Variation';
            Xarray = [0 1 3 28];
            XUnit = 'Time(h)';
    end
    
    % Creat folder to store results
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder)
    end
    
    
    % Extract key words
    NameList = dir(MainFolfer);
    ExtractIndex = 1;
    for NameIndex = 3:length(NameList)
        SampleFilderList{ExtractIndex} = NameList(NameIndex).name; % Names of all samples
        KeywordsCell = strsplit( NameList(NameIndex).name,{CelluloseType,Solvent});
        ExpNumArray{ExtractIndex} = KeywordsCell{1}(1:end-1); % Experiment munber
        VolumeFractionArray{ExtractIndex} = KeywordsCell{2}(1:end-1); % Volum fraction value
        ExpConditionArray{ExtractIndex} = KeywordsCell{3}(2:end); % Experiment condition
        ExtractIndex = ExtractIndex + 1; % Index
    end
    
    
    
    % Sub1. Loop for different Experiments
    for ExpNumIndex = 1:length(SampleFilderList)
        SubfolderName1 = SampleFilderList{ExpNumIndex};
        
        RunIndex = 0;
        % Sub2. Loop for Att
        for AttIndex = 1:length(AttArray)
            Att = AttArray{AttIndex};
            SubfolderName2 = [SubfolderName1,'_Att',num2str(Att)];
            
            % Sub3. Loop for Eps
            for EspIndex = 1:length(EpsArray)
                Esp = EpsArray{EspIndex};
                
                
                % Dose filter
                clear Dose
                Dose = str2num(Att).*str2num(Esp);
                if Dose < Doselim
                    continue
                end
                
                RunIndex = RunIndex + 1;
                SubfolderName3 = [SubfolderName2,'_Eps',num2str(Esp),' ms'];
                
                
                SaveNameCell{RunIndex,ExpNumIndex} = [SubfolderName3,'.mat'];
                MatSavePath = [MainFolfer,'\',SubfolderName1,'\',SaveNameCell{RunIndex,ExpNumIndex}];
                % Check if data storing is already done
%                 if isfile([MatSavePath,'.mat']) ==1 % Mat file exist, continue
%                     continue
%                 end
                
                % Record dose number
                DoseMatrix(RunIndex,ExpNumIndex) = Dose;
                
                % Sub4. Open 2 uid folders
                UidFolders = dir([MainFolfer,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\*uid=*']);
                for UidIndex = 1:length(UidFolders)
                    Uid{UidIndex} = UidFolders(UidIndex).name;
                    Csvfiles = dir([MainFolfer,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\',UidFolders(UidIndex).name]);
                    CsvSize = length(Csvfiles);
                    for CsvIndex = 3:CsvSize
                        CsvName = Csvfiles(CsvIndex).name;
                        switch CsvName(end-5:end-4)
                            case 'g2'
                                
                                g2CSVDataCell{UidIndex} = readmatrix([MainFolfer,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\',UidFolders(UidIndex).name,'\',CsvName]);
                                
                            case '2b'
                                g2bCSVDataCell{UidIndex} = readmatrix([MainFolfer,'\',SubfolderName1,'\',SubfolderName2,'\',SubfolderName3,'\',UidFolders(UidIndex).name,'\',CsvName]);
                                
                            otherwise
                                continue
                        end
                    end
                end
                
                
                % Store data in mat files
                
                save(MatSavePath,'g2CSVDataCell','g2bCSVDataCell');
                
            end
        end
        
    end
    
    clear UidIndex g2CSVDataCell g2bCSVDataCell DArray G2LengendArray
    
    %% Data Processing
    figure;
    for qNum = 3
        for PlotIndex = 1:length(Xarray)
            ExpNumIndex = PlotArray(PlotIndex);
            SavePath = [MainFolfer,'\',SampleFilderList{ExpNumIndex}];
            [SegmentSize,ExpSize] = size(SaveNameCell);
            DoseArray = DoseMatrix(:,PlotIndex);
            [DoseArray,DoseIndex] = sort(DoseArray); % Sort index by dose
            for SacnIndex = 1:RunIndex
                ScanNameCell{SacnIndex} = SaveNameCell{DoseIndex(SacnIndex),ExpNumIndex};
            end
            % Read mat file in folders
            for SegmentIndex = 1:SegmentSize
                load([SavePath,'\',ScanNameCell{SegmentIndex}]);
                
                for UidIndex = 1:length(g2CSVDataCell)
                    g2CSVData = g2CSVDataCell{UidIndex};
                    g2bCSVData = g2bCSVDataCell{UidIndex};
                    qValue = g2CSVData(1,qNum+2);
                    
                    % calculate all segments data
                    g2Temp = g2CSVData(2:end,qNum+2);
                    g2bTemp = g2bCSVData(2:end,qNum+2);
                    TauTemp = g2CSVData(2:end,2);
                    
                    % Store all srgments data
                    TauCell{SegmentIndex,UidIndex} = TauTemp;
                    g2ValueCell{SegmentIndex,UidIndex} = g2Temp;
                    g2bValueCell{SegmentIndex,UidIndex} = g2bTemp;
                    
                end
            end
            
            %% Merge multiple uid data
            [~,ScanSize] = size(g2ValueCell);
            RefineIndex = 0;
            for SegmentIndex = 1:SegmentSize
                g2ValueCat = [];
                Tau_RefinedCat = [];
                for UidIndex = 1:ScanSize
                    
                    %% Refine segment data by getFinedata function and spline fir R-square
                    DataRefined = getFinedata2(TauCell{SegmentIndex,UidIndex},g2ValueCell{SegmentIndex,UidIndex});
                    
                    % Fit model to data.
%                     [SegmnetFitresult, gof] = fit( DataRefined(:,1), DataRefined(:,2), ft );
%                     Rsqure(SegmentIndex,UidIndex) = gof.rsquare;
%                     FitTau{SegmentIndex,UidIndex} = DataRefined(:,1);
%                     FitG2{SegmentIndex,UidIndex} = SegmnetFitresult(DataRefined(:,1));
                    
                    
%                     if length(DataRefined(:,1)) < RefinedLength 
%                         Tau_RefinedRaw{SegmentIndex,UidIndex} = [];
%                         g2Value_RefinedRaw{SegmentIndex,UidIndex} = [];
%                     else
                        RefineIndex = RefineIndex + 1;
                        Tau_RefinedRaw{SegmentIndex,UidIndex} = DataRefined(:,1);
                        g2Value_RefinedRaw{SegmentIndex,UidIndex} = DataRefined(:,2);
                        
%                     end % end if
                   
                end
            end
            
            
            
            %% G2 Data Baseline Refine Process
            [SegmentSizeBase,UidSizeBase] = size(g2Value_RefinedRaw);
            for SegmentIndex = 1:SegmentSizeBase
                for UidIndex = 1:UidSizeBase
                    
                    SegmentxData = Tau_RefinedRaw{SegmentIndex,UidIndex};
                    SegmentyData = g2Value_RefinedRaw{SegmentIndex,UidIndex};
                    
                    % Set up fittype and options.
                    ft = fittype( 'a1*exp(a2*x)+a3*exp(a4*x)+a5', 'independent', 'x', 'dependent', 'y' );
                    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts.Display = 'Off';
                    opts.Lower = [0 -Inf 0 -Inf 1];
                    opts.StartPoint = [0.255095115459269 0.505957051665142 0.699076722656686 0.890903252535798 0.959291425205444];
                    opts.Upper = [Inf 0 Inf 0 Inf];
                    %                 opts.Upper = [Inf 0 Inf 0 2];
                    
                    % Fit model to data.
                    [SegmnetFitresult, gof] = fit( SegmentxData, SegmentyData, ft, opts );
                    CoefCell{SegmentIndex,UidIndex} = coeffvalues(SegmnetFitresult);
                    
                    
                end
            end
            
            
            % Remove g2 baseline
            
            
            for SegmentIndex = 1:SegmentSizeBase
                for UidIndex = 1:UidSizeBase
                    SegmentCoeffArray = CoefCell{SegmentIndex,UidIndex};
                    Baseline = SegmentCoeffArray(5);
                    g2ValueFitRefined{SegmentIndex,UidIndex} = g2Value_RefinedRaw{SegmentIndex,UidIndex} / Baseline;
                    
                end
            end
            
            
            % Merge multiple scans
            for SegmentIndex = 1:SegmentSizeBase
                Tau_RefinedCat = [];
                g2ValueCat = [];
                for UidIndex = 1:UidSizeBase
                    
                    g2ValueCat = cat(1,g2ValueFitRefined{SegmentIndex,UidIndex},g2ValueCat);
                    Tau_RefinedCat = cat(1,Tau_RefinedRaw{SegmentIndex,UidIndex},Tau_RefinedCat);
                end
                g2Value{SegmentIndex} = g2ValueCat;
                Tau_Refined{SegmentIndex} = Tau_RefinedCat;
            end
            
            % Assamble segments data
            
            CatG2 = [];
            CatTau = [];
            CatIndex = 1;
            for SegmentIndex = 1:SegmentSizeBase
                CatTau = cat(1,CatTau,Tau_Refined{SegmentIndex});
                CatG2 = cat(1,CatG2,g2Value{SegmentIndex});
                CatIndex = CatIndex + 1;
                
            end

            [TAU,TAUSortIndex] = sort(CatTau);
            g2 = CatG2(TAUSortIndex);


            
            %% Fit Refined Data with Double Exponetial Decay
            % Set up fittype and options.
            ft = fittype( 'a1*exp(a2*x)+a3*exp(a4*x)+a5', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [0 -Inf 0 -Inf 1];
            opts.StartPoint = [0.255095115459269 0.505957051665142 0.699076722656686 0.890903252535798 0.959291425205444];
            opts.Upper = [Inf 0 Inf 0 Inf];
            %                 opts.Upper = [Inf 0 Inf 0 2];
            
            % Fit model to data.
            [G2Fitresult, gof] = fit( TAU, g2, ft, opts );
            
            G2Fit = feval(G2Fitresult,TAU);
            CoefArray = coeffvalues(G2Fitresult);
            
            %% Plot G2 and G2 sgement fit
            ColorVec = 0.8*jet(length(PlotArray));
            
            
            % Compute D_t and D_r
            % Gamma_0 = q^2 * D_t
            % Gamma_1 = q^2 * D_t + 6*D_r
            Gamma_0 = -CoefArray(2)/CoefArray(5);
            Gamma_1 = -CoefArray(4)/CoefArray(5);
            if Gamma_0 > Gamma_1
                Gamma_temp = Gamma_0;
                Gamma_0 = Gamma_1;
                Gamma_1 = Gamma_temp;
            end
            D_t = Gamma_0/((qValue)^2); % Unit nm^2/s
            D_r = (Gamma_1-Gamma_0)/6; % Unit 1/s
            
            DArray(PlotIndex,:) = [D_t,D_r];
            
            G2LineArray(PlotIndex) = line(TAU,g2,'LineStyle','none','Marker','o','LineWidth',2,'Color',ColorVec(PlotIndex,:));
            G2LengendArray{PlotIndex} = SampleFilderList{ExpNumIndex};
            G2FitLineArray(PlotIndex) = line(TAU,G2Fit,'LineWidth',2,'Color',ColorVec(PlotIndex,:));
            
        end

    end
    G2Legend = legend(G2LineArray,G2LengendArray,'Interpreter', 'none');
    set(gca,'XScale','log','XTick',10.^(linspace(-3,2,6)))
    xlabel('\tau [s]')
    xlim([10^-3 10^2])
    ylabel('g_{2}')
    ylim([0.9 1.2])
    title('g2 vs \tau');
    box off
    
    FigureTitle = ExpName;
    
    SaveTitle = [FigureTitle,'.png'];
    SavePath = strcat(SaveFolder,SaveTitle);
    saveas(gcf,SavePath);
    
    % Store D Array
    DCell{CaseNum} = DArray;
end



%% Plot D_t and D_r

for CaseNum = 1:2
    figure;
    switch CaseNum
        case 1
            ExpName = 'Concentratioin Variation';
            Xarray = [1 2 3 4 5 6 7 8];
            XUnit = 'Concentration(wt%)';
        case 3
            ExpName = 'CaCl2 Variation';
            Xarray = [1 2 0.5];
            XUnit = 'CaCl_2 Amount';
        case 2
            ExpName = 'MagneticfieldExposureTime Variation';
            Xarray = [0 1 3 28];
            XUnit = 'Time(h)';
    end
    

        DArrayPlot = DCell{CaseNum};
%         display(DArrayPlot);
        D_tArray = DArrayPlot(:,1);
        D_rArray = DArrayPlot(:,2);
        
        yyaxis left
        D_tLine = line(Xarray,D_tArray,'LineStyle','--','Marker','o','LineWidth',2,'Color','blue');
        ylabel('D_t')
        
        yyaxis right
        D_rLine = line(Xarray,D_rArray,'LineStyle','--','Marker','o','LineWidth',2,'Color','red');
        legend({'D_tLine', 'D_rLine'})
        % Dlegend = legend([D_tLine,D_rLine],DLengendArray);
        xlabel(XUnit)
        ylabel('D_r')
        title(ExpName)
        FigureTitle = ['DPlot' ExpName];
        SaveTitle = [FigureTitle,'.png'];
        SavePath = strcat(SaveFolder,SaveTitle);
        saveas(gcf,SavePath);
   
end

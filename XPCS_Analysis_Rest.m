%% CNC XPCS Analysis - Resting Time Study
% Distinguish 2 detectors
% Loop Sample Number
% Change q number to search q and angle text
% Combine 2 detector data and use for CONTIN analysis
% KWW fit

clc,clear
close all
%%
MainFolfer = 'D:\CHX-June2021\Results_Renamed';
ResultFolder = 'D:\CHX-June2021\Analysis\';
SaveFolderName = 'Resting Time\';
MatFolder = 'D:\CHX-June2021\MatFile\';
SaveFolder = [ResultFolder,SaveFolderName];

% Read q and angle list 
DeRound = 4;
% q_Angle_List88 = readmatrix('D:\CHX-June2021\q_angle_88.csv');
% q_Angle_List180 = readmatrix('D:\CHX-June2021\q_angle_180.csv');
% qList88 = unique(round(10^DeRound*q_angle_List88(:,1))/10^DeRound);
qList_Low = [0.0009    0.0023    0.0037    0.0051    0.0065    0.0079    0.0093    0.0107    0.0121    0.0135];
% angleList88 = unique(floor(q_angle_List88(:,2)));
% qList180 = unique(round(10^DeRound*q_angle_List180(:,1))/10^DeRound);
qList_High = [ 0.0023    0.0037    0.0051    0.0065    0.0079    0.0093    0.0108    0.0122    0.0136    0.0150...
    0.0164    0.0178    0.0192    0.0206   0.0220];
% angleList180 = unique(floor(q_angle_List180(:,2)));
AngleList = [  -173  -158  -143  -128  -113   -98   -83   -68   -53   -38   -23    -8];
qNum_High = 4;% 1 to 9
qNum_Low = qNum_High + 1;
AngleNum = 1;

Angle = AngleList(AngleNum);
q_Low = qList_Low(qNum_Low);
q_High = qList_High(qNum_High);

% Creat folder to store results
if ~exist(SaveFolder, 'dir')
    mkdir(SaveFolder)
end

if ~exist(MatFolder, 'dir')
    mkdir(MatFolder)
end


SampleNameList = dir(MainFolfer);
ExtractIdx = 1;
for NameIndex = 3:length(SampleNameList)
    SampleFilderList{ExtractIdx} = SampleNameList(NameIndex).name;
    ExtractIdx = ExtractIdx + 1;
end
%% Data Processing

PlotArray = [14 40 42 46];

ColorVec = 0.8*jet(length(PlotArray));

for ArrayIndex = 1:length(PlotArray)
    ExpNumIndex = PlotArray(ArrayIndex);
    
        MatDataFolder = [MatFolder,SampleFilderList{ExpNumIndex}];
        
        % Read mat file in folders
        MatDataList = dir([MatDataFolder,'\*.mat']);
        
        MatDataIdx_Lowq = 1;
        MatDataIdx_Highq = 1;
        for MatDataIdx = 1:length(MatDataList)
            load([MatDataFolder,'\',MatDataList(MatDataIdx).name]);

%             % High bond check
%             HBArray = G2 > HighBond;
%             if sum(HBArray)/length(HBArray)<0.9
%                 continue
%             end
            
            if length(g2CSVData(1,3:end))>90
                
                % Calculate qNum for high q data
                q_Angle_List180 = g2CSVData(1,3:end);
                for qIdx_High = 1:length(q_Angle_List180)
                    TextCell_High = split(q_Angle_List180(qIdx_High),'_');
                    QValue_High = str2num(TextCell_High{1});
                    QValue_High_Round = round(10^DeRound*QValue_High)/10^DeRound;
                    Angle_High = str2num(TextCell_High{2});
                    if abs((QValue_High_Round - q_High)/q_High) < 0.01 &&...
                            abs((Angle_High - Angle)/Angle) < 0.01
                        qListNum_High = qIdx_High;
                        break
                    end
                end
                G2 = cell2mat(g2CSVData(3:end,qListNum_High+2));
                Q_AngleText = cell2mat(g2CSVData(1,qListNum_High+2));
                
                G2CellTemp_Highq{MatDataIdx_Highq} = G2;
                
                % Base line removal for high q data
                G2_Inf = mean(G2(end-10:end));
                
                TAUCellTemp_Highq{MatDataIdx_Highq} = cell2mat(g2CSVData(3:end,2));
                MatDataIdx_Highq = MatDataIdx_Highq + 1;
            else
                % Calculate qNum for low q data
                q_Angle_List88 = g2CSVData(1,3:end);
                for qIdx_Low = 1:length(q_Angle_List88)
                    TextCell_Low = split(q_Angle_List88(qIdx_Low),'_');
                    QValue_Low = str2num(TextCell_Low{1});
                    QValue_Low_Round = round(10^DeRound*QValue_Low)/10^DeRound;
                    Angle_Low = str2num(TextCell_Low{2});
                    if abs((QValue_Low_Round - q_Low)/q_Low) < 0.01 &&...
                            abs((Angle_Low - Angle)/Angle) < 0.01
                        qListNum_Low = qIdx_Low;
                        break
                    end
                end
                G2 = cell2mat(g2CSVData(3:end,qListNum_High+2));
                Q_AngleText = cell2mat(g2CSVData(1,qListNum_Low+2));

                G2CellTemp_Lowq{MatDataIdx_Lowq} = G2;                
                TAUCellTemp_Lowq{MatDataIdx_Lowq} = cell2mat(g2CSVData(3:end,2));
                MatDataIdx_Lowq = MatDataIdx_Lowq + 1;
                
            end
        end
            % Merge multi scans
            
            MergeData_Lowq = SgPackaging(TAUCellTemp_Lowq,G2CellTemp_Lowq);
            TAUCell_Lowq{ArrayIndex} = MergeData_Lowq(:,1);
            G2Cell_Lowq{ArrayIndex} = MergeData_Lowq(:,2);
            
            MergeData_Highq = SgPackaging(TAUCellTemp_Highq,G2CellTemp_Highq);
            TAUCell_Highq{ArrayIndex} = MergeData_Highq(:,1);
            G2Cell_Highq{ArrayIndex} = MergeData_Highq(:,2);
            
            % Merge 2 detector data
            MergeData_Total = SgPackaging({TAUCell_Lowq{ArrayIndex},TAUCell_Highq{ArrayIndex}},{G2Cell_Lowq{ArrayIndex},G2Cell_Highq{ArrayIndex}});
            TAUCell_Total{ArrayIndex} = MergeData_Total(:,1);
            G2Cell_Total{ArrayIndex} = MergeData_Total(:,2);

            % Fit data with KWW model
            FitResult_Lowq = KWWFit(TAUCell_Lowq{ArrayIndex}, G2Cell_Lowq{ArrayIndex});
            G2Fit_Lowq{ArrayIndex} = FitResult_Lowq{1};
            CoefArray_Lowq{ArrayIndex} = FitResult_Lowq{2};
            
            FitResult_Highq = KWWFit(TAUCell_Highq{ArrayIndex}, G2Cell_Highq{ArrayIndex});
            G2Fit_Highq{ArrayIndex} = FitResult_Highq{1};
            CoefArray_Highq{ArrayIndex} = FitResult_Highq{2};
            
            FitResult_Total = DoubleKWWFit(TAUCell_Total{ArrayIndex}, G2Cell_Total{ArrayIndex});
            G2Fit_Total{ArrayIndex} = FitResult_Total{1};
            CoefArray_Total{ArrayIndex} = FitResult_Total{2};
            
            x = TAUCell_Total{ArrayIndex};
            y = G2Cell_Total{ArrayIndex};


            figure;
            LineArray_Total = line(TAUCell_Total{ArrayIndex},G2Cell_Total{ArrayIndex},'LineStyle','none','Marker','o','LineWidth',1,'Color',ColorVec(ArrayIndex,:));
            LineFit_Total = line(TAUCell_Total{ArrayIndex},G2Fit_Total{ArrayIndex},'LineStyle','-','Marker','none','LineWidth',1,'Color','black');
            legend(LineArray_Total,Q_AngleText,'Interpreter','none');
            set(gca,'XScale','log')
            xlabel('\tau [s]')
            xlim([10^-4 10^2])
            ylabel('g_{2}')
            ylim([0.95 1.2])
            
            box off
            FigureTitle = [SampleFilderList{ExpNumIndex}];
            title(FigureTitle,'Interpreter','none');
            SaveTitle = [FigureTitle,'.png'];
            SavePath = strcat(SaveFolder,SaveTitle);
            saveas(gcf,SavePath);
            LegendArray{ArrayIndex} = SampleFilderList{ExpNumIndex};
%            
            
    
end

% Plot g2 curve with concentration variations
figure;
for PlotIdx = 1:length(TAUCell_Total)
    
    LineArray_Total(PlotIdx) = line(TAUCell_Total{PlotIdx},G2Cell_Total{PlotIdx},'LineStyle','none','Marker','o','LineWidth',1,'Color',ColorVec(PlotIdx,:));
    LineFit_Total(PlotIdx) = line(TAUCell_Total{PlotIdx},G2Fit_Total{PlotIdx},'LineStyle','-','Marker','none','LineWidth',1,'Color',ColorVec(PlotIdx,:));
end
legend(LineArray_Total,LegendArray,'Interpreter','none');

set(gca,'XScale','log')
xlabel('\tau [s]')
xlim([10^-4 10^2])
ylabel('g_{2}')
ylim([0.95 1.2])

box off
FigureTitle = [num2str(Angle),'deg_q=',num2str(q_High),'_Rest'];
title(FigureTitle,'Interpreter','none');
SaveTitle = [FigureTitle,'.png'];
SavePath = strcat(SaveFolder,SaveTitle);
saveas(gcf,SavePath);
    






%% this code is to simulate community culture of user-defined species in batch condition %%
function main_comm_batch

clc
clear
close all

currentDepth = 1; % get the supper path of the current path
currPath = fileparts(mfilename('fullpath'));% get current path
cd(currPath);

load('SPC_Co_Bth.mat');
load('SPC_Co_Bad.mat');
SPC{1}=SPC_Co_Bth;
SPC{2}=SPC_Co_Bad;
n_species=size(SPC,2);

SPC_EM=[8 8 8;
        14 6 6];
Coculture_Data=importdata('BB_COCULTURE_data.txt');    %%????????????????glucose??maltose??starch??????????????????????
Exp_co_data=Coculture_Data.data;
Exp_co_mets=Coculture_Data.textdata;

biom_ini=[0.04 0.02];
Time=Exp_co_data(:,1);
Data=Exp_co_data(:,2:end);
x0=Exp_co_data(1,3:end);
x0=[biom_ini x0]';
[T Y]=common_com_sim_3sub_unit(SPC,n_species,SPC_EM,x0);
Ymodel=Y(:,1:10);
Y_model_co=Y(:,1:10);
Ymodel(:,2)=Ymodel(:,1)+Ymodel(:,2);
Ymodel=Ymodel(:,2:10);
for mn=2:4
    Ymodel(:,mn)=Ymodel(:,mn)./Ymodel(1,mn).*2;   %%???????????????? 2%, 2g/100mL
    Data(:,mn)=Data(:,mn)./Data(1,mn).*2;
end
bmflag=1;
met_udf=Exp_co_mets(2:10);
% plot_simulation(T,Ymodel,Time,Data,bmflag,met_udf);
Tmodel=T;
Texp=Time;
Yexp=Data;

figure13 = figure;
% Create axes
axes1 = axes('Parent',figure13,'YTick',[0 2 4]);
hold(axes1,'on');

% Activate the left side of the axes
yyaxis(axes1,'left');
% Create multiple lines using matrix input to plot
set(gcf,'unit','centimeters','position',[17 27 25 7]);
plot1 = plot(Tmodel,Ymodel(:,2:4),'Linestyle','-','LineWidth',6.5,'Parent',axes1);
set(plot1(1),...
    'Color',[0.7137    0.4431    0.5569]);
set(plot1(2),...
    'Color',[0.9686    0.5333    0.5529]);
set(plot1(3),...
    'Color',[0.6314    0.2392    0.2314]);

% Create multiple lines using matrix input to plot
plot2 = plot(Texp,Yexp(:,2:4),'MarkerSize',18,'LineStyle','none','Parent',axes1);
set(plot2(1),...
    'MarkerFaceColor',[0.7137    0.4431    0.5569],...
    'MarkerEdgeColor',[0.7137    0.4431    0.5569],...
    'Marker','>','linewidth',4);
set(plot2(2),...
    'MarkerFaceColor',[0.9686    0.5333    0.5529],...
    'MarkerEdgeColor',[0.9686    0.5333    0.5529],...
    'Marker','^','linewidth',4);
set(plot2(3),...
    'MarkerFaceColor',[0.6314    0.2392    0.2314],...
    'MarkerEdgeColor',[0.6314    0.2392    0.2314],...
    'Marker','v','linewidth',4);


% Create ylabel
ylabel('Biomass(g/L)','FontSize',20,'FontName','arial');

% Set the remaining axes properties
set(axes1,'YColor',[0 0.447 0.741],'YTick',[0 2 4]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-0.5 4.6]);
% Activate the right side of the axes
yyaxis(axes1,'right');
set(gcf,'unit','centimeters','position',[17 27 25 7]);
% Create plot
plot(Tmodel,Ymodel(:,1),'LineWidth',6.5,...
    'Color',[0.3569    0.5020    0.5843]);

% Create plot
plot(Texp,Yexp(:,1),...
    'MarkerEdgeColor',[0.3569    0.5020    0.5843],...
    'MarkerSize',18,...
    'Marker','^',...
    'LineWidth',4,...
    'LineStyle','none');

% Create ylabel
% ylabel('[logCFU/mL]','FontSize',24,'FontName','Helvetica','Rotation',270);

% Set the remaining axes properties
set(axes1,'YColor',[0.85 0.325 0.098],'YTick',[0 1 2 3],'YScale','linear');
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-0.5 3.5]);
% Create xlabel
xlabel('Time [h]','FontSize',24,'FontName','Helvetica');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[-2 27]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'XTick',[0 5 10 15 20 25]);
% Create textbox
annotation(figure13,'textbox',...
    [0.566824644549765 0.806522123893805 0.329383886255924 0.0862831858407079],...
    'String',{'Coculture-Substrates'},...
    'LineStyle','none',...
    'FontSize',17,...
    'FontAngle','normal',...
    'FitBoxToText','off');

figure14 = figure;

% Create axes
axes1 = axes('Parent',figure14,'YTick',[0 20 40]);
hold(axes1,'on');

% Activate the left side of the axes
% Create multiple lines using matrix input to plot
set(gcf,'unit','centimeters','position',[17 18 25 7]);
plot1 = plot(Tmodel,Ymodel(:,5:end),'Linestyle','-','LineWidth',6.5,'Parent',axes1);
set(plot1(1),...
    'Color',[ 0.8627    0.5569    0.7098]);
set(plot1(2),...
    'Color',[0.0706    0.4196    0.6196]);
set(plot1(3),...
    'Color',[0.9725    0.8039    0.5922]);
set(plot1(4),...
    'Color',[0.9569    0.6784    0.6863]);
set(plot1(5),...
    'Color',[0    0.5294    0.3216]);

% Create multiple lines using matrix input to plot
plot2 = plot(Texp,Yexp(:,5:end),'MarkerSize',18,'LineStyle','none','Parent',axes1);
set(plot2(1),...
    'MarkerFaceColor',[ 0.8627    0.5569    0.7098],...
    'MarkerEdgeColor',[ 0.8627    0.5569    0.7098],...
    'Marker','>','linewidth',4);
set(plot2(2),...
    'MarkerFaceColor',[0.0706    0.4196    0.6196],...
    'MarkerEdgeColor',[0.0706    0.4196    0.6196],...
    'Marker','^','linewidth',4);
set(plot2(3),...
    'MarkerFaceColor',[0.9725    0.8039    0.5922],...
    'MarkerEdgeColor',[0.9725    0.8039    0.5922],...
    'Marker','v','linewidth',4);
set(plot2(4),...
    'MarkerFaceColor',[0.9569    0.6784    0.6863],...
    'MarkerEdgeColor',[0.9569    0.6784    0.6863],...
    'Marker','v','linewidth',4);
set(plot2(5),...
    'MarkerFaceColor',[0    0.5294    0.3216],...
    'MarkerEdgeColor',[0    0.5294    0.3216],...
    'Marker','v','linewidth',4);

% Create ylabel
ylabel('Conc (g/L)','FontSize',20,'FontName','arial');

% Set the remaining axes properties
set(axes1,'YColor',[0 0.447 0.741],'YTick',[0 25 50 75]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-12 85]);

% Set the remaining axes properties
% set(axes1,'YColor',[0.85 0.325 0.098],'YTick',[0 1 2 3],'scale','linear');
% Uncomment the following line to preserve the Y-limits of the axes
% Create xlabel
xlabel('Time [h]','FontSize',24,'FontName','Helvetica');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[-2 27]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'XTick',[0 5 10 15 20 25]);
% Create textbox
annotation(figure14,'textbox',...
    [0.566824644549765 0.806522123893805 0.329383886255924 0.0862831858407079],...
    'String',{'Coculture-Products'},...
    'LineStyle','none',...
    'FontSize',17,...
    'FontAngle','normal',...
    'FitBoxToText','off');



figure15 = figure;
% Create axes
axes1 = axes('Parent',figure15,'YTick',[0 1 2 3]);
hold(axes1,'on');

% Activate the left side of the axes
% Create multiple lines using matrix input to plot
set(gcf,'unit','centimeters','position',[17 10 25 7]);
plot1 = plot(Tmodel,Y(:,1:2),'Linestyle','-','LineWidth',6.5,'Parent',axes1);
set(plot1(1),...
    'Color',[0.3412    0.6078    0.7882]);
set(plot1(2),...
    'Color',[0.9255    0.7333    0.3255]);
ylabel('Conc (mM)','FontSize',20,'FontName','arial');

% Set the remaining axes properties
set(axes1,'YColor',[0 0.447 0.741],'YTick',[0 1 2 3]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-0.25 3]);

% Set the remaining axes properties
% set(axes1,'YColor',[0.85 0.325 0.098],'YTick',[0 1 2 3],'scale','linear');
% Uncomment the following line to preserve the Y-limits of the axes
% Create xlabel
xlabel('Time [h]','FontSize',24,'FontName','Helvetica');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[-2.5 27]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'XTick',[0 5 10 15 20 25]);
% Create textbox
annotation(figure15,'textbox',...
    [0.566824644549765 0.806522123893805 0.329383886255924 0.0862831858407079],...
    'String',{'Coculture-Growth'},...
    'LineStyle','none',...
    'FontSize',17,...
    'FontAngle','normal',...
    'FitBoxToText','off');


figure16=figure;
axes1 = axes('Parent',figure16);

Coculture_Data=importdata('Bth_and_Bado_Co_qPCR.csv');    
Exp_qPCR=Coculture_Data.data;
Exp_qPCR_mets=Coculture_Data.textdata;
Time_qPCR=Exp_qPCR(:,1);
met_qPCR=Exp_qPCR(:,2:end);

a=met_qPCR(:,[1 3]);
b=met_qPCR(:,[2 4]);
hBar = bar(Time_qPCR, a, 1);
for k1 = 1:size(a,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end

hold(axes1,'on');
% hBar = bar(Time_qPCR, a);
set(gcf,'unit','centimeters','position',[17 1 25 7]);

for k1 = 1:size(a,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end

errorbar1 = errorbar(ctr', ydt', b,'linestyle','none','Marker','.');
hold off
set(hBar(1), 'FaceColor', [0.3412    0.6078    0.7882],'BarWidth',1,'linestyle','none')       % Colour First Bar Set
set(hBar(2), 'FaceColor', [0.9255    0.7333    0.3255],'BarWidth',1,'linestyle','none')       % Colour First Bar Set
set(errorbar1(1),'LineWidth',2,...
    'Color',[0.5569    0.5569    0.5569]);
set(errorbar1(2),'LineWidth',2,...
    'Color',[0.5569    0.5569    0.5569]);
set(gca, 'YLim', [-1.5*10^6 1.8*10^8]);
set(axes1,'XTick',[0 5 10 15 20 25]);

xlabel('Time [h]','FontSize',24,'FontName','Helvetica');
box(axes1,'on');
ylabel('qPCR','FontSize',20,'FontName','arial');


% % addpath('D:\Research\BMGF project\Gut_infant\Dynami_Gut_Model\Species_Model\For_Manu\cocluture');
% % SpeciesPath='D:\Research\BMGF project\Gut_infant\Dynami_Gut_Model\Species_Model\For_Manu\cocluture';
% % CurrentPath=pwd;
% % cd(SpeciesPath)
% % ResultFile=[SpeciesPath '\' 'Bth_and_Bado_Co' '_plot_1002' '.xlsx'];
% % %     plot_pure_sim(Tmodel,Ymodel,Texp,Yexp,met);
% % %     plot_pure_simulation(Tmodel,Ymodel,Texp,Yexp,met,name);
% % %     print(gcf,'-dpng',[name '.png']);    %%??????????????????
% %     Single_Simulationi_Result.Model_result=num2cell([T,Ymodel]);
% %     Single_Simulationi_Result.Exp_result=num2cell([Time,Data]);
% % %     Single_Simulationi_Result{i}.growth=num2cell([Tmodel,miu]);
% % %     xlswrite(ResultFile, Single_Simulationi_Result{i}.Model_result, 'Sheet1');
% %     col_header=['Time', met_udf];    %% first row of the file
% %     col_header2={'Time','Model_co_Bth_OD', 'Model_co_Bado_OD'};
% %     col_header3=Exp_qPCR_mets;
% %     Model_output=[col_header;Single_Simulationi_Result.Model_result];
% %     Exp_output=[col_header;Single_Simulationi_Result.Exp_result];
% %     Model_OD_output=[col_header2;num2cell([T, Y_model_co(:,1:2)])];
% %     Exp_qPCR_output=[col_header3;num2cell(Exp_qPCR)];
% % % %     miu_output=[col_header2;Single_Simulationi_Result{i}.growth];
% % %     xlswrite(ResultFile,Model_output , 'Model_data');   %%????????????????????????????
% % %     xlswrite(ResultFile,Exp_output , 'Exp_data');
% % %     xlswrite(ResultFile,Model_OD_output , 'Model_Bth_Bado_OD');
% % %     xlswrite(ResultFile,Exp_qPCR_output , 'Exp_qPCR');
% % 
% % %%  
    

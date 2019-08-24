%% this code is to simulate community culture of user-defined species in batch condition 
function main_pure_batch_mulsub
clc
clear
close all

currentDepth = 1; % get the supper path of the current path
currPath = fileparts(mfilename('fullpath'));% get current path
cd(currPath);

CurrentPath=pwd;

SPC_EM=[14;6;6];
load('SPC_Co_Bad.mat');
[Tmodel,Ymodel,Texp,Yexp,metname]=common_pure_sim_mulsub(SPC_Co_Bad,SPC_EM,1);   


figure11 = figure;
% Create axes
axes1 = axes('Parent',figure11,'YTick',[0 2 4]);
hold(axes1,'on');

% Activate the left side of the axes
yyaxis(axes1,'left');
% Create multiple lines using matrix input to plot
set(gcf,'unit','centimeters','position',[17 26 25 8]);
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
set(gcf,'unit','centimeters','position',[17 26 25 8]);
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
annotation(figure11,'textbox',...
    [0.566824644549765 0.806522123893805 0.329383886255924 0.0862831858407079],...
    'String',{'B.adolescentis-Substrates'},...
    'LineStyle','none',...
    'FontSize',17,...
    'FontAngle','normal',...
    'FitBoxToText','off');

figure12 = figure;

% Create axes
axes1 = axes('Parent',figure12,'YTick',[0 20 40]);
hold(axes1,'on');

% Activate the left side of the axes
% Create multiple lines using matrix input to plot
set(gcf,'unit','centimeters','position',[17 17 25 8]);
plot1 = plot(Tmodel,Ymodel(:,5:end),'Linestyle','-','LineWidth',6.5,'Parent',axes1);
set(plot1(1),...
    'Color',[ 0.0706    0.4196    0.6196]);
set(plot1(2),...
    'Color',[ 0.9569    0.6784    0.6863]);
set(plot1(3),...
    'Color',[0    0.5294    0.3216]);

% Create multiple lines using matrix input to plot
plot2 = plot(Texp,Yexp(:,5:end),'MarkerSize',18,'LineStyle','none','Parent',axes1);
set(plot2(1),...
    'MarkerFaceColor',[ 0.0706    0.4196    0.6196],...
    'MarkerEdgeColor',[ 0.0706    0.4196    0.6196],...
    'Marker','<','linewidth',4);
set(plot2(2),...
    'MarkerFaceColor',[ 0.9569    0.6784    0.6863],...
    'MarkerEdgeColor',[ 0.9569    0.6784    0.6863],...
    'Marker','^','linewidth',4);
set(plot2(3),...
    'MarkerFaceColor',[0    0.5294    0.3216],...
    'MarkerEdgeColor',[0    0.5294    0.3216],...
    'Marker','v','linewidth',4);

% Create ylabel
ylabel('Conc (mM)','FontSize',20,'FontName','arial');

% Set the remaining axes properties
set(axes1,'YColor',[0 0.447 0.741],'YTick',[0 25 50 75 100]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-12 112]);

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
annotation(figure12,'textbox',...
    [0.566824644549765 0.806522123893805 0.329383886255924 0.0862831858407079],...
    'String',{'B.adolescentis-Products'},...
    'LineStyle','none',...
    'FontSize',17,...
    'FontAngle','normal',...
    'FitBoxToText','off');





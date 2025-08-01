%% plot the dynamic with therapy, read "Vtherapy.txt"
%% ========================================================================
Num = 1000;
str = "TherapyNoDeath6/";
set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',2);
%% read simulation data of more patients
  c = cell(Num,1);
OS_simulation = zeros(Num,1);
  tcr    = zeros(Num,1);
for i = 1:Num
    c{i} = load(strcat(str,"VTherapy-",num2str(i-1),".txt"));
    OS_simulation(i) = c{i}(end,1);
    R(i) = c{i}(end,10);
    OS_therapy(i) = c{i}(end,12);
        tcr(i) = max(c{i}(:,14));
end
%%     
  T = c{1}(:,1)/365;
MH1 = zeros(Num,1);
MH2 = zeros(Num,1);
MH3 = zeros(Num,1);
MB1 = zeros(Num,1);
MB2 = zeros(Num,1);
ML1 = zeros(Num,1);
ML2 = zeros(Num,1);
MR  = zeros(Num,1);
MD  = zeros(Num,1);
for i = 1:length(T)
    for j = 1:Num
        MH1(j) = c{j}(i,2);
        MH2(j) = c{j}(i,3);
        MH3(j) = c{j}(i,4);
        MB1(j) = c{j}(i,5);
        MB2(j) = c{j}(i,6);
        MB3(j) = c{j}(i,7);
        ML1(j) = c{j}(i,8);
        ML2(j) = c{j}(i,9);
         MR(j) = c{j}(i,10);
         MD(j) = c{j}(i,12);
    end
    TH1(:,i) = MH1;
    TH2(:,i) = MH2;   
    TH3(:,i) = MH3;
    TB1(:,i) = MB1;    
    TB2(:,i) = MB2;
    TB3(:,i) = MB3;
    TL1(:,i) = ML1;   
    TL2(:,i) = ML2;
    TR(:,i)  = MR;
    TD(:,i)  = MD;
end
%% ------------------------------------------------------------------------
Para = load('vpatients.txt');

Ilive   = zeros(Num,1);
Ideath1 = zeros(Num,1);
Ideath  = zeros(Num,1);

Ilive   = (R<0.01 & tcr'>0);
Ideath1 = (0.80 < R & R < 0.996);
Ideath  = (R>0.80 & tcr'>0);
%% choose the index of chemotherapy beginning
Dindex = zeros(Num,length(T));
index = 1:length(T);
for i = 1:sum(Ilive)
    Dindex(i,:) = (TD(i,:)==0);
    idx = index((TD(i,:)==0));
    DindexLive(i) = idx(end); %  save the index of beginning
end
for i = 1:sum(Ideath)
    Dindex(i,:) = (TD(i,:)==0);
    idx = index((TD(i,:)==0));
    DindexDeath(i) = idx(end); %  save the index of beginning
end
%% boxchart with Remission patients
month = {'0',...
    '','','','','','','','','','','','1',...
    '','','','','','','','','','','','5',...
    '','','','','','','','','','','','10',...
    '','','','','','','','','','','','15',...
    '','','','','','','','','','','','20',...
    };

month2 = {'0',...
    '','','','','','','','','','','','3',...
    '','','','','','','','','','','','15',...
    '','','','','','','','','','','','30',...
    '','','','','','','','','','','','45',...
    '','','','','','','','','','','','60',...
    };

month3 = {'0',...
    '','','','','','','','','','','','2',...
    '','','','','','','','','','','','4',...
    '','','','','','','','','','','','6',...
    '','','','','','','','','','','','8',...
    '','','','','','','','','','','','10',...
    };
%% compared Remission and Released
hCom = figure(7);
k = 1.0;
subplot(2,3,1)
boxchart(TB2(Ilive,1:100/2:floor(end/6)),'MarkerStyle','none','WhiskerLineStyle','none')
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month3,'XTickLabelRotation',0)
xlabel('Time (Month)')
ylabel('B_2')
ylim([0,7e8])
text(-10,7e8,'A','fontweight','bold','FontSize',30);

subplot(2,3,2)
boxchart(TR(Ilive,1:100/2:floor(end/6)),'MarkerStyle','none','WhiskerLineStyle','none')  
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month3,'XTickLabelRotation',0)
xlabel('Time (Month)')
ylabel('r_{LC}')    
ylim([0,0.1])
text(-10,0.1,'B','fontweight','bold','FontSize',30);

subplot(2,3,4)
boxchart(TB2(Ideath,1:100/2:floor(end/6)),'MarkerStyle','none','WhiskerLineStyle','none')  
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month3,'XTickLabelRotation',0)
xlabel('Time (Month)')
ylabel('B_2')
ylim([0,7e8])
text(-10,7e8,'D','fontweight','bold','FontSize',30);


subplot(2,3,5)
boxchart(TR(Ideath,1:100/2:floor(end/6)),'MarkerStyle','none','WhiskerLineStyle','none')  
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month3,'XTickLabelRotation',0)
xlabel('Time (Month)')
ylabel('r_{LC}')
text(-10,1,'E','fontweight','bold','FontSize',30);


TB2D = TB2(Ideath,:);
TRD  = TR(Ideath,:);
subplot(2,3,6)
arrowPlot(TB2D(1,:),TRD(1,:),'number',3,'LineWidth',2,'scale',2);
hold on
plot(TB2D(1,end),TRD(1,end),'k.','MarkerSize',30);
plot(126150000,0.0666737,'r.','MarkerSize',30);
hold off
xlim([min(TB2(8,:)),7e8])
set(gca,'box','on'); 
set(gca,'linewidth',k);
ylim([0,1])
xlabel('B_2')
ylabel('L_{LC}')
text(-1e8,1,'F','fontweight','bold','FontSize',30);

TB2L = TB2(Ilive,:);
TRL  = TR(Ilive,:);
subplot(2,3,3)
arrowPlot(TB2L(1,:),TRL(1,:),'number',10,'LineWidth',2,'scale',2);
hold on
plot(TB2L(1,end),TRL(1,end),'k.','MarkerSize',30);
plot(172853000,0.0617496,'r.','MarkerSize',30);
hold off
xlim([0,6.5e8]);
set(gca,'box','on'); 
set(gca,'linewidth',k);
xlabel('B_2')
ylabel('L_{LC}')
ylim([0,0.1])
text(-1.5e8,0.1,'C','fontweight','bold','FontSize',30);

%% save the figure as pdf
% set(hCom,'Units','centimeters');
%     screenposition = get(hCom,'Position');
%     set(hCom,'PaperPosition',[0 0 1.2*screenposition(3) 0.9*screenposition(4)],...
%             'PaperSize',       [1.2*screenposition(3) 0.9*screenposition(4)]);
% 
% print(hCom,'-dpdf','Fig6.pdf');

%%
Ilive   = (R<0.01);
Ideath  = (R>0.80);
h1 = figure;
set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',2);
xlabels = {'\beta_{L_1}', '\beta_{L_2}', 'd_{L_2}'};
letters = {'A', 'B', 'C'};

for i = 1:3
    subplot(2,3,i);  % 一行三列图的第i幅
    
    % 提取第 i 个参数
    param = Para(:, i);
    
    % No-Remission 组
    death_param = param(Ideath == 1);
    [f_death, x_death] = ksdensity(death_param);
    fill(x_death, f_death, [0.8500, 0.3250, 0.0980], ...
        'FaceAlpha', 0.3, 'EdgeColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 2); hold on;
    
    % Remission 组
    live_param = param(Ilive == 1);
    [f_live, x_live] = ksdensity(live_param);
    fill(x_live, f_live, [0.3010, 0.7450, 0.9330], ...
        'FaceAlpha', 0.3, 'EdgeColor', [0.3010, 0.7450, 0.9330], 'LineWidth', 2);

    % 设置标签和图例
    xlabel(['$', xlabels{i}, '$'], 'Interpreter', 'latex', 'FontSize', 24);
    ylabel('Distribution');
    legend('No-Remission', 'Remission');

    % 添加子图角标 A/B/C
    text(-0.20, 1.0, letters{i}, 'Units', 'normalized', ...
         'FontWeight', 'bold', 'FontSize', 30);
end
subplot(2,3,1)
xlim([0,0.9])
subplot(2,3,2)
xlim([0,0.9])
subplot(2,3,3)
xlim([0,0.6])

% -------------------------------------------------------------------------
% set(h1,'Units','centimeters');
%     screenposition = get(h1,'Position');
%     set(h1,'PaperPosition',[0 0 1.2*screenposition(3) 0.9*screenposition(4)],...
%            'PaperSize',[1.2*screenposition(3) 0.9*screenposition(4)]);
% 
% print(h1,'-dpdf','Fig7.pdf');


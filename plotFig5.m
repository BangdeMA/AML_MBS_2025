%% plot the dynamic of more patient, read "Vpatient-Num.txt"
%% ========================================================================
Num = 1000;
str = "vpatient/";
set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',2);

%% read simulation data of more patients
  c = cell(Num,1);
OS_simulation = zeros(Num,1);
for i = 1:Num
    c{i} = load(strcat(str,"Vpatient-",num2str(i-1),".txt"));
    r(i) = c{i}(end,10);
end
%% read clinical survival data and simulation survival data
        clinical    = csvread('/Users/mabangde/Downloads/科研汇报/Project_AML_Ma_Lai/code/TCGA_os.csv',1,0);
   OS_clinical(:)   = round(clinical(:)/1);
[fc, xc, flow, fup] = ecdf(OS_clinical,  'Function','survivor','Alpha',0.01,'Bounds','on');
              [f,x] = ecdf(OS_simulation,'Function','survivor','Alpha',0.01,'Bounds','on');    
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
    end
    %% remove noise data
    q = quantile(MH1,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    MH1(MH1<lower_bound|MH1>upper_bound) = NaN;
    TH1(:,i) = MH1;

    q = quantile(MH2,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    MH2(MH2<lower_bound|MH2>upper_bound) = NaN;
    TH2(:,i) = MH2;

    q = quantile(MH3,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    MH3(MH3<lower_bound|MH3>upper_bound) = NaN;    
    TH3(:,i) = MH3;

    q = quantile(MB1,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    MB1(MB1<lower_bound|MB1>upper_bound) = NaN;  
    TB1(:,i) = MB1;

    q = quantile(MB2,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    MB2(MB2<lower_bound|MB2>upper_bound) = NaN;     
    TB2(:,i) = MB2;

    q = quantile(MB3,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    MB3(MB3<lower_bound|MB3>upper_bound) = NaN;
    TB3(:,i) = MB3;

    q = quantile(ML1,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    ML1(ML1<lower_bound|ML1>upper_bound) = NaN;
    TL1(:,i) = ML1;

    q = quantile(ML2,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    ML2(ML2<lower_bound|ML2>upper_bound) = NaN;   
    TL2(:,i) = ML2;

    q = quantile(MR,[0.25,0.50,0.75]);
    IQR = q(3)-q(1);
    lower_bound = q(1)-1.5*IQR;
    upper_bound = q(3)+1.5*IQR;
    MR(MR<lower_bound|MR>upper_bound) = NaN; 
    TR(:,i) = MR;
end

%% boxchart which needs MATLAB 2023A
month = {'0','','','','1','','','','2','','','','3','','','','4','','','','5','','','','6'};
k = 1.0;
h1star = figure(5);
subplot(3,3,1)
boxchart(TH1(:,1:50:floor(end/3)),'MarkerStyle','none')
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('H_1')
text(-4,3e7,'A','fontweight','bold','FontSize',30);

subplot(3,3,2)
boxchart(TH2(:,1:50:floor(end/3)),'MarkerStyle','none')
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('H_2')
text(-4,4e6,'B','fontweight','bold','FontSize',30);

subplot(3,3,3)
boxchart(TH3(:,1:50:floor(end/3)),'MarkerStyle','none')
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('H_3')
text(-4,2.5e8,'C','fontweight','bold','FontSize',30);

subplot(3,3,4)
boxchart(TB1(:,1:50:floor(end/3)),'MarkerStyle','none')
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('B_1')
text(-4,6e4,'D','fontweight','bold','FontSize',30);

subplot(3,3,5)
boxchart(TB2(:,1:50:floor(end/3)),'MarkerStyle','none')
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('B_2')
text(-4,8e8,'E','fontweight','bold','FontSize',30);

subplot(3,3,6)
boxchart(TB3(:,1:50:floor(end/3)),'MarkerStyle','none') 
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('B_3')
text(-4,6e7,'F','fontweight','bold','FontSize',30);

subplot(3,3,7)
boxchart(TL1(:,1:50:floor(end/3)),'MarkerStyle','none')   
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('L_1')
text(-4,6e6,'G','fontweight','bold','FontSize',30);

subplot(3,3,8)
boxchart(TL2(:,1:50:floor(end/3)),'MarkerStyle','none')  
set(gca,'box','on'); 
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('L_2')
text(-4,15e7,'H','fontweight','bold','FontSize',30);

subplot(3,3,9)
boxchart(TR(:,1:50:floor(end/3)),'MarkerStyle','none')  
set(gca,'box','on');  
set(gca,'linewidth',k);
set(gca,'XTickLabel',month,'XTickLabelRotation',0)
xlabel('Time (month)')
ylabel('r_{LC}')    
text(-4,0.8,'I','fontweight','bold','FontSize',30);
%% ========================================================================
% save the figure as pdf
set(h1star,'Units','centimeters');
    screenposition = get(h1star,'Position');
    set(h1star,'PaperPosition',[0 0 1.2*screenposition(3) 0.9*screenposition(4)],...
            'PaperSize',       [1.2*screenposition(3) 0.9*screenposition(4)]);

print(h1star,'-dpdf','Fig3.pdf');
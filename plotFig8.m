%% plot the dynamic with therapy with death, read "Vtherapy.txt"
%% ========================================================================
Num = 1000;
str = "VTherapy/";
set(0,'defaultaxesfontsize',12,'defaultlinelinewidth',2);
%% plot more patient for overall survival =================================   
            c = cell(Num,1);
OS_simulation = zeros(Num,1);
OS_therapy    = zeros(Num,1);
       tcr    = zeros(Num,1);
for i = 1:Num
    c{i} = load(strcat(str,"VTherapy-",num2str(i-1),".txt"));
    OS_simulation(i) = c{i}(end,1);
       OS_therapy(i) = c{i}(end,13);
          R_death(i) = c{i}(end,10);
              tcr(i) = max(c{i}(:,14));
end
% the patient number without therapy
sum(tcr == 0)
%% read clinical survival data and simulationm survival data
        clinical    = csvread('TCGA_os.csv',1,0);
   OS_clinical(:)   = round(clinical(:)/1);
              [f,x] = ecdf(OS_simulation,'Function','survivor','Alpha',0.01,'Bounds','on');
            [f2,x2] = ecdf(OS_therapy,   'Function','survivor','Alpha',0.01,'Bounds','on'); 
%%
% read data
data = readtable('TCGA_clinical_OS.csv');

% save 'DeathTime'
death_times = data.DeathTime;

censoring = isnan(death_times); %  recognize 'NA' by 'isnan()'
death_times(censoring) = max(death_times, [], 'omitnan') + 1;    
censoring = double(~censoring); 

[fc, xc, flow, fup] = ecdf(death_times, 'censoring', ~censoring); 
%%
Para = load("vpatients.txt");
Ilive  = [];
Ideath = [];
k = 0;
j = 0;
for i = 1:Num
       subplot(2,2,2)
       hold on;
       plot(c{i}(:,1)/365,c{i}(:,9),'linewidth',2.0);
       hold off;
       xlabel('Time (months)')
       ylabel('L_2')       
       box on;
   if c{i}(end,10)<0.0001 && tcr(i)>0
       k = k+1;
       Ilive = [Ilive,i];
       subplot(2,2,3)
       hold on;
       plot(c{i}(:,1)/365,c{i}(:,10),'linewidth',2.0);
       hold off;
       aa(k) = c{i}(end,10); 
       xlabel('Time (months)')
       ylabel('r_{LC}')
       ylim([0,1])
       box on;
   else if c{i}(end,10)>0.0001 && tcr(i)>0
       j = j+1;
       Ideath = [Ideath,i];
       subplot(2,2,4)
       hold on;
       plot(c{i}(:,1)/365,c{i}(:,10),'linewidth',2.0);
       hold off;
       bb(j) = c{i}(end,10);
       xlabel('Time (months)')
       ylabel('r_{LC}')
       ylim([0,1])
       box on;
    end
   end
end
%%
set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',2);
k = 1.0;
h16 = figure(16);
subplot(2,2,1)
[f,x] = ecdf(OS_simulation,'Function','survivor','Alpha',0.01,'Bounds','on');
h0 = plot(x/30,f,'k');
hold on;
h1 = stairs(xc/30, 1 - flow,'b:');
h2 = stairs(xc/30, 1 - fc,  'b-');
h3 = stairs(xc/30, 1 - fup, 'b:');
hold off;
xlabel('Time (months)')
ylabel('Survival probability')
xlim([0,30])
ylim([0,1])
legend([h0, h2],'Simulation','Clinical','Box','off')
set(gca,'box','on');
set(gca,'linewidth',k);
text(-6,1,'A','fontweight','bold','FontSize',30);

subplot(2,2,2)
[fl,xl] = ecdf(OS_simulation(Ilive),'Function','survivor','Alpha',0.01,'Bounds','on');
[fd,xd] = ecdf(OS_simulation(Ideath),'Function','survivor','Alpha',0.01,'Bounds','on');
% stairs(xl/30,fl);
plot(xl/30,fl);
hold on
% stairs(xd/30,fd);
plot(xd/30,fd);
hold off
xlim([0,30])
ylim([0,1])
xlabel('Time (months)')
ylabel('Survival probability')
legend('Remission','Relapse','Box','off')
set(gca,'box','on');
set(gca,'linewidth',k);
text(-6,1,'B','fontweight','bold','FontSize',30);

% -------------------------------------------------------------------------
% set(h16,'Units','centimeters');
%     screenposition = get(h16,'Position');
%     set(h16,'PaperPosition',[0 0 1.2*screenposition(3) 0.8*screenposition(4)],...
%            'PaperSize',[1.2*screenposition(3) 0.8*screenposition(4)]);
% 
% print(h16,'-dpdf','FigOS.pdf');

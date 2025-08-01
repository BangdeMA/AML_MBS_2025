%% plot the dynamic with therapy with death, read "Vtherapy.txt"
%% ========================================================================
Num = 1000;
set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',2);
%% read data
OS_simulation = zeros(Num,1);
A = load('CombTherapyOSComb.txt');
B = load('CombTherapyOSCB.txt');
C = load('CombTherapyOSBC.txt');
D = load('CombTherapyOSC.txt');
OS_simulationA  = A(:,2);
OS_simulationB  = B(:,2);
OS_simulationC  = C(:,2);
OS_simulationD  = D(:,2);
%% read data
        clinical    = csvread('/Users/mabangde/Downloads/科研汇报/Project_AML_Ma_Lai/code/TCGA_os.csv',1,0);
   OS_clinical(:)   = round(clinical(:)/1);
[f, x, flow, fup]   = ecdf(OS_clinical,    'Function','survivor','Alpha',0.01,'Bounds','on');
            [fa,xa] = ecdf(OS_simulationA, 'Function','survivor','Alpha',0.01,'Bounds','on');
            [fb,xb] = ecdf(OS_simulationB, 'Function','survivor','Alpha',0.01,'Bounds','on');
            [fc,xc] = ecdf(OS_simulationC, 'Function','survivor','Alpha',0.01,'Bounds','on');
            [fd,xd] = ecdf(OS_simulationD, 'Function','survivor','Alpha',0.01,'Bounds','on');
%% ========================================================================
h5 = figure(5);
k  = 1.0;
% ---- subplot(1) ----
subplot(2,3,1)
plot(xa/30,fa);
hold on
plot(xd/30,fd,'r')
hold off
xlim([0,36])
ylim([0,1])
xlabel('Time (months)')
ylabel('Survival probability')
legend('(low-chemo,high targeted)','Chemotherapy','box','off')
set(gca,'box','on');
set(gca,'linewidth',k);
text(-6,1,'A','fontweight','bold','FontSize',30);

% ---- subplot(2) ----
subplot(2,3,2)
plot(xb/30,fb);
hold on
plot(xd/30,fd,'r')
hold off
xlim([0,36])
ylim([0,1])
xlabel('Time (months)')
ylabel('Survival probability')
legend('high-chemo+high targeted','Chemotherapy','box','off')
set(gca,'box','on');
set(gca,'linewidth',k);
text(-6,1,'B','fontweight','bold','FontSize',30);

% ---- subplot(3) ----
subplot(2,3,3)
plot(xc/30,fc);
hold on
plot(xd/30,fd,'r')
hold off
xlim([0,36])
ylim([0,1])
xlabel('Time (months)')
ylabel('Survival probability')
legend('high targeted+low-chemo','Chemotherapy','box','off')
set(gca,'box','on');
set(gca,'linewidth',k);
text(-6,1,'C','fontweight','bold','FontSize',30);


nColor = 256;
blue_to_red = [linspace(0, 1, nColor)', ... % R: 0 → 1
               zeros(nColor,1), ...
               linspace(1, 0, nColor)'];    % B: 1 → 0

% ---- subplot(4) ----
subplot(2,3,4)
X = A(:,3); Y = A(:,6);
data = [X Y];
density = ksdensity(data, data);
scatter(X, Y, 40, density, 'filled');
hold on
plot(min(A(:,3)),min(A(:,6)),'r.','MarkerSize',25)
hold off
xlim([min(A(:,3)),max(A(:,3))])
ylim([0,0.5])
colormap(blue_to_red);
clim([min(density), max(density)]);
xlabel('B_2'); ylabel('r_{LC}');
set(gca,'box','on','linewidth',1);
text(-0.6e8,0.5,'D','fontweight','bold','FontSize',30);


ax = gca;
cb = colorbar('northoutside');

% ---- subplot(5) ----
subplot(2,3,5)
X = B(:,3); Y = B(:,6); data = [X Y];
density = ksdensity(data, data);
density_size = rescale(density, 20, 100);
scatter(X, Y, 40, density, 'filled');
hold on
plot(min(A(:,3)),min(A(:,6)),'r.','MarkerSize',25)
hold off
xlim([min(B(:,3)),max(B(:,3))])
ylim([0,0.5])
colormap(blue_to_red);
clim([min(density), max(density)]);
xlabel('B_2'); ylabel('r_{LC}');
set(gca,'box','on','linewidth',1);
ax = gca;
cb = colorbar('northoutside');
text(-0.6e8,0.5,'E','fontweight','bold','FontSize',30);

% ---- subplot(6) ----
subplot(2,3,6)
X = C(:,3); Y = C(:,6); data = [X Y];
density = ksdensity(data, data);
density_size = rescale(density, 20, 100);
scatter(X, Y, 40, density, 'filled');
hold on
plot(min(A(:,3)),min(A(:,6)),'r.','MarkerSize',25)
hold off
xlim([min(C(:,3)),max(C(:,3))])
ylim([0,0.5])
colormap(blue_to_red);
clim([min(density), max(density)]);
xlabel('B_2'); ylabel('r_{LC}');
set(gca,'box','on','linewidth',1);
ax = gca;
cb = colorbar('northoutside');
text(-0.6e8,0.5,'F','fontweight','bold','FontSize',30);


%% ========================================================================
set(h5,'Units','centimeters');
    screenposition = get(h5,'Position');
    set(h5,'PaperPosition',[0 0 1.2*screenposition(3) 0.8*screenposition(4)],...
           'PaperSize',[1.2*screenposition(3) 0.8*screenposition(4)]);

print(h5,'-dpdf','Fig10.pdf');
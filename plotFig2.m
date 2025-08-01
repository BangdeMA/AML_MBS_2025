%% plot the dynamic of single patient, read  "patient.txt"
%% ========================================================================
Num = 1000;
set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',3);
 A = load('Patient.txt'); 
 T = A(:,1)/30;
H1 = A(:,2);
H2 = A(:,3);
H3 = A(:,4);
B1 = A(:,5);
B2 = A(:,6);
B3 = A(:,7);
L1 = A(:,8);
L2 = A(:,9);
 R = A(:,10);    
%   lambda  = A(:,11);
% mutation  = A(:,12);
%% simulation data of single patient    
h2 = figure(3);
k  = 1.0;
subplot(3,3,1)
plot(T,H1);
set(gca,'linewidth',k);
xlabel('Time (month)');
xlim([0,15])
ylabel('H_1');
text(-2,3e7,'A','fontweight','bold','FontSize',30);
% text(-2*20,10e7,'A','fontweight','bold','FontSize',30);

subplot(3,3,2)
plot(T,H2)
set(gca,'linewidth',k);
xlabel('Time (month)');
xlim([0,15])
ylabel('H_2')
text(-2,4e6,'B','fontweight','bold','FontSize',30);

subplot(3,3,3) 
plot(T,H3)
set(gca,'linewidth',k);
xlabel('Time (month)');
xlim([0,15])
ylabel('H_3')
text(-2,3e8,'C','fontweight','bold','FontSize',30);

subplot(3,3,4)
plot(T,B1)
set(gca,'linewidth',k);
xlabel('Time (month)');
%     xlim([0,1])
ylabel('B_1')
xlim([0,15])
text(-2,max(B1),'D','fontweight','bold','FontSize',30);

subplot(3,3,5)
plot(T,B2)
set(gca,'linewidth',k);
xlabel('Time (month)');
ylabel('B_2')
xlim([0,15])
text(-2,8e8,'E','fontweight','bold','FontSize',30);

subplot(3,3,6)
plot(T,B3)
set(gca,'linewidth',k);
xlabel('Time (month)');
ylabel('B_3')
xlim([0,15])
text(-2,6e7,'F','fontweight','bold','FontSize',30);

subplot(3,3,7)
plot(T,L1)
set(gca,'linewidth',k);
xlabel('Time (month)');
ylabel('L_1')
xlim([0,15])
text(-2,3e6,'G','fontweight','bold','FontSize',30);

subplot(3,3,8)
plot(T,L2)
set(gca,'linewidth',k);
xlabel('Time (month)');
ylabel('L_2')
xlim([0,15])
text(-2,10.5e7,'H','fontweight','bold','FontSize',30);

subplot(3,3,9)
plot(T,R)
hold on
plot(T,0.2*ones(length(T),1),'--');
hold off
set(gca,'linewidth',k);
ylim([0,1])
xlabel('Time (month)');
xlim([0,15])
ylabel('r_{LC}')  
text(-2,1,'I','fontweight','bold','FontSize',30);

%% save the figure as pdf
set(gcf,'Units','centimeters');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 1.2*screenposition(3) 0.9*screenposition(4)],...
           'PaperSize',[1.2*screenposition(3) 0.9*screenposition(4)]);

print(h2,'-dpdf','Fig2.pdf');

%%
% figure(5)
% subplot(2,2,1)
% plot(T,lambda)
% 
% figure(5)
% subplot(2,2,2)
% plot(T,mutation)
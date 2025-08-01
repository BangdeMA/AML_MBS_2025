%% plot the sensitive analysis result, read  "patient-Num.txt"
%% ========================================================================
Num = 1000;
str = "patient/";
set(0,'defaultaxesfontsize',12,'defaultlinelinewidth',2);

%% read simulation data of more patients
  c = cell(Num,1);
OS_simulation = zeros(Num,1);
for i = 1:Num
    c{i} = load(strcat(str,"patient-",num2str(i-1),".txt"));
    OS_simulation(i) = c{i}(end,1);
end
  YY = load(strcat(str,"parameter.txt"));  
   Y = [YY(:,1:11),YY(:,14:19+4)];
H1SA = zeros(Num,1);
H2SA = zeros(Num,1);
H3SA = zeros(Num,1);
B1SA = zeros(Num,1);
B2SA = zeros(Num,1);
B3SA = zeros(Num,1);
L1SA = zeros(Num,1);
L2SA = zeros(Num,1);
 RSA = zeros(Num,1);
for i=1:Num
    H1SA(i) = c{i}(end,2);
    H2SA(i) = c{i}(end,3);
    H3SA(i) = c{i}(end,4);
    B1SA(i) = c{i}(end,5);
    B2SA(i) = c{i}(end,6);
    B3SA(i) = c{i}(end,7);
    L1SA(i) = c{i}(end,8);
    L2SA(i) = c{i}(end,9);
     RSA(i) = c{i}(end,10);
end

%% scatter plot
figure(6)
for i = 1:size(Y,2)
   subplot(5,5,i);
   scatter(Y(:,i),RSA,'.');
end

%% ------------------------------------------------------------------------
  Num = 1000;
Ilive  = zeros(Num,1);
Ideath = zeros(Num,1);

R = RSA;
Ilive  = (R<0.0001);
Ideath = (R>0.9999 & R<1);

h1 = figure;
set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',2);
xlabels = { '\beta_{H_1}','\beta_{H_2}','\beta_{B_1}','\beta_{B_2}','\beta_{L_1}','\beta_{L_2}',...
'\kappa_{H_1}','\kappa_{H_2}','\kappa_{B_1}','\kappa_{B_2}','\kappa_{L_1}',...
'd_{H_3}','d_{B_1}','d_{B_2}','d_{B_3}','d_{L_1}','d_{L_2}', 'K_{H1B1}','K_{H2B3}','K_{B1L2}','K_{L2B2}'};  

for i = 1:17+4
    subplot(5,5,i);  
    
    param = Y(:, i);
    
    % No-Remission
    death_param = param(Ideath == 1);
    [f_death, x_death] = ksdensity(death_param);
    fill(x_death, f_death, [0.8500, 0.3250, 0.0980], ...
        'FaceAlpha', 0.3, 'EdgeColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 2); hold on;
    
    % Remission
    live_param = param(Ilive == 1);
    [f_live, x_live] = ksdensity(live_param);
    fill(x_live, f_live, [0.3010, 0.7450, 0.9330], ...
        'FaceAlpha', 0.3, 'EdgeColor', [0.3010, 0.7450, 0.9330], 'LineWidth', 2);

    % KS test
    [~, pval] = kstest2(death_param, live_param);

    xlabel(['$', xlabels{i}, '$'], 'Interpreter', 'latex', 'FontSize', 24);
    
    text(0.05, 0.85, ['p = ', num2str(pval, '%.3g')], ...
        'Units', 'normalized', 'FontSize', 14, 'Color', 'k');

    if pval < 0.05
        text(0.05, 0.70, '*', 'Units', 'normalized', ...
            'FontSize', 24, 'Color', 'r', 'FontWeight', 'bold');
    end

end

%% ========================================================================
% save the figure as pdf
set(gcf,'Units','centimeters');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 1.0*screenposition(3) 1.2*screenposition(4)],...
           'PaperSize',[1.0*screenposition(3) 1.2*screenposition(4)]);

print(h1,'-dpdf','Fig3.pdf');

%% ------------------------------------------------------------------------ 
selected_idx = [2, 5, 6, 8, 11, 17];
param_labels = { '\beta_{H_2}', '\beta_{L_1}', '\beta_{L_2}', ...
                 '\kappa_{H_2}', '\kappa_{L_1}', 'd_{L_2}' };

Y_sub = Y(:, selected_idx);

Ilive  = (RSA < 0.0001);
Ideath = (RSA > 0.9999 & RSA<1);

h1  = figure;
set(0, 'defaultaxesfontsize', 24, 'defaultlinelinewidth', 2.0);

n = length(selected_idx);  

for i = 1:n
    for j = 1:n
        subplot(n, n, (i-1)*n + j);

        if i == j
            % KDE distribution
            param = Y_sub(:, i);

            % No-Remission
            dp = param(Ideath == 1);
            [f_d, x_d] = ksdensity(dp);
            fill(x_d, f_d, [0.8500, 0.3250, 0.0980], ...
                'FaceAlpha', 0.3, 'EdgeColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5); hold on;

            % Remission
            lp = param(Ilive == 1);
            [f_l, x_l] = ksdensity(lp);
            fill(x_l, f_l, [0.3010, 0.7450, 0.9330], ...
                'FaceAlpha', 0.3, 'EdgeColor', [0.3010, 0.7450, 0.9330], 'LineWidth', 1.5);

            axis tight; axis square;

            % KS test
            [~, pval] = kstest2(dp, lp);
            
            text(0.05, 0.85, ['p = ', num2str(pval, '%.3g')], ...
                'Units', 'normalized', 'FontSize', 14, 'Color', 'k');
        
            if pval < 0.05
                text(0.05, 0.70, '*', 'Units', 'normalized', ...
                    'FontSize', 24, 'Color', 'r', 'FontWeight', 'bold');
            end

        elseif j < i
            % sub：Remission heatplot
            x = Y_sub(Ilive == 1, j);
            y = Y_sub(Ilive == 1, i);

            xi = linspace(min(x), max(x), 100);
            yi = linspace(min(y), max(y), 100);
            [X, Yg] = meshgrid(xi, yi);
            [f, ~] = ksdensity([x, y], [X(:), Yg(:)]);
            F = reshape(f, size(X));

            imagesc(xi, yi, F);
            axis xy square;
            xlim([min(xi), max(xi)]);
            ylim([min(yi), max(yi)]);

        elseif j > i
            % up：No-Remission heatplot
            x = Y_sub(Ideath == 1, j);
            y = Y_sub(Ideath == 1, i);

            xi = linspace(min(x), max(x), 100);
            yi = linspace(min(y), max(y), 100);
            [X, Yg] = meshgrid(xi, yi);
            [f, ~] = ksdensity([x, y], [X(:), Yg(:)]);
            F = reshape(f, size(X));

            imagesc(xi, yi, F);
            axis xy square;
            xlim([min(xi), max(xi)]);
            ylim([min(yi), max(yi)]);
        end

        xlabel(['$', param_labels{j}, '$'], 'Interpreter', 'latex');
        ylabel(['$', param_labels{i}, '$'], 'Interpreter', 'latex');
        if i==j 
            ylabel('');
        end
    end
end
%% ------------------------------------------------------------------------

legend_center_y = 0.03;
box_w = 0.015;
box_h = 0.015;

y_box = legend_center_y - box_h/2;
y_text = y_box;

% --- Remission legend ---
annotation('rectangle', [0.18, y_box, box_w, box_h], ...
    'FaceColor', [0.3010, 0.7450, 0.9330], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
annotation('textbox', [0.18 + box_w + 0.005, y_text, 0.1, box_h], ...
    'String', 'Favourable', 'FontSize', 18, ...
    'EdgeColor', 'none', 'VerticalAlignment', 'middle');

% --- No-Remission legend ---
annotation('rectangle', [0.34, y_box, box_w, box_h], ...
    'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
annotation('textbox', [0.34 + box_w + 0.005, y_text, 0.15, box_h], ...
    'String', 'Unfavourable', 'FontSize', 18, ...
    'EdgeColor', 'none', 'VerticalAlignment', 'middle');

% --- colorbar ---
cmap = parula(256);
cbar_img = reshape(cmap, [1, size(cmap,1), 3]);

% colorbar location
cbar_x = 0.665;
cbar_y = 0.02;
cbar_w = 0.18;
cbar_h = 0.025;
cbar_center_y = cbar_y + cbar_h / 2;

axes('Position', [cbar_x, cbar_y, cbar_w, cbar_h]);
image(cbar_img);
axis off;

% colorbar location
y_textbox = cbar_center_y - box_h / 2;

annotation('textbox', [0.63, y_box, 0.03, box_h], ...
    'String', 'Low', 'EdgeColor', 'none', ...
    'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

annotation('textbox', [0.85, y_box, 0.05, box_h], ...
    'String', 'High', 'EdgeColor', 'none', ...
    'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%% ------------------------------------------------------------------------
% save the figure as pdf
set(gcf,'Units','centimeters');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 1.0*screenposition(3) 1.2*screenposition(4)],...
           'PaperSize',[1.0*screenposition(3) 1.2*screenposition(4)]);

print(h1,'-dpdf','Fig3.pdf');
%% 
   c = cell(Num,1);
OS_simulation = zeros(Num,1);
for i = 1:Num
    c{i} = load(strcat(str,"patient-",num2str(i-1),".txt"));
    OS_simulation(i) = c{i}(end,1);
end
  YY = load(strcat(str,"parameter.txt"));  
   Y = [YY(:,1:11),YY(:,14:19)];
H1SA = zeros(Num,1);
H2SA = zeros(Num,1);
H3SA = zeros(Num,1);
B1SA = zeros(Num,1);
B2SA = zeros(Num,1);
B3SA = zeros(Num,1);
L1SA = zeros(Num,1);
L2SA = zeros(Num,1);
 RSA = zeros(Num,1);
for i=1:Num
    H1SA(i) = c{i}(floor(end/36),2);
    H2SA(i) = c{i}(floor(end/36),3);
    H3SA(i) = c{i}(floor(end/36),4);
    B1SA(i) = c{i}(floor(end/36),5);
    B2SA(i) = c{i}(floor(end/36),6);
    B3SA(i) = c{i}(floor(end/36),7);
    L1SA(i) = c{i}(floor(end/36),8);
    L2SA(i) = c{i}(floor(end/36),9);
     RSA(i) = c{i}(floor(end/36),10);
end

set(0,'defaultaxesfontsize',24,'defaultlinelinewidth',2);
h1 = figure(4);
subplot(3,1,1)
SA(Y(Ideath,:),H3SA(Ideath,:)) 
text(-0.07, 1.0, 'A', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

subplot(3,1,2)
SA(Y(Ideath,:),L2SA(Ideath,:))
text(-0.07, 1.0, 'B', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

subplot(3,1,3)
SA(Y(Ideath,:),RSA(Ideath,:))
set(gca,'box','on'); 
text(-0.07, 1.0, 'C', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

% save the figure as pdf
% set(gcf,'Units','centimeters');
%     screenposition = get(gcf,'Position');
%     set(gcf,'PaperPosition',[0 0 1.0*screenposition(3) 1.2*screenposition(4)],...
%            'PaperSize',[1.0*screenposition(3) 1.2*screenposition(4)]);
% 
% print(h1,'-dpdf','Fig4.pdf');
%%
% figure(6)
% for i = 1:size(Y,2)
%    subplot(5,5,i);
%    scatter(Y(Ideath,i),RSA(Ideath,:),'.');
% end
% title('r_{LC}')


%%
function SA(P,R)
    % P means parameters, while R means Results
    p = {'$\beta_{H_1}$', '$\beta_{H_2}$', '$\beta_{B_1}$', '$\beta_{B_2}$', '$\beta_{L_1}$', '$\beta_{L_2}$', ...
     '$\kappa_{H_1}$', '$\kappa_{H_2}$', '$\kappa_{B_1}$', '$\kappa_{B_2}$', '$\kappa_{L_1}$', ...
     '$d_{H_3}$', '$d_{B_1}$', '$d_{B_2}$', '$d_{B_3}$', '$d_{L_1}$', '$d_{L_2}$'};

    S = zeros(numel(p),1); 
    pval = zeros(numel(p),1);
    for i = 1:size(P,2)    
       [S(i),pval(i)] = corr(P(:,i),R);  
    end
    I_pval = find(abs(pval)<0.05);
    S = S(I_pval)
    p = p(I_pval)
    [~, I_sort] = sort(S,'descend');
    S_sorted = S(I_sort);         
    p = p(I_sort);
%     p1=bar(S_sorted);
%     set(p1,'FaceColor',[100 100 100]./255);
%     hold on
%     ax = gca;
%     ax.XTickLabel = p;
%     ax.XTickLabel = ax.XTickLabel(I_sort);
%     ax = gca;
%     ax.XTickLabelRotation=0;
    bar(S_sorted,'LineWidth',1.5);
% -------------------------------------------------------------------------
    % % h = my_xticklabels(gca,1:1:length(p),p,'HorizontalAlignment','left','Fontsize',24);
    % xticks(1:1:length(p));
    % xticklabels(p);
    % % set(gca, 'FontSize', 24, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    % set(gca, 'FontSize', 24, 'TickLabelInterpreter', 'latex'); 
    % % axis([0 9 -0.08 0.08])
    % ylabel('PCCs','Fontsize',24)
%--------------------------------------------------------------------------
xticks(1:length(p));
xticklabels(p);

ax = gca; 
ax.FontSize = 24;
ax.XTickLabelRotation = 0;


ax.XAxis.TickLabelInterpreter = 'latex';  
ax.YAxis.TickLabelInterpreter = 'none';   

ylabel('PCCs', 'FontSize', 24); 
end
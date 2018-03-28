function construct_subplots(Y1,pr,savename, method)

% Create figure
figure;

% Create axes
axes1 = subplot(1,2,1);
hold(axes1,'on');

% Create plot
for i = 1:size(Y1,2)
    plot(pr.mspan1,Y1(1:length(pr.mspan1),i),'-s','DisplayName',['s=',num2str(pr.s_span(i))],'LineWidth',2);
end
% Create xlabel
xlabel({'\textbf{Number of measurements} $\mathbf{(m)}$'},...
    'LineWidth',1,...
    'Interpreter','latex');

% Create title
% switch method
%     case 'cosamp'
% title(['\textbf{Relative reconstruction error vs number of measurements; for CoSaMP with} $\mathbf{R=',num2str(pr.R), ',n=',num2str(pr.n),'}$'],...
%     'Interpreter','latex');
%     case 'robust-cosamp'
%     title(['\textbf{Relative reconstruction error vs number of measurements; for robust CoSaMP with} $\mathbf{R=',num2str(pr.R),',n=',num2str(pr.n),'}$'],...
%     'Interpreter','latex');
% end

% Create ylabel
ylabel({'\textbf{Reconstruction Error}; $\mathbf{\frac{||x^*-x||}{||x^*||}}$'},...
    'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',20,'FontName','MS Sans Serif');
grid on;
grid minor;
ylim([0,1]);
%yticks(0:0.1:1.5)

hold(axes1,'off')

axes2 = subplot(1,2,2);
hold(axes2,'on');

% Create plot
for i = 1:size(Y1,2)
    plot(pr.mspan2,Y1(length(pr.mspan1)+1:end,i),'-s','DisplayName',['s=',num2str(pr.s_span(i))],'LineWidth',2);
end
% Create xlabel
xlabel({'\textbf{Number of measurements} $\mathbf{(m)}$'},...
    'LineWidth',1,...
    'Interpreter','latex');

% Create title
% switch method
%     case 'cosamp'
% title(['\textbf{Relative reconstruction error vs number of measurements; for CoSaMP with} $\mathbf{R=',num2str(pr.R), ',n=',num2str(pr.n),'}$'],...
%     'Interpreter','latex');
%     case 'robust-cosamp'
%     title(['\textbf{Relative reconstruction error vs number of measurements; for robust CoSaMP with} $\mathbf{R=',num2str(pr.R),',n=',num2str(pr.n),'}$'],...
%     'Interpreter','latex');
% end




% Create ylabel
%ylabel({'\textbf{Reconstruction Error}; $\mathbf{\frac{||x^*-x||}{||x^*||}}$'},...
%    'Interpreter','latex');

box(axes2,'on');
% Set the remaining axes properties
set(axes2,'FontSize',16);
% Create legend
legend1 = legend(axes2,'show');
set(legend1,'FontSize',20,'FontName','MS Sans Serif');
grid on;
grid minor;
ylim([0,1]);
%yticks(0:0.1:1.5)

hold(axes2,'off');
switch method
    case 'cosamp'
 p = mtit(['\textbf{Relative reconstruction error vs number of measurements; for CoSaMP with} $\mathbf{||x^*||=',num2str(pr.del),'}$'],...
    'Interpreter','latex','FontSize',20);
    case 'robust-cosamp'
  p =  mtit(['\textbf{Relative reconstruction error vs number of measurements; for robust CoSaMP with} $\mathbf{||x^*||=',num2str(pr.del),'}$'],...
    'Interpreter','latex');
end

savefig(['./results/mod_recovery_results/',savename,'_double.fig'])


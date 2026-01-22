%% Define universal parameters

delta_thetas = -90:0.01:90;
c = sqrt(2)/exp(-0.5);

%% Baseline DoG

A = 6;
w = 1/20;
b = 0;

y = (delta_thetas * A * w * c) .* exp(-(w * delta_thetas).^2) + b;

%%% Plot baseline DoG %%%

figure('Color', 'w');

plot(delta_thetas, y,'LineWidth', 1, 'Color', 'k');

xlim([-90, 90]);
ylim([-10, 10]);
ylabel('Bias (°)');
xlabel('\Delta\theta (°)');
xticks([-90:45:90]);
yticks([-10:5:10]);
set(gca,'TickDir','out', 'TickLength', [0.020, 0.020], 'FontSize', 12, 'FontName', 'Helvetica');
box off;
axis square;

hold on;
line([0, 0], [-10, 10], 'LineWidth', 1, 'Color', 'k');
line([-90, 90], [0, 0], 'LineWidth', 1, 'Color', 'k');


%% Calculate DoG with different amplitudes

A = [2 4 6]';
w = 1/20;
b = 0;

y = (delta_thetas .* A .* w .* c) .* exp(-(w * delta_thetas).^2) + b;


%%% Plot DoG %%%

figure('Color', 'w');

plot(delta_thetas, y,'LineWidth', 1);

xlim([-90, 90]);
ylim([-8, 8]);
ylabel('Bias (°)');
xlabel('\Delta\theta (°)');
xticks([-90:45:90]);
yticks(-8:4:8);
set(gca,'TickDir','out', 'TickLength', [0.020, 0.020], 'FontSize', 12, 'FontName', 'Helvetica');
box off;
axis square;

hold on;
line([0, 0], [min(ylim), max(ylim)], 'LineWidth', 1, 'Color', 'k');
line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');

legend({'2°', '4°', '6°'});

%% Calculate DoG with different widths

A = 4;
w = [1/10 1/20 1/30]';
b = 0;

y = (delta_thetas .* A .* w .* c) .* exp(-(w * delta_thetas).^2) + b;

%%% Plot DoG %%%

figure('Color', 'w');

plot(delta_thetas, y,'LineWidth', 1);

xlim([-90, 90]);
ylim([-8, 8]);
ylabel('Bias (°)');
xlabel('\Delta\theta (°)');
xticks([-90:45:90]);
yticks(-8:4:8);
set(gca,'TickDir','out', 'TickLength', [0.020, 0.020], 'FontSize', 12, 'FontName', 'Helvetica');
box off;
axis square;

hold on;
line([0, 0], [min(ylim), max(ylim)], 'LineWidth', 1, 'Color', 'k');
line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');

legend({'10°', '20°', '30°'});

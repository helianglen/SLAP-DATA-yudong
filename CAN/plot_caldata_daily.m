
% plot individual calibration data as straigtlines. 

clear 

cal=readtable('Cal_table.csv'); 
x = cal{:, 2}; 
yv = cal{:, 3}; 
yh = cal{:, 4}; 

%scatter plots
sv=scatter(x, yv, 100)
sv.MarkerEdgeColor = 'b';
sv.MarkerFaceColor = [0 0.5 0.5];
sv.LineWidth=1.8;
axis([0  350  2e6 4e6]);
hold on
sh=scatter(x, yh, 100)
sh.LineWidth=1.8;
sh.MarkerEdgeColor = [0.2 0.6 0.2];
hold on

for i=2:length(yv) 
  
sv=plot([x(1) x(i)], [yv(1) yv(i)], 'b'); 
sv.MarkerEdgeColor = 'b';
sv.MarkerFaceColor = [0 0.5 0.5];
sv.LineWidth=1.8; 
axis([0  350  2e6 4e6]); 
title('Calibration data, blue: V-pol; green: H-pol')
xlabel('Temperature (K)')
ylabel('Count')
hold on
plot([x(1) x(i)], [yh(1) yh(i)], 'g') 

end 
%plot lake lines

plot([50, 150], [2.78e6, 2.78e6], 'Color', 'blue')
plot([50, 150], [2.86e6, 2.86e6], 'Color', 'green')

print('scatter_caldata_daily', '-dpng');


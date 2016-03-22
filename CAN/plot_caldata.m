
clear 

cal=readtable('Cal_table.csv'); 
x = cal{:, 2}; 
yv = cal{:, 3}; 
yh = cal{:, 4}; 

% linear fit 
fitv = polyfit(x, yv, 1); 
fith = polyfit(x, yh, 1); 

sv=scatter(x, yv, 100) 
sv.MarkerEdgeColor = 'b';
sv.MarkerFaceColor = [0 0.5 0.5];
sv.LineWidth=1.8; 
axis([0  350  2e6 4e6]); 
title('Calibration data, blue: V-pol; green: H-pol')
xlabel('Temperature (K)')
ylabel('Count')
hold on
plot(x, fitv(1)*x+fitv(2), 'b', 'Linewidth', 2); 

sh=scatter(x, yh, 100) 
sh.LineWidth=1.8; 
sh.MarkerEdgeColor = [0.2 0.6 0.2]; 
%sh.MarkerFaceColor = 'g';
hold on
plot(x, fith(1)*x+fith(2), 'g', 'Linewidth', 2); 

print('scatter_caldata', '-dpng');



clear all

T1=readtable('buoy1_all_data.csv', 'ReadVariableNames', false); 
date1= datenum(T1{:, 1}, 'mm/dd/yyyy HH:MM');
t0cm = T1{:, 2}; 
tb0cm = T1{:, 3}; 
t70cm = T1{:, 4}; 
tb70cm = T1{:, 5}; 
plot(date1, t0cm+273.16); 
hold on 
plot(date1, tb0cm); 
plot(date1, t70cm+273.16); 
plot(date1, tb70cm); 
datetick('x', 6) 
axis([-Inf Inf 80 300]) 
title('Lake temperature at buoys and Tb (K)') 
xlabel('Date') 
ylabel('Temperature (K)') 

T2=readtable('buoy2_all_data.csv', 'ReadVariableNames', false); 
date2= datenum(T2{:, 1}, 'mm/dd/yyyy HH:MM');

t0cm = T2{:, 2};
tb0cm = T2{:, 3};
t70cm = T2{:, 4};
tb70cm = T2{:, 5};
plot(date2, t0cm+273.16);
hold on
plot(date2, tb0cm);
plot(date2, t70cm+273.16);
plot(date2, tb70cm);

print('lake_t_vs_tb', '-dpng'); 


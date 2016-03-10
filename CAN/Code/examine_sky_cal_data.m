%% sky cal examiner
clear
files = dir('RAD*FB.mat');
timeall = [];
h2all = [];
v2all = [];
for i = 13:17%1:length(files)
    load(files(i).name)
    h2all = [h2all; squeeze(FullMomGroup.m2_ant(:,1,1))]; 
    v2all = [v2all; squeeze(FullMomGroup.m2_ant(:,1,1))]; 
    timeall = [timeall; time];
end

figure; plot(timeall, h2all, '.'), datetick(gca)
figure; plot(timeall, v2all, '.'), datetick(gca)

% preliminary sky cal data
% skyH = 2.76e6;
% skyV = 2.73e6;

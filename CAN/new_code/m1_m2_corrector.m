clear
files = dir('RAD*FB.mat');

for i = 1:length(files)
    files(i).name
    load(files(i).name)
    
    % H Pol
    %%
    sample_period = 500e-6; %sec
    sample_period_days = sample_period/86400;
    
    M2_I = squeeze(FullMomGroup.m2_ant(:,1,:));
    M2_Q = squeeze(FullMomGroup.m2_ant(:,2,:));
    
    M2_I = reshape(M2_I', length(time)*4,1);
    M2_Q = reshape(M2_Q', length(time)*4,1);
    
    M1_I = squeeze(FullMomGroup.m1_ant(:,1,:));
    M1_Q = squeeze(FullMomGroup.m1_ant(:,2,:));
    
    M1_I = reshape(M1_I', length(time)*4,1);
    M1_Q = reshape(M1_Q', length(time)*4,1);
    
    
    raw_M2 = (M2_I + M2_Q)/2;
    raw_M1 = (M1_I + M1_Q)/2;
    
    M2_subselect = 12;
    M1_subselect = 6;
    samples = 43968; % = 458 µs * 96 MHZ
    
    % Then scale with the subselect pointer
    
    actual_M2 = (raw_M2./samples) .* 2.^(M2_subselect);
    
    % The same has to be done with M1
    
    actual_M1 = (raw_M1./samples) .* 2.^(M1_subselect);
    
    % After this is done, the M1 component can be subtracted out with the equation
    
    h2ant = actual_M2 - actual_M1.^2 ;
    
    % % V Pol
    
    M2_I = squeeze(FullMomGroup.m2_ant(:,3,:));
    M2_Q = squeeze(FullMomGroup.m2_ant(:,4,:));
    
    M2_I = reshape(M2_I', length(time)*4,1);
    M2_Q = reshape(M2_Q', length(time)*4,1);
    
    M1_I = squeeze(FullMomGroup.m1_ant(:,3,:));
    M1_Q = squeeze(FullMomGroup.m1_ant(:,4,:));
    
    M1_I = reshape(M1_I', length(time)*4,1);
    M1_Q = reshape(M1_Q', length(time)*4,1);

    raw_M2 = (M2_I + M2_Q)/2;
    raw_M1 = (M1_I + M1_Q)/2;
    
    M2_subselect = 12;
    M1_subselect = 6;
    samples = 43968; % = 458 µs * 96 MHZ
  
    % Then scale with the subselect pointer
    
    actual_M2 = (raw_M2./samples) .* 2.^(M2_subselect);
    
    % The same has to be done with M1
    
    actual_M1 = (raw_M1./samples) .* 2.^(M1_subselect);
    
    % After this is done, the M1 component can be subtracted out with the equation
    
    v2ant = actual_M2 - actual_M1.^2 ;
    
    time_file = [time, time + sample_period_days, time + sample_period_days*2, time + sample_period_days*3];
    time_file =  reshape(time_file', numel(time_file), 1);
    
    save([files(i).name(1:end-4),'_m2data.mat'], 'h2ant', 'v2ant', 'time_file', 'pcm')
   
end


quit


%% import then combine all step attenuator data
filesDaq = dir('ADAQSTEPATTN*dlm');
[step_time, step_h, step_v] = deal([]); 
for i = 1:length(filesDaq)

  if filesDaq(i).bytes  > 0   % skip empty files
    filesDaq(i).name
    [tmp_time, tmp_h, tmp_v] = import_step(filesDaq(i).name);
    step_time = [step_time; tmp_time]; 
    step_h = [step_h; tmp_h]; 
    step_v = [step_v; tmp_v]; 
  end
end
    save step_all step_time step_h step_v 
    
quit


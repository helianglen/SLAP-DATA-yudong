%% import then combine all daq data
filesDaq = dir('ADAQTELEM*dlm');
daq_all= [];
for i = 1:length(filesDaq)

  if filesDaq(i).bytes  > 0   % skip empty files
    temp = import_ADAQ_file(filesDaq(i).name);
    daq_all = [daq_all; temp];
  end
end
    save daq_all daq_all
    

quit


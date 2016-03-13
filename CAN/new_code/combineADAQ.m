%% import then combine all daq data
filesDaq = dir('ADAQTELEM*dlm');
[timeDaq, IMA_H_SWITCH,IMA_V_SWITCH, TempDiplexerH, TempDiplexerV, TempLNAH, TempLNAV, TempAntennaH, TempAntennaV ] = deal([]); 
for i = 1:length(filesDaq)

  if filesDaq(i).bytes  > 0   % skip empty files
    filesDaq(i).name
    [timeDaq_t, IMA_H_SWITCH_t,IMA_V_SWITCH_t, TempDiplexerH_t, TempDiplexerV_t, TempLNAH_t, TempLNAV_t, TempAntennaH_t, TempAntennaV_t ] = import_ADAQ_file(filesDaq(i).name);
    timeDaq = [timeDaq; timeDaq_t]; 
    IMA_H_SWITCH = [IMA_H_SWITCH; IMA_H_SWITCH_t]; 
    IMA_V_SWITCH = [IMA_V_SWITCH; IMA_V_SWITCH_t]; 
    TempDiplexerH = [TempDiplexerH; TempDiplexerH_t]; 
    TempDiplexerV = [TempDiplexerV; TempDiplexerV_t]; 
    TempLNAH = [TempLNAH; TempLNAH_t]; 
    TempLNAV = [TempLNAV; TempLNAV_t]; 
    TempAntennaH = [TempAntennaH; TempAntennaH_t]; 
    TempAntennaV = [TempAntennaV; TempAntennaV_t]; 
  end
end
    save daq_all timeDaq IMA_H_SWITCH IMA_V_SWITCH TempDiplexerH TempDiplexerV TempLNAH TempLNAV TempAntennaH TempAntennaV
    
quit


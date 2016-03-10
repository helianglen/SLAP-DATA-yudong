%% This script manually removes 1.09 Hz noise in SLAP data 
% This is required from Feb 22nd, 2014 flight back throughout history of SLAP flights. 

filesRAD = dir('RAD*m2data.mat');

ind = 1822;
window_length = 4;
for i = 1:6%24:length(filesRAD)
    filesRAD(i).name
    load(filesRAD(i).name)
    h2new = h2ant;
    v2new = v2ant;
    first_ind = find(~isnan(v2ant),1);
    ind_cycle = first_ind:(first_ind + ind);

    while ~isempty(first_ind) && ind_cycle(end) < length(v2ant)
        % find indices of first cycle of real data and 1.09 Hz noise
        cycle = v2ant(ind_cycle);

        ind_noise = find(cycle<nanmean(cycle)-nanstd(cycle));
        if ~isempty(ind_noise)
            ind_first = [];
            ind_last = [];
            ind_diff = diff(ind_noise);
            % find indices of at least 3 consecutive values below the threshold
            % which would signal a period of 1.09 Hz noise
            for j = 1:length(ind_diff)-window_length+1
                ind_window = ind_diff(j:j+window_length-1);
                ind_consec = find(ind_window == 1);
                if length(ind_consec) == window_length
                    ind_first = j;
                    break
                end
            end
            for j = length(ind_diff):-1:window_length
                ind_window = ind_diff(j-window_length+1:j);
                ind_consec = find(ind_window > 15);
                if isempty(ind_consec)
                    ind_last = j;
                    break
                end
            end
            % add an index to ind_last because these indices are deduced using
            % the diff function which subtracts an index value
            h2new(ind_cycle(ind_noise(ind_first:ind_last+1))) = nan;
            v2new(ind_cycle(ind_noise(ind_first:ind_last+1))) = nan;
        end
        ind_cycle = ind + ind_cycle;
    end
    h2new_v2 = h2new;
    v2new_v2 = v2new;
    save([filesRAD(i).name], '*new_v2', '-append')
end

%%
figure
subplot(211)
plot(v2ant, '.')
hold on
plot(v2new_v2, 'r.')

subplot(212)
plot(h2ant, '.')
hold on
plot(h2new_v2, 'r.')
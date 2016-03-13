function timeRad = repair_time_tags(timeRad)

% The timetags in the radiometer file  have some peculiar traits. They 
% jump around ephemerally every once in a while, both to lower
% and higher values. Additionally, they jump up permanently every once in a
% while so the final time is larger in value than what would be
% expected if the time was steadily increasing every 2 ms as it should be.

% to counteract this, the following code attempts to 'repair' the rad
% timetags to make them smoother, while maintaining the occasional 'leaps'
% in time because it is evident they correspond to leaps in the m2 data 

% first create a rad timetag that linearly increases every 0.5 ms from the
% start time of timeRad
timeRadLinear = timeRad(1):0.0005/86400:(timeRad(1) + (length(timeRad)-1)*0.0005/86400);
timeRadLinear = timeRadLinear';

% initialize jump variable because upcoming while-loop needs jump to not be
% an empty vector
jump = 1;
jumpind = 0;
indnext = 1;
while 1
    interval = (timeRad(indnext:end) - timeRadLinear(indnext:end));
    jump = find(interval > 0.3/86400, 4000);
    jumpdbl = find(interval > 1.5/86400, 50);
    if ~isempty(jumpdbl) && jumpdbl(1) > jump(end)
        jumpdbl = [];
    end
    jumpind = find(diff(jump)>1);
        % when there is a permanent jump (usually of 1 second)  in the data, add 1 second to timeRadLinear.   
    if (~isempty(jumpind) && jumpind(1) > 300) || (~isempty(jumpind) && length(jumpdbl)>15)  ||  (isempty(jumpind) && length(jump)> 400)  
        % The first component in the
        % if-statement makes sure the jump is not ephemeral,
        % meaning the timetags dont go back down to the original
        % increments. This is found by checking whether
        % the jump vector is larger than  400 points meaning it
        % isn't ephemeral and that the time indices dont go back down to the
        % original value within 400 data points.
        if ~isempty(jumpind) && length(jumpdbl) < 15 && jumpind(1) > 300
            newint = indnext + jump(1);
            jumptime = nanmean(interval(jump(1:jumpind(1))));
            indjump = jump(jumpind(1));
            % The second component checks whether there is
            % a gap between the jumps in timetags within the 4000
            % data points as well and if there is a double jump. If so,
            % add 1 second to the timetags at the beginning of
            % the data span then reset to the end of the first data span before the
            % second jump so another second could be added on.
        elseif ~isempty(jumpind) && length(jumpdbl) > 15
            jumpvalue = find(jumpind>400,1);
            jumpval = jumpind(jumpvalue);
            newint = indnext + jump(jumpval);
            jumptime =  nanmean(interval(jump(jumpval:end)));
            indjump = jump(jumpval);
            % The third component is the case where
            % jumpind is empty, meaning the jump in timetags is steady for at least 400
            % data points, so just add 1 second to the beginning of the data span.
        elseif length(jump)>400
            newint = indnext + jump(1);
            jumptime = nanmean(interval(jump));
            indjump = jump(end);
        end
        timeold = timeRadLinear;
        timeRadLinear(newint:end) = timeRadLinear(newint:end) + jumptime;
        timenew = timeRadLinear;
    else
        % otherwise, this is a jump in the rad timetags that comes back to
        % the orig value. In this case, don't increase the new timetag by any amount
        timenew = timeRadLinear;
        timeold = timeRadLinear;
    end
    % break out of while loop if there is no more jumps or if the jumps are
    % ephemeral. also break out if there is very little (<0.1 seconds)
    % change in the time index file.
    
    if abs(diff([timenew(end), timeold(end)])) < 0.1/86400 ||  length(jump) < 400
        break
    end
    
    indnext = indnext + indjump;
end

% reset timeRad to be new time series created in last while loop
timeRad = timeRadLinear;
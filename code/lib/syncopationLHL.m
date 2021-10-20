function syncopationScore = syncopationLHL(pattern, meter, eventsInCycle, varargin)

% "If N is a note that precedes a rest, R, and R has a metric weight greater than or equal to N, 
% then the pair (N, R) is said to constitute a monophonic syncopation."
% 
% If N < R
% Syncopation = R - N

% metric template (pulse periods) is specified by a string: 
%     '2_4'
%     x . . . .
%     x . x . x 
% 
%     '2_6'
%     x . . . . . 
%     x . x . x . 
% 
%     '3_6'
%     x . . . . . 
%     x . . x . . 
% 
% by default, the function will return syncopation score summed across the
% whole input pattern (even when longer than eventsInCycle)
% 
% with "perbar" option (varargin), the function will sum syncopation score 
% separately over each cycle in the input



if strcmpi(meter,'2_4')
    salience = [0,-2,-1,-2]; 
    weightgrid = repmat(salience, 1, ceil(size(pattern,2)/length(salience))); 
elseif strcmpi(meter,'2_6')
    salience = [0,-2,-1,-2,-1,-2]; 
    weightgrid = repmat(salience, 1, ceil(size(pattern,2)/length(salience))); 
elseif strcmpi(meter,'3_6')
    salience = [0,-2,-2,-1,-2,-2]; 
    weightgrid = repmat(salience, 1, ceil(size(pattern,2)/length(salience))); 
else
    error('meter string specified incorrectly')
end
    
    
if any(strcmpi(varargin, 'perbar'))
    
    % get number of cycles (bars) in the input
    nBars = length(pattern)/eventsInCycle; 

    % allocate, dim: [file x bar]
    syncopationScore = zeros(size(pattern,1),nBars); 
    
    for filei = 1:size(pattern,1)
        synidx = 0; 
        bar = 1; 
        for i=1:size(pattern,2)
            if pattern(filei,i)
                c=1; 
                tmpx = []; 
                while 1
                    if (i+c)>size(pattern,2)
                        break
                    elseif pattern(filei,i+c)==1
                        break
                    elseif pattern(filei,i+c)==0
                        tmpx = [tmpx, weightgrid(i+c)]; 
                    end
                    c = c+1; 
                end
                if any(tmpx>weightgrid(i))
                    synidx = synidx + (max(tmpx)-weightgrid(i)); 
                end
            end
            if mod(i+1,eventsInCycle)==0
               syncopationScore(filei, bar) = synidx; 
               synidx = 0; 
               bar = bar+1;  
            end
        end
    end
        
else
       
    syncopationScore = zeros(1,size(pattern,1)); 

    for filei=1:size(pattern,1)
        synidx = 0; 
        for i=1:size(pattern,2)
            if pattern(filei,i)
                c=1; 
                tmpx = []; 
                while 1
                    if (i+c)>size(pattern,2)
                        break
                    elseif pattern(filei,i+c)==1
                        break
                    elseif pattern(filei,i+c)==0
                        tmpx = [tmpx, weightgrid(i+c)]; 
                    end
                    c = c+1; 
                end
                if any(tmpx>weightgrid(i))
                    synidx = synidx + (max(tmpx)-weightgrid(i)); 
                end
            end
        end

        syncopationScore(filei) = synidx; 
    end

end




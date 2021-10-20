function obj = constructStim(obj)
% construct rhythmic sequences based on requested number of on-beat events
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get beat frequency based on the requested period
obj.meterFrex = 1/(obj.beatPeriod*obj.gridIOI); 

% get which grid positions are on-the-beat
onbeatPos = [1:obj.beatPeriod:obj.nEvents]; 

% get which grid positions are off-the-beat
offbeatPos = [1:obj.nEvents]; 
offbeatPos(ismember(offbeatPos,onbeatPos)) = []; 

% if requested more onbeat events than available beat positions, issue a
% warning and correct to the largest possible value
if obj.nSoundsOnbeat>length(onbeatPos)
    warning(sprintf('you want more sounds to be on-beat than there are available on-beat positions...\n\ncorrecting to maximum posible number on-beat = %d',length(onbeatPos)));
    obj.nSoundsOnbeat = length(onbeatPos); 
end

% get number of sounds that will be offbeat
obj.nSoundsOffbeat = round(obj.nSounds-obj.nSoundsOnbeat); 

% if requested more offbeat events than available positions, give error
if obj.nSoundsOnbeat>length(onbeatPos)
    error(sprintf('you want more sounds to be off-beat than there are available off-beat positions...')); 
end


obj.patterns = zeros(obj.nTrials,obj.nEvents); 

% go over trials
for triali=1:obj.nTrials

    % construct the pattern in this trial
    while 1

        % randomly select onbeat and offbeat positions where sounds
        % will be placed
        soundPositions = sort([randsample(onbeatPos,obj.nSoundsOnbeat), ...
                           randsample(offbeatPos,obj.nSoundsOffbeat)]); 

        % continue if there is: 
        %     1) no silent gaps longer than allowed 
        %     2) first event is sound 
        % repeat the random selection of sound positions otherwise 
        if ~any(diff(soundPositions)>obj.maxEventsSilence) & any(soundPositions==1)
            % assign sounds to sound positions in grid representation
            obj.patterns(triali,soundPositions) = 1; 
            break
        end
    end

end


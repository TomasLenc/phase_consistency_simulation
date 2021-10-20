function C = syncopationPE(in, grouping, varargin)
% C = omega*O+U
% 
% counterevidence (C) against the induction of a meter is determined as 
% the weighted (omega) sum of omitted (O) and unaccented (U) downbeats
% 
% omega = 4:1 (weighting of silent vs. unaccented downbeats)
%

W = 4; 

C = zeros(1,size(in,1)); 

for filei=1:size(in,1)
    
    pattern         = in(filei,:); 
    accents         = zeros(size(pattern)); 
    N               = length(pattern);     
    
    % pad the pattern from front and back (this will only be used to
    % estimate accents)
    if any(strcmpi(varargin,'zeropad')) 
        % pad with zeroes (same as thinking of the pattern as not being
        % preceded or followed by anything)
        padding = zeros(1,size(in,2)); 
        padded_pattern  = [padding, pattern, padding]; 
    else
        % assume circularity (default)
        padded_pattern  = [pattern, pattern, pattern]; 
    end
    
    
    for i = [1:length(pattern)]+N
        
        % if there is sound, consider it an accent if: 
        if padded_pattern(i)
            
            % it's followed by silence (i.e. last in a group or alone)
            if padded_pattern(i+1)==0
                accents(i-N)=1; 
                continue
            end
            
            % or it's first in group of >=3 tones 
            % (silence at i-1, sound at i+1 and at i+2)
            if (padded_pattern(i-1)==0 & padded_pattern(i+1)==1 & padded_pattern(i+2)==1)
                accents(i-N)=1; 
                continue
            end
            
        end
    end

    coinc_accent    = sum(accents(1:grouping:end)==1); 
    coinc_unaccent  = sum(pattern(1:grouping:end)==1 & accents(1:grouping:end)==0); 
    coinc_silence   = sum(pattern(1:grouping:end)==0); 

    C(filei) = (W * coinc_silence) + (coinc_unaccent); 
end




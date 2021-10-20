function plot2subSimul(sub1, sub2, varargin)

    % this is a main function that plots multipanel figure with
    % kernels, time-domain, and frequency-domain representation

    %%%%%% figure
    f = figure('color','white','position', [436 559 1072 208]); 

    pnl = panel(f);

    pnl.pack('h',[10,20,35,35]); 

    pnl(1).pack('v',2); % time on top, freq on bottom
    
    pnl(2).pack('v',2); 
    
    pnl(4).pack('h',2); 

    
    % -------------------------------------------------

    %%%%%% IR
    if numel(sub1.ir)==1 % if only one kernel, plot with neutral color
        cols = {sub1.irCol}; 
    elseif numel(sub1.ir)==2 % if two kernels, plot with base/target colors
        cols = {sub1.col_f_target, sub1.col_f_base}; 
    end

    iri = 1; 

    % time
    y = sub1.ir.ir; 

    xLims = [floor(min([sub1.ir.T])*sub1.prec)/sub1.prec, ...
                 ceil(max([sub1.ir.T])*sub1.prec)/sub1.prec]; 

    yLims = [floor(min([y])*sub1.prec)/sub1.prec, ...
                 ceil(max([y])*sub1.prec)/sub1.prec]; 


    ax = pnl(1,1).select(); 
    plot(ax, sub1.ir.T, y,'color',cols{iri},'linew',sub1.linew_t)
    box off
    ax.XLim = xLims; 
    ax.YLim = yLims; 
    ax.XTick = xLims; 
    ax.YTick = yLims; 
    ax.YAxis.Visible = 'off'; 

    % freq
    ax = pnl(1,2).select(); 
    stem(ax, sub1.ir.freq, sub1.ir.mX, ...
         'color', cols{iri},...
         'linew',sub1.linew_f, ...
         'marker','none')
    box off
    ax.XLim = [0, sub1.maxFreqLim]; 
    ax.YLim = [0,max(sub1.ir.mX)]; 
    ax.XTick = [0, sub1.maxFreqLim]; 
    ax.YTick = [0,max(sub1.ir.mX)]; 
    ax.YAxis.Visible = 'off'; 
        
        



    
 %%   time domain
      
 
    %%%%%% TD
    triali = 1; 
    colSub1 = sub1.col_t; 
    colSub2 = [157, 18, 222]/255; 

    % sub 1
    if ~isfield(sub1.eeg,'erp')
        sub1 = sub1.analEEG('erp'); 
    end         
    y1 = squeeze(mean(mean(sub1.eeg.erp.erpCycle,1),2));
    y1 = y1-mean(y1); 
    
    t = [0:length(y1)-1]/sub1.fs; 

    xLims = [0, sub1.nEventsCycle*sub1.gridIOI]; 
    xLimsIdx = [1, round(sub1.nEventsCycle*sub1.gridIOI*sub1.fs)]; 
    yLims1 = [floor(min(y1(xLimsIdx(1):xLimsIdx(2))) * sub1.prec)/sub1.prec, ...
             ceil(max(y1(xLimsIdx(1):xLimsIdx(2))) *sub1.prec)/sub1.prec]; 

    soundPos    = (find(sub1.patterns(triali,:))-1) * sub1.gridIOI; 
    onBeatPos   = (find(sub1.beatPos(triali,:))-1) * sub1.gridIOI; 
    offBeatPos  = (find(~sub1.beatPos(triali,:))-1) * sub1.gridIOI; 

    if isempty(yLims1)
        yRange = max(y1)-min(y1); 
        yLims1 = [min(y1)-yRange*0.3, ...
                  max(y1)+yRange*0.3]; 
    end
        
    % sub 2    
    if ~isfield(sub2.eeg,'erp')
        sub2 = sub2.analEEG('erp'); 
    end   
    y2 = squeeze(mean(mean(sub2.eeg.erp.erpCycle,1),2)); 
    y2 = y2-mean(y2); 
    
    yLims2 = [floor(min(y2(xLimsIdx(1):xLimsIdx(2))) * sub2.prec)/sub2.prec, ...
             ceil(max(y2(xLimsIdx(1):xLimsIdx(2))) *sub2.prec)/sub2.prec]; 
    if isempty(yLims2)
        yRange = max(y2)-min(y2); 
        yLims2 = [min(y2)-yRange*0.3, ...
                  max(y2)+yRange*0.3]; 
    end
        
    % get overall ylimits
    yLims = [min([yLims1,yLims2]), max([yLims1,yLims2])]; 
          
    % plot
    pnl(3).select(); 
    ax = gca; 
    plot(t, y1, 'color', colSub1, 'linew', sub1.linew_t)
    hold on
    plot(t, y2, 'color', colSub2, 'linew', sub1.linew_t)

    % plot sound positions
    plot([soundPos;soundPos]', yLims', ':', ...
       'color',sub1.col_f_neutral, ...
       'linew',1.5)

    %%%%%% plot meter
    colMeter = [235, 149, 52]/255; 
    % for plotting pulse positions, leave some vertical and horizontal space for the 
    % rectangles and arches
    archSpaceX = 0.1 * sub1.gridIOI; 
    meterPhase = 0; 
    pulseArch = cell(1,length(sub1.meterPeriods)); 
    for archi=1:length(sub1.meterPeriods)
        % get number of samples for each pulse arch
        pulseArchN = round( (sub1.meterPeriods(archi) .* sub1.gridIOI - archSpaceX) .* sub1.fs ); 
        % make it odd
        pulseArchN = pulseArchN + mod(pulseArchN+1,2); 
        % get time vector centered around 0
        tPulseArch = [ -(floor(pulseArchN/2)) : (floor(pulseArchN/2))]; 
        % we choose sigma that's too large, so within the window the gaussian
        % doens't go all the way to 0. It's ad hoc hack, looks fine. 
        sigma = round(sub1.meterPeriods(archi)*0.9 * sub1.gridIOI * sub1.fs); 
        % make gaussian
        pulseArch{archi} = exp(-tPulseArch.^2./(2*sigma^2)); 
        % scale its amplitude
        pulseArch{archi} = 1 * pulseArch{archi}; 
    end
    % move arch to 0 on y axis
    pulseArch = cellfun(@(x) x-min(cellfun(@min,pulseArch)), pulseArch, 'uni',0); 
    archSpaceY = yLims(2)*1 + max(cellfun(@range,pulseArch))*[0:length(sub1.meterPeriods)-1]; 
    % go over each pulse period
    for periodi=1:length(sub1.meterPeriods)
        % find pulse positions based on requested period and phase 
        pulsePos = [meterPhase : sub1.meterPeriods(periodi) : sub1.nEvents-1]; 
        pulsePosTime = pulsePos * sub1.gridIOI; 
        % trim at maxTlim seconds for visualisation
        pulsePosTime(pulsePosTime>t(end)) = []; 
        % go over each pulse
        for pulsi=1:length(pulsePosTime)
            % find the index of the pulse position
            pulsePosIdx = round(pulsePosTime(pulsi)*sub1.fs); 
            % plot arch starting at pulse position
            plot( t(pulsePosIdx+1:pulsePosIdx+length(pulseArch{periodi})), ...
                  archSpaceY(periodi) + pulseArch{periodi}, ...
                      'LineWidth',2, ...
                      'LineStyle','-', ...
                      'Color', colMeter)            
        end
    end    
    yLims(2) = max(archSpaceY(periodi) + pulseArch{periodi}); 

    box off
    ax = gca; 
    ax.YAxis.Visible = 'off'; 
    ax.XLim = xLims; 
    ax.XTick = ax.XLim; 
    ax.YLim = yLims; 

    
    leg = legend({'condition 1','condition 2'}); 
    leg.Box = 'off'; 
    leg.Position = [0.4 0.7990 0.1268 0.2036]; 
    
    
    
    
    
    
    
    %%
    %% plot stim 
    %%
    %%
    

    % sub 1
    y1 = sub1.stim.x(1,:); 
    t = [0:length(y1)-1]/sub1.fs; 

    xLims = [0, sub1.nEventsCycle*sub1.gridIOI]; 
    xLimsIdx = [1, round(sub1.nEventsCycle*sub1.gridIOI*sub1.fs)]; 
    yLims = [floor(min(y1(xLimsIdx(1):xLimsIdx(2))) * sub1.prec)/sub1.prec, ...
             ceil(max(y1(xLimsIdx(1):xLimsIdx(2))) *sub1.prec)/sub1.prec]; 
    if isempty(yLims)
        yRange = max(y1)-min(y1); 
        yLims = [min(y1)-yRange*0.3, ...
                  max(y1)+yRange*0.3]; 
    end
    
    % plot
    pnl(2,1).select(); 
    ax = gca; 
    plot(t, y1, 'color', colSub1, 'linew', sub1.linew_t)
    hold on 

    %%%%%% plot meter
    colMeter = [235, 149, 52]/255; 
    % for plotting pulse positions, leave some vertical and horizontal space for the 
    % rectangles and arches
    archSpaceX = 0.1 * sub1.gridIOI; 
    meterPhase = 0; 
    pulseArch = cell(1,length(sub1.meterPeriods)); 
    for archi=1:length(sub1.meterPeriods)
        % get number of samples for each pulse arch
        pulseArchN = round( (sub1.meterPeriods(archi) .* sub1.gridIOI - archSpaceX) .* sub1.fs ); 
        % make it odd
        pulseArchN = pulseArchN + mod(pulseArchN+1,2); 
        % get time vector centered around 0
        tPulseArch = [ -(floor(pulseArchN/2)) : (floor(pulseArchN/2))]; 
        % we choose sigma that's too large, so within the window the gaussian
        % doens't go all the way to 0. It's ad hoc hack, looks fine. 
        sigma = round(sub1.meterPeriods(archi)*0.9 * sub1.gridIOI * sub1.fs); 
        % make gaussian
        pulseArch{archi} = exp(-tPulseArch.^2./(2*sigma^2)); 
        % scale its amplitude
        pulseArch{archi} = 1 * pulseArch{archi}; 
    end
    % move arch to 0 on y axis
    pulseArch = cellfun(@(x) x-min(cellfun(@min,pulseArch)), pulseArch, 'uni',0); 
    archSpaceY = yLims(2)*1.1 + max(cellfun(@range,pulseArch))*[0:length(sub1.meterPeriods)-1]; 
    % go over each pulse period
    for periodi=1:length(sub1.meterPeriods)
        % find pulse positions based on requested period and phase 
        pulsePos = [meterPhase : sub1.meterPeriods(periodi) : sub1.nEvents-1]; 
        pulsePosTime = pulsePos * sub1.gridIOI; 
        % trim at maxTlim seconds for visualisation
        pulsePosTime(pulsePosTime>t(end)) = []; 
        % go over each pulse
        for pulsi=1:length(pulsePosTime)
            % find the index of the pulse position
            pulsePosIdx = round(pulsePosTime(pulsi)*sub1.fs); 
            % plot arch starting at pulse position
            plot( t(pulsePosIdx+1:pulsePosIdx+length(pulseArch{periodi})), ...
                  archSpaceY(periodi) + pulseArch{periodi}, ...
                      'LineWidth',2, ...
                      'LineStyle','-', ...
                      'Color', colMeter)            
        end
    end    
    yLims(2) = max(archSpaceY(periodi) + pulseArch{periodi}); 

    box off
    ax = gca; 
    ax.YAxis.Visible = 'off'; 
    ax.XLim = xLims; 
    ax.XTick = ax.XLim; 
    ax.YLim = yLims; 
    
    
    
    
    
    %%
    
    
    % sub 2    
    y2 = sub2.stim.x(1,:); 
    yLims = [floor(min(y2(xLimsIdx(1):xLimsIdx(2))) * sub2.prec)/sub2.prec, ...
             ceil(max(y2(xLimsIdx(1):xLimsIdx(2))) *sub2.prec)/sub2.prec]; 
    if isempty(yLims)
        yRange = max(y2)-min(y2); 
        yLims = [min(y2)-yRange*0.3, ...
                  max(y2)+yRange*0.3]; 
    end
    
    % plot
    pnl(2,2).select(); 
    ax = gca; 
    plot(ax, t, y2, 'color', colSub2, 'linew', sub1.linew_t)
    hold on 
    
    %%%%%% plot meter
    colMeter = [235, 149, 52]/255; 
    % for plotting pulse positions, leave some vertical and horizontal space for the 
    % rectangles and arches
    archSpaceX = 0.1 * sub1.gridIOI; 
    meterPhase = 0; 
    pulseArch = cell(1,length(sub1.meterPeriods)); 
    for archi=1:length(sub1.meterPeriods)
        % get number of samples for each pulse arch
        pulseArchN = round( (sub1.meterPeriods(archi) .* sub1.gridIOI - archSpaceX) .* sub1.fs ); 
        % make it odd
        pulseArchN = pulseArchN + mod(pulseArchN+1,2); 
        % get time vector centered around 0
        tPulseArch = [ -(floor(pulseArchN/2)) : (floor(pulseArchN/2))]; 
        % we choose sigma that's too large, so within the window the gaussian
        % doens't go all the way to 0. It's ad hoc hack, looks fine. 
        sigma = round(sub1.meterPeriods(archi)*0.9 * sub1.gridIOI * sub1.fs); 
        % make gaussian
        pulseArch{archi} = exp(-tPulseArch.^2./(2*sigma^2)); 
        % scale its amplitude
        pulseArch{archi} = 1 * pulseArch{archi}; 
    end
    % move arch to 0 on y axis
    pulseArch = cellfun(@(x) x-min(cellfun(@min,pulseArch)), pulseArch, 'uni',0); 
    archSpaceY = yLims(2)*1.1 + max(cellfun(@range,pulseArch))*[0:length(sub1.meterPeriods)-1]; 
    % go over each pulse period
    for periodi=1:length(sub1.meterPeriods)
        % find pulse positions based on requested period and phase 
        pulsePos = [meterPhase : sub1.meterPeriods(periodi) : sub1.nEvents-1]; 
        pulsePosTime = pulsePos * sub1.gridIOI; 
        % trim at maxTlim seconds for visualisation
        pulsePosTime(pulsePosTime>t(end)) = []; 
        % go over each pulse
        for pulsi=1:length(pulsePosTime)
            % find the index of the pulse position
            pulsePosIdx = round(pulsePosTime(pulsi)*sub1.fs); 
            % plot arch starting at pulse position
            plot( t(pulsePosIdx+1:pulsePosIdx+length(pulseArch{periodi})), ...
                  archSpaceY(periodi) + pulseArch{periodi}, ...
                      'LineWidth',2, ...
                      'LineStyle','-', ...
                      'Color', colMeter)            
        end
    end    
    yLims(2) = max(archSpaceY(periodi) + pulseArch{periodi}); 

    box off
    ax = gca; 
    ax.YAxis.Visible = 'off'; 
    ax.XLim = xLims; 
    ax.XTick = ax.XLim; 
    ax.YLim = yLims; 
    
        
    


    
    
    
    
    %% frequency domain
    
    % find ylimits 
    
    mX = sub1.eeg.magn.mXavgTime; 
    freq = sub1.eeg.magn.freq;
    idxMinFrex = dsearchn(freq',0.3); 
    xLims = [0, sub1.maxFreqLim]; 
    yLims1 = [floor(min(mX(idxMinFrex:end))*sub1.prec)/sub1.prec, ...
             ceil(max(mX(idxMinFrex:end))*sub1.prec)/sub1.prec]; 
    
    mX = sub2.eeg.magn.mXavgTime; 
    yLims2 = [floor(min(mX(idxMinFrex:end))*sub1.prec)/sub1.prec, ...
             ceil(max(mX(idxMinFrex:end))*sub1.prec)/sub1.prec]; 
   
    yLims = [min([yLims1,yLims2]), max([yLims1,yLims2])]; 
    
    
    
    % -------------------------------------------------

    %%%%%% FD
    if ~isfield(sub1.eeg,'magn')
        sub1 = sub1.analEEG('magnitude'); 
    end
    if ~isfield(sub2.eeg,'magn')
        sub2 = sub2.analEEG('magnitude'); 
    end

    
    %%%% sub 1
    mX = sub1.eeg.magn.mXavgTime; 

    pnl(4,1).select(); 
    ax = gca; 

    % plot whole spectrum
    h = stem(freq, mX, 'color', sub1.col_f_neutral, 'linew',sub1.linew_f, 'marker','none'); 
    hold on 

    % 12 frequencies of interest (just a sanity check, all frex
    % should be assigned to meterRel or meterUnrel anyway)
    h = stem(freq(sub1.frexIdx), mX(sub1.frexIdx), ...
             'color', 'k',...
             'linew',sub1.linew_f, ...
             'marker','none'); 

    % meter-related frequecnies
    h = stem(freq(sub1.frexIdx(sub1.frexMeterRel_idx)), ...
               mX(sub1.frexIdx(sub1.frexMeterRel_idx)), ...
             'color', sub1.col_f_target, ...
             'linew',sub1.linew_f, ...
             'marker','none'); 

    % meter-related frequecnies
    h = stem(freq(sub1.frexIdx(sub1.frexMeterUnrel_idx)), ...
               mX(sub1.frexIdx(sub1.frexMeterUnrel_idx)), ...
             'color', sub1.col_f_base, ...
             'linew',sub1.linew_f, ...
             'marker','none'); 

    xlim(xLims)
    ylim(yLims)
    box off
    ax = gca; 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
    ax.XTick = ax.XLim; 

    
    
    %%%% sub 2    
    mX = sub2.eeg.magn.mXavgTime; 

    pnl(4,2).select(); 
    ax = gca; 

    % plot whole spectrum
    h = stem(freq, mX, 'color', sub2.col_f_neutral, 'linew',sub2.linew_f, 'marker','none'); 
    hold on 

    % 12 frequencies of interest (just a sanity check, all frex
    % should be assigned to meterRel or meterUnrel anyway)
    h = stem(freq(sub2.frexIdx), mX(sub2.frexIdx), ...
             'color', 'k',...
             'linew',sub2.linew_f, ...
             'marker','none'); 

    % meter-related frequecnies
    h = stem(freq(sub2.frexIdx(sub2.frexMeterRel_idx)), ...
               mX(sub2.frexIdx(sub2.frexMeterRel_idx)), ...
             'color', sub2.col_f_target, ...
             'linew',sub2.linew_f, ...
             'marker','none'); 

    % meter-related frequecnies
    h = stem(freq(sub2.frexIdx(sub2.frexMeterUnrel_idx)), ...
               mX(sub2.frexIdx(sub2.frexMeterUnrel_idx)), ...
             'color', sub2.col_f_base, ...
             'linew',sub2.linew_f, ...
             'marker','none'); 

    xlim(xLims)
    ylim(yLims)
    box off
    ax = gca; 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
    ax.XTick = ax.XLim; 

    
    
    
   %% 
    
    
    
    %%%%%% labels
    pnl.fontsize = sub1.fontsize; 

%     pnl(1,1).title('response kerenel'); 

%     ax = pnl(2).select(); 
%     text(ax, mean(ax.XLim),ax.YLim(1)-abs(diff(ax.YLim))*0.1, 'time (s)', ...
%         'fontsize',sub1.fontsize, ...
%         'HorizontalAlignment','center'); 

%     ax = pnl(3).select(); 
%     text(ax, mean(ax.XLim),ax.YLim(1)-abs(diff(ax.YLim))*0.1, 'frequency', ...
%         'fontsize',sub1.fontsize, ...
%         'HorizontalAlignment','center'); 

    pnl(4,1).title(sprintf('condition 1\nzMeter %.3g', ...
                            sub1.eeg.magn.zMeterRelAvgTime)); 

    pnl(4,2).title(sprintf('condition 2\nzMeter %.3g', ...
                            sub2.eeg.magn.zMeterRelAvgTime)); 


    %%%%%% margins
%     pnl(1,1).de.margin = [3,1,2,10]; 

    pnl.margin = [15,15,5,8]; 
    pnl.marginbottom = 10; 
    pnl.margintop = 15; 



end      


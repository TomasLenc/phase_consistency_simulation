function simulateEegExp(patterns,fs, DO_PLOT,DO_SAVE_FIG)





















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if DO_PLOT
    

    mX_ylim = [-0.1,1]; 
    TD_ylim = max(abs([eeg_cond1(:);eeg_cond2(:)])); 

    col_cond1 = [102, 153, 0]/255; 
    col_cond2 = [102, 0, 204]/255; 

    fontsize = 16; 
    linew = 1.7; 






    %% plot example trial env (transform vs. no-transform)

    % plot envelope from one example trial for EEG without any transformation
    % (one-to-one stimulus tracking) and also for EEG with neural empahsis at
    % on-beat positions
    f = figure('color','white','position',[127 146 861 480]); 
    pnl = panel(f); 
    pnl.pack({40,60}); 
    pnl(1).pack('h',3); 


    pnl(1,1).select()
    plot(t,stim(1,:),'color','k', 'linew',linew)
    tit = pnl(1,1).title(sprintf('stimulus \nrepresentation')); 
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,t(end)], 'ylim',[0,TD_ylim], 'fontsize',fontsize)

    pnl(1,2).select()
    plot(irT,ir,'color',[255, 153, 0]/255,'linew',linew)
    tit = pnl(1,2).title('impulse response'); 
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,irT(end)], 'fontsize',fontsize)

    pnl(1,3).select()
    plot(t,stim_cond1(1,:),'color','r', 'linew',linew)
    tit = pnl(1,3).title(sprintf('stimulus \nrepresentation \n(with emphasis)')); 
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,t(end)], 'ylim',[0,TD_ylim], 'fontsize',fontsize)



    pnl(2).select()
    plot(t,eeg_cond0(1,:),'color','k', 'linew',linew)
    hold on
    plot(t,eeg_cond1(1,:),'color','r', 'linew',linew)
    idx = find(stim(1,:)); 
    plot([t(idx);t(idx)],[-100;100],':','color',[145,125,145]/255,'linew',linew); 
    tit = pnl(2).title('example trial'); 
    l=legend({'no transformation','emphasis on beat'}); 
    l.Position = [0.7950 0.4611 0.1823 0.0823]; 
    l.Box = 'off'; 
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,t(end)], 'ylim',[-TD_ylim,TD_ylim], 'fontsize',fontsize)


    xlab = pnl(1).xlabel('Time (s)'); 
    ylab = pnl(1).ylabel('Amplitude (a.u.)'); 

    xlab = pnl(2).xlabel('Time (s)'); 
    ylab = pnl(2).ylabel('Amplitude (a.u.)'); 

    pnl.de.margin = 2; 
    pnl(2).margintop = 25;
    pnl.margin = [13 10 2 15];

    pnl.fontsize = fontsize; 






    %% plot main figure


    fig_name = sprintf('%d sounds onbeat | %d sounds offbeat | %d events total | brain emphasis %.1f',nSoundsOnbeat, nSoundsOffbeat, nEvents, onbeatBrainEmphasis); 

    f = figure('color','white','position',[154 192 1444 799],'name',fig_name); 


    % ================================================================================
    % plot envelope averaged across trials


    subplot(3,4,1)
    plot(t,mean(eeg_cond1,1),'color',col_cond1, 'linew',linew)
    xlabel('time (s)')
    title('cond1 mean eeg across trials')
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,t(end)], 'ylim',[-TD_ylim,TD_ylim], 'fontsize',fontsize)

    subplot(3,4,2)
    plot(t,mean(eeg_cond2,1),'color',col_cond2, 'linew',linew)
    xlabel('time (s)')
    title('cond2 mean eeg across trials')
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,t(end)], 'ylim',[-TD_ylim,TD_ylim], 'fontsize',fontsize)



    subplot(3,4,3)
    hold on
    plot(t,mean(eeg_cond1,1),'color',col_cond1, 'linew',linew)
    plot(t,mean(eeg_cond2,1),'color',col_cond2, 'linew',linew)
    xlabel('time (s)')
    % title('cond1 vs. cond2')
    l = legend({'cond 1 (consistent)','cond 2 (random)'}); 
    l.Box = 'off';  
    l.Position = [0.5758    0.8867    0.1139    0.0494]; 
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,t(end)], 'ylim',[-TD_ylim,TD_ylim], 'fontsize',fontsize)





    % ================================================================================
    % plot magnitude spectra (first averaged across trials envelope in time-domain)
    subplot(3,4,5)
    stem(freq,mXavgTime_cond1,'k','marker','none','linew',linew)
    hold on
    stem(freq(beatFrexIdx),mXavgTime_cond1(beatFrexIdx),'r','marker','none','linew',linew)
    box off
    ax = gca; 
    set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,maxFreqLim],'ylim',mX_ylim, 'fontsize',fontsize)
    ax.XAxis.Label.Color=[0 0 0];
    ax.XAxis.Label.Visible='on';
    title(sprintf('cond1 mX time-domain avg')); 
    xlabel('frequency (Hz)')
    ylabel('magnitude')

    subplot(3,4,6)
    stem(freq,mXavgTime_cond2,'k','marker','none','linew',linew)
    hold on
    stem(freq(beatFrexIdx),mXavgTime_cond2(beatFrexIdx),'r','marker','none','linew',linew)
    box off
    ax = gca; 
    ax.XAxis.Visible = 'off'; 
    set(ax,'xtick',[],'ytick',[], 'xlim',[0,maxFreqLim],'ylim',mX_ylim, 'fontsize',fontsize)
    title(sprintf('cond2 mX time-domain avg')); 
    xlabel('frequency (Hz)')
    ylabel('magnitude')






    % ================================================================================
    % plot magnitude spectra (first FFT taken for each trial, magnitudes averaged across trials)
    subplot(3,4,9)
    stem(freq,mXavgFreq_cond1,'k','marker','none','linew',linew)
    hold on
    stem(freq(beatFrexIdx),mXavgFreq_cond1(beatFrexIdx),'r','marker','none','linew',linew)
    box off
    ax = gca; 
    set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,maxFreqLim],'ylim',mX_ylim, 'fontsize',fontsize)
    ax.XAxis.Label.Color=[0 0 0];
    ax.XAxis.Label.Visible='on';
    title(sprintf('cond1 mX freq-domain avg')); 
    xlabel('frequency (Hz)')
    ylabel('magnitude')

    subplot(3,4,10)
    stem(freq,mXavgFreq_cond2,'k','marker','none','linew',linew)
    hold on
    stem(freq(beatFrexIdx),mXavgFreq_cond2(beatFrexIdx),'r','marker','none','linew',linew)
    box off
    ax = gca; 
    set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,maxFreqLim],'ylim',mX_ylim, 'fontsize',fontsize)
    ax.XAxis.Label.Color=[0 0 0];
    ax.XAxis.Label.Visible='on';
    title(sprintf('cond2 mX freq-domain avg')); 
    xlabel('frequency (Hz)')
    ylabel('magnitude')




    % ================================================================================
    % plot the impulse response
    subplot(3,4,4)
    plot(irT,ir,'color',[255, 153, 0]/255,'linew',linew)
    xlabel('time (s)')
    title('impulse response')
    box off
    set(gca,'xtick',[],'ytick',[],'xlim',[0,irT(end)], 'fontsize',fontsize)




    % ================================================================================
    % plot phase at beat frequency and the mean vector
    subplot(3,4,[7])
    polarplot(aXbeat_cond1, repmat(1,1,length(aXbeat_cond1))', 'ko')
    hold on
    polarplot([meanAngle_cond1(beatFrexIdx),meanAngle_cond1(beatFrexIdx)]', ...
              [0,itpc_cond1(beatFrexIdx)]', 'r','linew',linew)
    set(gca,'thetaAxisUnits','radians','thetatick',[0,pi/2,pi,3*pi/2],'thetaticklabel',{},'rtick',[],'rlim',[0,1.1], 'fontsize',fontsize)
    title(sprintf('cond1')); 

    subplot(3,4,[8])
    polarplot(aXbeat_cond2, repmat(1,1,length(aXbeat_cond2))', 'ko')
    hold on
    polarplot([meanAngle_cond2(beatFrexIdx),meanAngle_cond2(beatFrexIdx)]', ...
              [0,itpc_cond2(beatFrexIdx)]', 'r','linew',linew)
    set(gca,'thetaAxisUnits','radians','thetatick',[0,pi/2,pi,3*pi/2],'thetaticklabel',{},'rtick',[],'rlim',[0,1.1], 'fontsize',fontsize)
    title(sprintf('cond2')); 






    % ================================================================================
    % plot ITPC 
    subplot(3,4,11)
    stem(freq,itpc_cond1,'k','marker','none','linew',linew)
    hold on
    stem(freq(beatFrexIdx),itpc_cond1(beatFrexIdx),'r','marker','none','linew',linew)
    box off
    set(gca,'xtick',[],'ytick',1,'ylim',[0,1.1], 'fontsize',fontsize)
    title(sprintf('cond1 ITPC=%.2f',itpc_cond1(beatFrexIdx))); 
    xlabel('frequency (Hz)')
    ylabel('ITPC')

    subplot(3,4,12)
    stem(freq,itpc_cond2,'k','marker','none','linew',linew)
    hold on
    stem(freq(beatFrexIdx),itpc_cond2(beatFrexIdx),'r','marker','none','linew',linew)
    box off
    set(gca,'xtick',[],'ytick',1,'ylim',[0,1.1], 'fontsize',fontsize)
    title(sprintf('cond2 ITPC=%.2f',itpc_cond2(beatFrexIdx))); 
    xlabel('frequency (Hz)')
    ylabel('ITPC')


    


    %% save the figure

    if DO_SAVE_FIG

        savepath = 'eeg_2conditions'; 
        if ~isdir(savepath); mkdir(savepath); end
        saveas(f, fullfile(savepath,'output.fig'))

        % save also README file with some parameters 
        fid = fopen(fullfile(savepath,'README.txt'),'w'); 
        fprintf(fid, '\n%d sounds in %d events\n',nSounds, nEvents); 
        fprintf(fid, '\n%d sounds on-beat\n%d sounds off-beat\n',nSoundsOnbeat, nSoundsOffbeat); 
        fprintf(fid, '\nbeat interval = %d gridpoints (gridIOI = %dms)\n\n',beatPeriod, round(gridIOI*1000)); 
        fclose(fid); 

    end

       
    
    
end




% ----------------------------------------------------------------------------
% This is a simulation of one SyncSweep experiment with all segments included
% in the plotted summary figures. 
% We can see how stimulus degradation interacts with SNR and neural emphasis, 
% and how it influences the different measures
%     - time-domain average across trials -> FFT 
%     - FFT for each trial -> average magnitudes
% This can be run for avergaged meter-related frequencies, but also
% individual frequencies. 
% ----------------------------------------------------------------------------

clc
clear 

addpath(genpath('lib'))

figSavePath = fullfile('..','figures'); 
if ~isdir(figSavePath)
    mkdir(figSavePath); 
end

colNoEmph           = [102,166,30]/255; 
colConsistentEmph 	= [255, 153, 0]/255; % [102, 153, 0]/255;
colRandomEmph       = [102, 0, 204]/255; 

% construct class
sub = subject(); 

% frequencies of interest
sub.meterFrex       = 1/2.4 * [3,6,12]; % meter-related frequencies 
sub.nonMeterFrex    = []; % meter-unrelated frequencies 

% load sequences used in SyncSweep experiment 
load('../stimuli/XPSyncSweep_allSeqEEG.mat');

% remove 2 padding tones at begining and end
sequences = sequences(:,3:end-2); 

% segmentation parameters
nEventsSeq  = size(sequences,2); 
overlap     = 0.5; 
segmL       = 6*12; 
jumpL       = segmL * overlap; 
nSegm       = (nEventsSeq-(1-overlap)*segmL) / jumpL; 

% simulation parameters
nSub            = 20; % number of subjects in the simulated experiment 
emphStrength    = [0.05, 0.1, 0.3]; % strength of neural emphasis
SNRs            = [Inf]; % signal-to-noise ratio (set to Inf to get 0 noise)


for snri=1:length(SNRs)
    
    fprintf('-----------------------\nSNR = %.2f\n\n',SNRs(snri)); 
    
    sub.noiseSNR = SNRs(snri); 
    
    % figure
    figName = sprintf('%d subjects, SNR %.2f, meter frequencies %s', ...
                      nSub, ...
                      sub.noiseSNR, ...
                      num2str(sub.meterFrex,'%.2fHz ')); 
                  
    f = figure('color','white','position',[286 223 835 682], 'name',figName); 
    pnl = panel(f); 
    pnl.pack('h',3); 
    pnl(1).pack('v',3); 
    pnl(2).pack('v',3); 
    pnl(3).pack('v',3); 

    pnl.select('all'); 

    pnl.fontsize = 12; 

    for emphi=1:length(emphStrength)

        fprintf('emphasis strength = %.2f\n\nSub: ',emphStrength(emphi)); 
        
        % set current on-beat emphasis
        sub.onbeatBrainEmphasis = emphStrength(emphi); 

        % allocate
        ampMeter_timeAvg_noEmphasis         = nan(nSub,nSegm); 
        ampMeter_timeAvg_consistentEmphasis = nan(nSub,nSegm); 
        ampMeter_timeAvg_randomEmphasis     = nan(nSub,nSegm); 

        ampMeter_freqAvg_noEmphasis         = nan(nSub,nSegm); 
        ampMeter_freqAvg_consistentEmphasis = nan(nSub,nSegm); 
        ampMeter_freqAvg_randomEmphasis     = nan(nSub,nSegm); 

        itpcMeter_noEmphasis         = nan(nSub,nSegm); 
        itpcMeter_consistentEmphasis = nan(nSub,nSegm); 
        itpcMeter_randomEmphasis     = nan(nSub,nSegm); 

        % go over subjects
        for subi=1:nSub
            
            fprintf('%d,',subi); 
            
            % go over segments
            for segmi=1:nSegm

                % find segment start
                idx = (segmi-1)*jumpL; 

                % get segment and assign pattern
                sub.patterns = sequences(:,idx+1:idx+segmL); 

                % simulate EEG 
                sub = sub.getEEG(); 

                % analyze EEG 
                sub = sub.analyzeEEG(); 

                ampMeter_timeAvg_noEmphasis(subi,segmi)          = mean(sub.eeg_noEmph.mXavgTime(sub.meterFrexIdx)); 
                ampMeter_timeAvg_consistentEmphasis(subi,segmi)  = mean(sub.eeg_consistentEmph.mXavgTime(sub.meterFrexIdx)); 
                ampMeter_timeAvg_randomEmphasis(subi,segmi)      = mean(sub.eeg_randomEmph.mXavgTime(sub.meterFrexIdx)); 

                ampMeter_freqAvg_noEmphasis(subi,segmi)          = mean(sub.eeg_noEmph.mXavgFreq(sub.meterFrexIdx)); 
                ampMeter_freqAvg_consistentEmphasis(subi,segmi)  = mean(sub.eeg_consistentEmph.mXavgFreq(sub.meterFrexIdx)); 
                ampMeter_freqAvg_randomEmphasis(subi,segmi)      = mean(sub.eeg_randomEmph.mXavgFreq(sub.meterFrexIdx)); 

                itpcMeter_noEmphasis(subi,segmi)                 = mean(sub.eeg_noEmph.itpc(sub.meterFrexIdx)); 
                itpcMeter_consistentEmphasis(subi,segmi)         = mean(sub.eeg_consistentEmph.itpc(sub.meterFrexIdx)); 
                itpcMeter_randomEmphasis(subi,segmi)             = mean(sub.eeg_randomEmph.itpc(sub.meterFrexIdx)); 

            end

        end
        fprintf('\n\n'); 
        
        % PLOT
        pnl(1,emphi).select(); 
        
        plotAcrossSegm(nSegm, ...
                       mean(ampMeter_timeAvg_noEmphasis,1), ...
                       std(ampMeter_timeAvg_noEmphasis,[],1) ./ sqrt(nSub), ...
                       colNoEmph,1.7); 
        plotAcrossSegm(nSegm, ...
                       mean(ampMeter_timeAvg_consistentEmphasis,1), ...
                       std(ampMeter_timeAvg_consistentEmphasis,[],1) ./ sqrt(nSub), ...
                       colConsistentEmph,1.7); 
        plotAcrossSegm(nSegm, ...
                       mean(ampMeter_timeAvg_randomEmphasis,1), ...
                       std(ampMeter_timeAvg_randomEmphasis,[],1) ./ sqrt(nSub), ...
                       colRandomEmph,1.7); 
        set(gca,'xlim',[0,nSegm+1],'box','off','ytick',[],'xtick',[]); 

        
        
        pnl(2,emphi).select(); 
        plotAcrossSegm(nSegm, ...
                       mean(ampMeter_freqAvg_noEmphasis,1), ...
                       std(ampMeter_freqAvg_noEmphasis,[],1) ./ sqrt(nSub), ...
                       colNoEmph,1.7); 
        plotAcrossSegm(nSegm, ...
                       mean(ampMeter_freqAvg_consistentEmphasis,1), ...
                       std(ampMeter_freqAvg_consistentEmphasis,[],1) ./ sqrt(nSub), ...
                       colConsistentEmph,1.7); 
        plotAcrossSegm(nSegm, ...
                       mean(ampMeter_freqAvg_randomEmphasis,1), ...
                       std(ampMeter_freqAvg_randomEmphasis,[],1) ./ sqrt(nSub), ...
                       colRandomEmph,1.7); 
        set(gca,'xlim',[0,nSegm+1],'box','off','ytick',[],'xtick',[]); 

        
        pnl(3,emphi).select(); 
        lineHandle = plotAcrossSegm(nSegm, ...
                       mean(itpcMeter_noEmphasis,1), ...
                       std(itpcMeter_noEmphasis,[],1) ./ sqrt(nSub), ...
                       colNoEmph,1.7); 
        lineHandle(2) = plotAcrossSegm(nSegm, ...
                       mean(itpcMeter_consistentEmphasis,1), ...
                       std(itpcMeter_consistentEmphasis,[],1) ./ sqrt(nSub), ...
                       colConsistentEmph,1.7); 
        lineHandle(3) = plotAcrossSegm(nSegm, ...
                       mean(itpcMeter_randomEmphasis,1), ...
                       std(itpcMeter_randomEmphasis,[],1) ./ sqrt(nSub), ...
                       colRandomEmph,1.7); 
        set(gca,'xlim',[0,nSegm+1],'xtick',[],'box','off','ylim',[0,1],'ytick',[0,1]); 
        
        % plot legend in the corner
        if emphi==1
            l = legend(lineHandle, {'no','consistent','random'}); 
            l.Title.String = 'emphasis'; 
            l.Box = 'off'; 
        end

        pnl(2,emphi).title(sprintf('on-beat emphasis %.2f',sub.onbeatBrainEmphasis)); 

    end

    % labs and titles
    tit1 = pnl(1).title('FFT time average'); 
    tit2 = pnl(2).title('FFT frequency average'); 
    tit3 = pnl(3).title('ITPC'); 

    xlab = pnl.xlabel('segment number (regular to degraded)'); 

    ylab1 = pnl(1,2).ylabel('mean magnitude (a.u.)'); 
    ylab2 = pnl(2,2).ylabel('mean magnitude (a.u.)'); 
    ylab3 = pnl(3,2).ylabel('mean ITPC'); 
    ylab3.Position(1) = ylab1.Position(1); 

    set(pnl(1,3).axis, 'xtick',[1:7]); 
    set(pnl(2,3).axis, 'xtick',[1:7]); 
    set(pnl(3,3).axis, 'xtick',[1:7]); 


    % margins
    pnl.de.margin = 2; 

    pnl.de.marginleft = 6; 

    pnl.de.margintop = 15; 

    pnl.margin = [10, 20, 25, 15]; 


    % tweak text positions (this needs to be done after margins are set
    % otherwise it doesn't work!)
    l.Position = [0.8713    0.8827    0.1293    0.1166]; 

    xlab.Position(2) = -0.05; 

    tit1.FontWeight = 'bold'; 
    tit1.Position(2) = 1.04; 
    tit2.FontWeight = 'bold'; 
    tit2.Position(2) = 1.04; 
    tit3.FontWeight = 'bold'; 
    tit3.Position(2) = 1.04; 


    % save filename
    outName = sprintf('meterFrex%sSNR%.2f_subtrDFTbins%d-%d_nSub%d', ...
                      num2str(sub.meterFrex,'%.2f_'), ...
                      sub.noiseSNR, ...
                      sub.snrmin,sub.snrmax, ...
                      nSub); 
                  
    if ~isdir(figSavePath)
        mkdir(figSavePath); 
    end
    
    % save figure
    saveas(f,fullfile(figSavePath,[outName,'.fig']))
    
    % save subject class instance
    save(fullfile(figSavePath,[outName,'.mat']), 'sub')
    
    
    % close figure
    close(f); 
    
end


%% plot main firugre for example segments 

for segmi=[1,7]
    
    % find segment start
    idx = (segmi-1)*jumpL; 

    % get segment and assign pattern
    sub.patterns = sequences(:,idx+1:idx+segmL); 

    % simulate EEG 
    sub = sub.getEEG(); 

    % analyze EEG 
    sub = sub.analyzeEEG(); 

    f = sub.plotMainFigure();

    outName = sprintf('../figures/mainFig_SNR%.2f_subtrDFTbins%d-%d_segm%d.png', ...
                      sub.noiseSNR, ...
                      sub.snrmin,sub.snrmax, segmi);

    saveas(f,outName); 
    
end


%% plot example trial for example segments 

for segmi=[1,7]
    
    % find segment start
    idx = (segmi-1)*jumpL; 

    % get segment and assign pattern
    sub.patterns = sequences(:,idx+1:idx+segmL); 

    % simulate EEG 
    sub = sub.getEEG(); 

    % analyze EEG 
    sub = sub.analyzeEEG(); 

    f = sub.plotExampleTrial();

    outName = sprintf('../figures/exampleTrial_SNR%.2f_subtrDFTbins%d-%d_segm%d.png', ...
                      sub.noiseSNR, ...
                      sub.snrmin,sub.snrmax, segmi);

    saveas(f,outName); 

end


%% plot impulse response 

f = sub.plotIR; 
saveas(f,'../figures/IR.png'); 









%% plot sequence (pattern) analysis
% 
% ----------------------------------------------------------------------------
% analyse the grid-representation of the sequences using PE, LHL, CHIU-FFT
% etc. 
% ----------------------------------------------------------------------------

% construct class
sub = subject(); 

% set some parameters
sub.meterFrex       = [1.25, 2.5]; 
sub.nonMeterFrex    = []; 

% load sequences used in SyncSweep
load('../stimuli/XPSyncSweep_allSeqEEG.mat');

% remove 2 padding tones at begining and end
sequences = sequences(:,3:end-2); 

% segmentation parameters
nTrials     = 15; 
nEventsSeq  = size(sequences,2); 
overlap     = 0.5; 
segmL       = 6*12; 
jumpL       = segmL * overlap; 
nSegm       = (nEventsSeq-(1-overlap)*segmL) / jumpL; 


PEmin = nan(nTrials,nSegm); 
PErange = nan(nTrials,nSegm); 

LHLmin = nan(nTrials,nSegm); 
LHLrange = nan(nTrials,nSegm); 

CHIUmagn = nan(nTrials,nSegm); 

% go over segments
for segmi=1:nSegm

    % find segment start
    idx = (segmi-1)*jumpL; 

    % get segment and assign pattern
    sub.patterns = sequences(:,idx+1:idx+segmL); 

    % analyze sequences 
    sub = sub.analyzePatterns(); 

    PEmin(:,segmi) = sub.PE.min; 
    PErange(:,segmi) = sub.PE.range; 
    LHLmin(:,segmi) = sub.LHL.min; 
    LHLrange(:,segmi) = sub.LHL.range; 
    
    CHIUmagn(:,segmi) = mean(sub.chiu.meterAmp,2); 
    
end



f = figure('color','white','position',[1150 386 485 549]); 
pnl = panel(f); 

pnl.pack('v',3); 
pnl(1).pack('h',2); 
pnl(2).pack('h',2); 
pnl(3).pack('h',2); 

pnl.fontsize = 16; 

markersize = 10;
col = 'b'; 


pnl(1,1).select(); 
plot([1:nSegm],CHIUmagn,'o','color',col)
hold on
plot([1:nSegm],mean(CHIUmagn,1),'o','color',col,'MarkerFaceColor',col,'markersize',markersize); 
set(gca,'xlim',[0,nSegm+1],'box','off','ytick',[],'xtick',[]); 


pnl(2,1).select(); 
plot([1:nSegm],PEmin,'o','color',col)
hold on
plot([1:nSegm],mean(PEmin,1),'o','color',col,'MarkerFaceColor',col,'markersize',markersize); 
set(gca,'xlim',[0,nSegm+1],'box','off','ytick',[],'xtick',[]); 

pnl(2,2).select(); 
plot([1:nSegm],PErange,'o','color',col)
hold on
plot([1:nSegm],mean(PErange,1),'o','color',col,'MarkerFaceColor',col,'markersize',markersize); 
set(gca,'xlim',[0,nSegm+1],'box','off','ytick',[],'xtick',[]); 


pnl(3,1).select(); 
plot([1:nSegm],LHLmin,'o','color',col)
hold on
plot([1:nSegm],mean(LHLmin,1),'o','color',col,'MarkerFaceColor',col,'markersize',markersize); 
set(gca,'xlim',[0,nSegm+1],'box','off','ytick',[],'xtick',[]); 

pnl(3,2).select(); 
plot([1:nSegm],LHLrange,'o','color',col)
hold on
plot([1:nSegm],mean(LHLrange,1),'o','color',col,'MarkerFaceColor',col,'markersize',markersize); 
set(gca,'xlim',[0,nSegm+1],'box','off','ytick',[],'xtick',[]); 


pnl.de.margin = 2; 
pnl.de.margintop = 10; 

pnl.margin = [10, 15, 5, 10]; 

pnl(1,1).title('Chiu DFT'); 
pnl(2,1).title('P&E minimum'); 
pnl(2,2).title('P&E range'); 
pnl(3,1).title('LHL minimum'); 
pnl(3,2).title('LHL range'); 

pnl(1,1).ylabel('mean meter amplitude'); 
pnl(2,1).ylabel('P&E score'); 
pnl(3,1).ylabel('LHL score'); 

set(pnl(3,1).axis, 'xtick',[1:7]); 
set(pnl(3,2).axis, 'xtick',[1:7]); 

xlab = pnl.xlabel('segment'); 
xlab.Position(2) = -0.05; 


% save figure
saveas(f,fullfile(figSavePath,['sequenceAnalysis.fig']))

% close figure
close(f); 


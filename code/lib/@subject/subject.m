classdef subject
    
    
    
    properties (Access = public)
        % -----------------------------------------------
        
        % beat period (specificed as number of smallest grid intervals)
        beatPeriod = 2; 

        % number of sound events that will be on/off-the-beat
        nSoundsOnbeat; 
        
        % number of sounds in one trial
        nSounds; 

        % total number of events in one trial
        nEvents; 

        % maximum allowed successive silence events
        maxEventsSilence; 

        % the strength of response that is added to on-beat positions in the EEG
        % response (can be interpreted as magnitude of neural process that
        % emphasizes on-beat timepoints)
        onbeatBrainEmphasis = 0.1; 

        % signal to noise ratio for the added white noise 
        % (no noise when set to Inf)
        noiseSNR = Inf;  

        % DFT bins used for noise subtraction in magnitude spectra
        snrmin = 2; 
        snrmax = 3; 

        % -----------------------------------------------
        
        % smallest grid interval (one event duration in seconds)
        gridIOI = 0.2; 

        % onset and offset ramp duration for sound events
        rampon = 0.010; 
        rampoff = 0.010; 

        % number of simulated trials
        nTrials = 30; 
        
        % sampling rate
        fs = 2000; 
        
        % number of samples
        N; 
        
        % trial duration 
        durTrial; 
        
        % time vector
        t; 
        
        % fequency vector
        freq; 
        
        % upper frequency limit for plotting and saving
        maxFreqLim = 6; 
        
        % meter-related frequencies
        meterFrex; 
        
        % meter-related frequency index
        meterFrexIdx; 

        % meter-unrelated frequencies
        nonMeterFrex; 
        
        % meter-unrelated frequency index
        nonMeterFrexIdx; 

        % Povel-Essens grouping
        PEgrouping = 4; 
        
        % LHL meter string
        LHLmeter = '2_4'; 
        
        
        % immpulse response
        ir;
        % time vector or IR
        irT; 
        % impulse response amplitude (apply each entry to the correcponding
        % F0 partial)
        irA = 1; 
        % impulse response latency
        irT0 = 0; 
        % impulse response time constant
        irTau = 0.050; 
        % impulse response frequency components (can be vector of
        % frequencies)
        irF0 = 4; 
        % impulse response duration
        irDur = 0.5; 
        
                
        % matrix of patterns (grid representation) 
        patterns; 
        
        stim_noEmph; 
        stim_consistentEmph; 
        stim_randomEmph; 
        
        eeg_noEmph; 
        eeg_consistentEmph; 
        eeg_randomEmph; 
        
        chiu; 
        PE; 
        LHL; 
    end
    
    
    
    
    properties (Access = private)

        % maximum frequency index
        maxFreqIdx; 
        
        % -----------------------------------------------
        % PLOTING PROPERTIES
        colNoEmph = [102,166,30]/255; 
        colConsistentEmph = [255, 153, 0]/255; % [102, 153, 0]/255;
        colRandomEmph   = [102, 0, 204]/255; 

        % impulse response color
        irCol = [217,217,217]/255; 
        
        % get separate color for each meter frequency
        % obj.meterCols = linspecer(length(obj.meterFrex));             
        meterCols    = [228,26,28
                        55,126,184
                        77,175,74
                        152,78,163
                        255,127,0
                        255,255,51
                        166,86,40]/255; 

        TD_ylim = 1; 
        mXavgTime_ylim = [-0.1,1]; 
        mXavgFreq_ylim = [-0.1,1]; 

        fontsize = 16; 
        linew = 1.7; 

    end
    
    
    
%%  METHODS    



    methods
            
    
        function obj = subject()
        % class constructor 
            if nargin>0
                
            end
            
        end

        
        
        
        
        function obj = getIR(obj)
            
            % IR for envelope-following process

            % (ERP-like kernel)
            irN         = round(obj.irDur*obj.fs); 
            obj.irT     = [0:irN-1]/obj.fs; 
            irAmp       = (obj.irT-obj.irT0)/obj.irTau .* exp(1-(obj.irT-obj.irT0)/obj.irTau); 
            obj.ir      = zeros(1,irN); 
            for fi=1:length(obj.irF0)
                obj.ir = obj.ir + obj.irA(fi) * sin( 2*pi*obj.irF0(fi)*(obj.irT-obj.irT0) ) ; 
            end
            obj.ir = obj.ir ./ length(obj.irF0); 
            obj.ir = irAmp .* obj.ir; 


%             % (ERP-like kernel with fast+slow component)
%             obj.irF0    = [1,7]; 
%             obj.irA     = [1,2]; 
%             obj.irT0    = [0,0];  
%             obj.irTau   = [0.200, 0.050]; 
%             obj.irDur   = 1; 
%             irN         = round(obj.irDur*obj.fs); 
%             obj.irT     = [0:irN-1]/obj.fs; 
%             obj.ir      = zeros(1,irN); 
%             for fi=1:length(obj.irF0)
%                 obj.ir = obj.ir + obj.irA(fi) * ...
%                          (obj.irT-obj.irT0(fi))/obj.irTau(fi) .* ...
%                          exp(1-(obj.irT-obj.irT0(fi))/obj.irTau(fi)) .* ...
%                          sin( 2*pi*obj.irF0(fi)*(obj.irT-obj.irT0(fi)) ) ; 
%             end
%             obj.ir = obj.ir ./ length(obj.irF0); 

   
%             % half gaussian
%             obj.irA = 1; 
%             obj.irT0 = 0; 
%             obj.irTau = 0.100; 
%             obj.irDur = 1; 
%             irN = round(obj.irDur*obj.fs); 
%             obj.irT = [0:irN-1]/obj.fs; 
%             obj.ir = obj.irA * exp(-(obj.irT-obj.irT0).^2/(2*obj.irTau^2)) ; 

%             % (sound-envelope kernel)
%             obj.irT = [0:round(obj.gridIOI*obj.fs)-1]/obj.fs; 
%             obj.ir = ones(1, length(obj.irT)); 
%             obj.ir(1:round(obj.rampon*obj.fs)) = linspace(0,1,round(obj.rampon*obj.fs)); 
%             obj.ir(end-round(obj.rampoff*obj.fs)+1:end) = linspace(1,0,round(obj.rampoff*obj.fs)); 

            
            % normalize IR to have RMS of 1
            obj.ir = obj.ir./rms(obj.ir); 
                        
        end
        
               
        
        
        
        function obj = getEEG(obj)

            % update number of trials
            obj.nTrials = size(obj.patterns,1); 
            
            % update nEvents
            obj.nEvents = size(obj.patterns,2); 
            
            % get trial duration and number of samples
            obj.durTrial = obj.nEvents*obj.gridIOI; 
            obj.N = round(obj.durTrial*obj.fs); 

            % make time vector for the trial
            obj.t = [0 : obj.N-1]/obj.fs; 

            % get the impulse response from another method
            obj = getIR(obj); 
            
            % allocate 
            stim_noEmph = zeros(obj.nTrials,obj.N); 
            stim_consistentEmph = zeros(obj.nTrials,obj.N); 
            stim_randomEmph = zeros(obj.nTrials,obj.N); 

            eeg_noEmph = zeros(obj.nTrials,obj.N); 
            eeg_consistentEmph = zeros(obj.nTrials,obj.N); 
            eeg_randomEmph = zeros(obj.nTrials,obj.N); 

            % go over trials
            for triali=1:obj.nTrials

                soundPos = find(obj.patterns(triali,:)); 

                % construct stimulus representation using impulse functions
                soundPositionsIdx = round((soundPos-1)*obj.gridIOI*obj.fs)+1; 

                % ------------------ eeg_noEmph ------------------
                stim_noEmph(triali, soundPositionsIdx) = 1; 

                % convolve with event envelope to get "no-transformation EEG
                % response"
                eeg_noEmph(triali,:) = stim_noEmph(triali,:); 
                tmp = conv( obj.ir, eeg_noEmph(triali,:), 'full'); 
                eeg_noEmph(triali,:) = tmp(1:size(eeg_noEmph,2)); 

                % -------------- eeg_consistentEmph --------------
                stim_consistentEmph(triali,:) = stim_noEmph(triali,:); 

                % make the beat phase identical (consistent across trials)
                beatPhase = 1; 
                onbeatPosCond1 = [beatPhase:obj.beatPeriod:obj.nEvents]; 
                onbeatPosIdx = round((onbeatPosCond1-1)*obj.gridIOI*obj.fs)+1; 
                stim_consistentEmph(triali,onbeatPosIdx) = stim_consistentEmph(triali,onbeatPosIdx) + obj.onbeatBrainEmphasis; 

                % convolve with event envelope
                tmp = conv( obj.ir, stim_consistentEmph(triali,:), 'full'); 
                eeg_consistentEmph(triali,:) = tmp(1:size(eeg_consistentEmph,2)); 

                % ----------------- eeg_randomEmph -----------------
                stim_randomEmph(triali,:) = stim_noEmph(triali,:); 

                % choose the beat phase randomly across trials
                beatPhase = randsample([1:obj.beatPeriod],1); 
                onbeatPosCond2 = [beatPhase:obj.beatPeriod:obj.nEvents]; 
                onbeatPosIdx = round((onbeatPosCond2-1)*obj.gridIOI*obj.fs)+1; 
                stim_randomEmph(triali,onbeatPosIdx) = stim_randomEmph(triali,onbeatPosIdx) + obj.onbeatBrainEmphasis; 

                % convolve with event envelope
                tmp = conv( obj.ir, stim_randomEmph(triali,:), 'full'); 
                eeg_randomEmph(triali,:) = tmp(1:size(eeg_randomEmph,2)); 
            end

            % ===============================================
            % add pink noise

            % baseline (no transformation)
            eeg_noEmph = addScaledNoise(eeg_noEmph,obj.noiseSNR,obj.N,obj.fs); 

            % condition 1 (consistent beat alignement)
            eeg_consistentEmph = addScaledNoise(eeg_consistentEmph,obj.noiseSNR,obj.N,obj.fs); 

            % condition 2 (random beat alignement across trials)
            eeg_randomEmph = addScaledNoise(eeg_randomEmph,obj.noiseSNR,obj.N,obj.fs); 

            
            % ===============================================
            % assign
            
            obj.stim_noEmph = []; 
            obj.stim_noEmph.x = stim_noEmph; 
            
            obj.stim_consistentEmph = []; 
            obj.stim_consistentEmph.x = stim_consistentEmph; 

            obj.stim_randomEmph = []; 
            obj.stim_randomEmph.x = stim_randomEmph; 
            
            
            obj.eeg_noEmph = []; 
            obj.eeg_noEmph.x = eeg_noEmph; 
            
            obj.eeg_consistentEmph = []; 
            obj.eeg_consistentEmph.x = eeg_consistentEmph; 

            obj.eeg_randomEmph = []; 
            obj.eeg_randomEmph.x = eeg_randomEmph; 
            
            
            % update ylims for plotting time-domain 
            obj.TD_ylim = max(abs([mean(obj.eeg_noEmph.x,1), ...
                               mean(obj.eeg_consistentEmph.x,1), ...
                               mean(obj.eeg_randomEmph.x,1)])); 

        end        
        
        
        
        
        
        function obj = analyzeEEG(obj)        
            % analyse data in all the conditions usign DFT
            
            % find maximum frequency index (for plotting and saving)
            obj.maxFreqIdx = round(obj.maxFreqLim/obj.fs*obj.N)+1; 

            % make vector of frequency bins in Hz
            obj.freq = [0:obj.maxFreqIdx-1]/obj.N*obj.fs; 

            % find meter-related frequency index
            obj.meterFrexIdx = obj.meterFrex ./ obj.fs * obj.N + 1; 
            if any(mod(obj.meterFrexIdx,1))
                warning(['Meter-related frequencies [',sprintf('%.2f, ',obj.meterFrex(logical(mod(obj.meterFrexIdx,1)))),'Hz] are not centered on an FFT bin! I will round them to the closest bin...']); 
                obj.meterFrexIdx = round(obj.meterFrexIdx); 
            end

            % find meter-unrelated frequency index
            obj.nonMeterFrexIdx = obj.nonMeterFrex ./ obj.fs * obj.N + 1; 
            if any(mod(obj.nonMeterFrexIdx,1))
                warning(['Meter-unrelated frequencies [',sprintf('%.2f, ',obj.nonMeterFrex(logical(mod(obj.nonMeterFrexIdx,1)))),'Hz] are not centered on an FFT bin! I will round them to the closest bin...']); 
                obj.nonMeterFrexIdx = round(obj.nonMeterFrexIdx); 
            end
            
            % analyse
            obj.eeg_noEmph          = analyseFFT(obj, obj.eeg_noEmph); 
            obj.eeg_consistentEmph  = analyseFFT(obj, obj.eeg_consistentEmph); 
            obj.eeg_randomEmph      = analyseFFT(obj, obj.eeg_randomEmph); 
               
            
            % update mX_ylim for plotting
            obj.mXavgTime_ylim = max([obj.eeg_noEmph.mXavgTime(obj.snrmax+1:end), ...
                                      obj.eeg_consistentEmph.mXavgTime(obj.snrmax+1:end), ...
                                      obj.eeg_randomEmph.mXavgTime(obj.snrmax+1:end)]); 
                                      
            obj.mXavgFreq_ylim = max([obj.eeg_noEmph.mXavgFreq(obj.snrmax+1:end), ...
                                      obj.eeg_consistentEmph.mXavgFreq(obj.snrmax+1:end), ...
                                      obj.eeg_randomEmph.mXavgFreq(obj.snrmax+1:end)]); 
        end
        
        
        

        
        
        
        function data = analyseFFT(obj, data)
            % helper function for analyzeEEG that carries out FFT analysis on
            % the passed data
            % ----------------------------------------------------
            % PHASE 

            % get complex-valued FFT
            data.X = fft(data.x,[],2) / obj.N * 2; 
            data.X = data.X(:,1:obj.maxFreqIdx); 

            % get angles
            data.aX = angle(data.X); 

            % get angles at meter frequencies
            data.aXmeter = data.aX(:,obj.meterFrexIdx); 

            % get mean vector, it's angle and length (aka itpc)
            data.meanVector  = mean(exp(1j*data.aX),1); 
            data.meanAngle   = angle(data.meanVector); 
            data.itpc        = abs(data.meanVector); 

            % ----------------------------------------------------
            % MAGNITUDE

            % average trials in frequency domain (i.e. average magnitudes)
            data.mXavgFreq = mean(abs(data.X),1); 
            % set DC to 0
            data.mXavgFreq(1) = 0; 
            % get empirical SNR from the magnitude spectra
            data.snrAvgFreq = nan(1,length(obj.meterFrex)); 
            for fi=1:length(obj.meterFrex)
                % get surrounding bin indices
                signalIdx       = obj.meterFrexIdx(fi); 
                noiseIdx        = [ signalIdx-obj.snrmax : signalIdx-obj.snrmin, signalIdx+obj.snrmin : signalIdx+obj.snrmax ]; 
                signalAmp       = data.mXavgFreq(signalIdx); 
                noiseAmp        = mean(data.mXavgFreq(noiseIdx)); 
                data.snrAvgFreq(fi)  = signalAmp / noiseAmp; 
            end
            % subtract surrounding bins
            data.mXavgFreq = subtractSNR(data.mXavgFreq,obj.snrmin,obj.snrmax); 


            
            % average trials in the time domain (i.e. average complex numbers)
            data.mXavgTime = abs(fft(mean(data.x,1)) / obj.N * 2); 
            data.mXavgTime = data.mXavgTime(:,1:obj.maxFreqIdx); 
            % set DC to 0
            data.mXavgTime(1) = 0;             
            % get empirical SNR from the magnitude spectra
            data.snrAvgTime = nan(1,length(obj.meterFrex)); 
            for fi=1:length(obj.meterFrex)
                % get surrounding bin indices
                signalIdx       = obj.meterFrexIdx(fi); 
                noiseIdx        = [ signalIdx-obj.snrmax : signalIdx-obj.snrmin, signalIdx+obj.snrmin : signalIdx+obj.snrmax ]; 
                signalAmp       = data.mXavgTime(signalIdx); 
                noiseAmp        = mean(data.mXavgTime(noiseIdx)); 
                data.snrAvgTime(fi)  = signalAmp / noiseAmp; 
            end
            % subtract surrounding bins
            data.mXavgTime = subtractSNR(data.mXavgTime,obj.snrmin,obj.snrmax); 

        end

        
        
        
        
        function obj = analyzePatterns(obj)
            
            % update number of trials
            obj.nTrials = size(obj.patterns,1); 
            
            % update nEvents
            obj.nEvents = size(obj.patterns,2); 
            
            % ----------------------------------------------------
            % CHIU DFT
            
            chiu = []; 
            
            chiu.N = size(obj.patterns,2);  
            chiu.hN = floor(chiu.N/2)+1; 
            chiu.X = fft(obj.patterns,[],2); 
            chiu.X = chiu.X(:,1:chiu.hN); 
            
            % meter-related
            chiu.meterEvents = (1 ./ obj.meterFrex) ./ obj.gridIOI;            
            chiu.meterCycles = chiu.N ./ chiu.meterEvents;            
            if any(chiu.meterCycles>=chiu.N)
                warning(sprintf('\n----- CHIU FFT ----- \nsome meterFrex are too high for Chiu-FFT, I will remove them...')); 
                chiu.meterEvents(chiu.meterCycles>=chiu.N) = []; 
                chiu.meterCycles(chiu.meterCycles>=chiu.N) = []; 
            end
            if any(mod(chiu.meterCycles,1))
                warning(sprintf('\n----- CHIU FFT ----- \nsome meterFrex don''t have integer number of cycles in the analysis window, I will remove them...')); 
                chiu.meterEvents(logical(mod(chiu.meterCycles,1))) = []; 
                chiu.meterCycles(logical(mod(chiu.meterCycles,1))) = []; 
            end            
            chiu.meterIdx = chiu.meterCycles+1; 

            % meter-unrelated
            chiu.nonMeterEvents = (1 ./ obj.nonMeterFrex) ./ obj.gridIOI;            
            chiu.nonMeterCycles = chiu.N ./ chiu.nonMeterEvents;            
            if any(chiu.nonMeterCycles>=chiu.N)
                warning(sprintf('\n----- CHIU FFT ----- \nsome nonMeterFrex are too high for Chiu-FFT, I will remove them...')); 
                chiu.nonMeterEvents(chiu.nonMeterCycles>=chiu.N) = []; 
                chiu.nonMeterCycles(chiu.nonMeterCycles>=chiu.N) = []; 
            end
            if any(mod(chiu.nonMeterCycles,1))
                warning(sprintf('\n----- CHIU FFT ----- \nsome nonMeterFrex don''t have integer number of cycles in the analysis window, I will remove them...')); 
                chiu.nonMeterEvents(logical(mod(chiu.nonMeterCycles,1))) = []; 
                chiu.nonMeterCycles(logical(mod(chiu.nonMeterCycles,1))) = []; 
            end            
            chiu.nonMeterIdx = chiu.nonMeterCycles+1; 
            
%             figure
%             b = bar(abs(chiu.X(5,:)),'k'); 
%             hold on
%             plot(chiu.meterIdx,abs(chiu.X(5,chiu.meterIdx)),'ro','MarkerFaceColor','r')
            
            chiu.meterAmp = abs(chiu.X(:,chiu.meterIdx)); 
            chiu.nonMeterAmp = abs(chiu.X(:,chiu.nonMeterIdx)); 
            
            obj.chiu = chiu; 
            
            % ----------------------------------------------------
            % Povel & Essens (PE)

            PE = []; 
            
            % original phase
            PE.origPhase = syncopationPE(obj.patterns,obj.PEgrouping,'zeropad'); 
            
            % move by +-2 events (assume circularity) 
            PE.circShifted = nan(obj.nTrials, 5); 
            for nShift=[-2:2]; 
                PE.circShifted(:,nShift+3) = syncopationPE(circshift(obj.patterns,nShift,2), obj.PEgrouping,'zeropad')';    
            end
            
            PE.range    = max(PE.circShifted,[],2) - min(PE.circShifted,[],2); 
            PE.min      = min(PE.circShifted,[],2); 
            
            obj.PE = PE; 
            % ----------------------------------------------------
            % LHL

            LHL = []; 
            
            % original phase
            LHL.origPhase = syncopationLHL(obj.patterns, obj.LHLmeter); 
            
            % move by +-2 events (assume circularity) 
            LHL.circShifted = nan(obj.nTrials, 5); 
            for nShift=[-2:2]; 
                LHL.circShifted(:,nShift+3) = syncopationLHL(circshift(obj.patterns,nShift,2), obj.LHLmeter);    
            end
            
            LHL.range    = max(LHL.circShifted,[],2) - min(LHL.circShifted,[],2); 
            LHL.min      = min(LHL.circShifted,[],2); 
            
            obj.LHL = LHL; 
        end
        

        
        
        
        
        



        function f = plotExampleTrial(obj)  
        % plot example trial env (transform vs. no-transform)

            % plot envelope from one example trial for EEG without any transformation
            % (one-to-one stimulus tracking) and also for EEG with neural empahsis at
            % on-beat positions
            f = figure('color','white','position',[127 146 861 480]); 
            pnl = panel(f); 
            pnl.pack({40,60}); 
            pnl(1).pack('h',3); 


            pnl(1,1).select()
            plot(obj.t, obj.stim_noEmph.x(1,:),'color','k', 'linew',obj.linew)
            tit = pnl(1,1).title(sprintf('stimulus \nrepresentation')); 
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.t(end)], 'ylim',[0,obj.TD_ylim], 'fontsize',obj.fontsize)

            pnl(1,2).select()
            plot(obj.irT,obj.ir,'color',[255, 153, 0]/255,'linew',obj.linew)
            tit = pnl(1,2).title('impulse response'); 
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.irT(end)], 'fontsize',obj.fontsize)

            pnl(1,3).select()
            plot(obj.t,obj.stim_consistentEmph.x(1,:),'color','r', 'linew',obj.linew)
            tit = pnl(1,3).title(sprintf('stimulus \nrepresentation \n(with emphasis)')); 
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.t(end)], 'ylim',[0,obj.TD_ylim], 'fontsize',obj.fontsize)


            pnl(2).select()
            plot(obj.t, obj.eeg_noEmph.x(1,:), 'color','k', 'linew',obj.linew)
            hold on
            plot(obj.t,obj.eeg_consistentEmph.x(1,:),'color','r', 'linew',obj.linew)
            idx = find(obj.stim_noEmph.x(1,:)); 
            plot([obj.t(idx); obj.t(idx)], [-100;100], ':','color',[145,125,145]/255,'linew',obj.linew); 
            tit = pnl(2).title('example trial'); 
            l=legend({'no transformation','emphasis on beat'}); 
            l.Position = [0.7950 0.4611 0.1823 0.0823]; 
            l.Box = 'off'; 
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.t(end)], 'ylim',[-obj.TD_ylim,obj.TD_ylim], 'fontsize',obj.fontsize)


            xlab = pnl(1).xlabel('Time (s)'); 
            ylab = pnl(1).ylabel('Amplitude (a.u.)'); 

            xlab = pnl(2).xlabel('Time (s)'); 
            ylab = pnl(2).ylabel('Amplitude (a.u.)'); 

            pnl.de.margin = 2; 
            pnl(2).margintop = 25;
            pnl.margin = [13 10 2 15];

            pnl.fontsize = obj.fontsize; 


        end



        function f = plotIR(obj)
            
            % impulse response

            f = figure('color','white','position',[665 480 570 193]); 
            pnl = panel(f); 
            
            pnl.pack('h',2); 
            
            pnl(1).select(); 
            plot(obj.irT,obj.ir,'color',obj.irCol,'linew',obj.linew)
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.irT(end)], 'fontsize',obj.fontsize)
            xlabel('time'); 
            
            pnl(2).select(); 
            mX = abs(fft(obj.ir,obj.N)); 
            mX = mX(1:obj.maxFreqIdx); 
            plot(obj.freq,mX,'linew',obj.linew)
            xlabel('frequency (Hz)'); 
            set(gca,'xtick',[0,obj.maxFreqLim]); 
            set(gca,'ytick',[]); 
            
            pnl.fontsize = obj.fontsize; 

        end
        
        
        
        
        function f = plotMainFigure(obj)
                        
            
            fig_name = sprintf('brain emphasis %.1f | noise SNR %.2f', ...
                                obj.onbeatBrainEmphasis, obj.noiseSNR); 

            f = figure('color','white','position',[84 213 1747 651],'name',fig_name); 
            pnl = panel(f); 
            
            pnl.pack('h',2); 
            
            pnl(1).pack('v',{1/3, 2/3}); 
            pnl(1,1).pack('h',3); 
            pnl(1,2).pack('v',2); 
            pnl(1,2,1).pack('h',3); 
            pnl(1,2,2).pack('h',3); 
            
            pnl(2).pack('v',{1/3, 2/3}); 
            pnl(2,1).pack('h',2); 
            pnl(2,2).pack('h',3); 
            pnl(2,2,1).pack('v',{50,[]}); 
            pnl(2,2,2).pack('v',{50,[]}); 
            pnl(2,2,3).pack('v',{50,[]}); 
                        
            pnl.fontsize = obj.fontsize; 
            pnl.fontweight = 'normal';             
            
            % ----------------------------------------------------
            % plot envelope averaged across trials
            
            % no emphasis
            pnl(1,1,1).select(); 
            plot(obj.t, mean(obj.eeg_noEmph.x,1),'color',obj.colNoEmph, 'linew',obj.linew)
            title('mean across trials')
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.t(end)], 'ylim',[-obj.TD_ylim,obj.TD_ylim], 'fontsize',obj.fontsize)

            % consitent emphasis
            pnl(1,1,2).select(); 
            plot(obj.t, mean(obj.eeg_consistentEmph.x,1),'color',obj.colConsistentEmph, 'linew',obj.linew)
            title('mean across trials')
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.t(end)], 'ylim',[-obj.TD_ylim,obj.TD_ylim], 'fontsize',obj.fontsize)

            % random emphasis
            pnl(1,1,3).select(); 
            plot(obj.t,mean(obj.eeg_randomEmph.x,1),'color',obj.colRandomEmph, 'linew',obj.linew)
            title('mean across trials')
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.t(end)], 'ylim',[-obj.TD_ylim,obj.TD_ylim], 'fontsize',obj.fontsize)
                        
            % overlay
            pnl(2,1,1).select(); 
            hold on
            plot(obj.t, mean(obj.eeg_consistentEmph.x,1),'color',obj.colConsistentEmph, 'linew',obj.linew)
            plot(obj.t,mean(obj.eeg_randomEmph.x,1),'color',obj.colRandomEmph, 'linew',obj.linew)
            l = legend({'consistent','random'}); 
            l.Box = 'off';  
            l.Position = [0.6644 0.9418 0.1139 0.0494]; 
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.t(end)], 'ylim',[-obj.TD_ylim,obj.TD_ylim], 'fontsize',obj.fontsize)
            
            
            % ----------------------------------------------------
            % magnitude spectra (averaged across trials envelope in time-domain)
            
            % no emphasis
            pnl(1,2,1,1).select(); 
            stem(obj.freq, obj.eeg_noEmph.mXavgTime, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_noEmph.mXavgTime(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[-0.1, obj.mXavgTime_ylim], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';

            % consitent emphasis
            pnl(1,2,1,2).select(); 
            stem(obj.freq, obj.eeg_consistentEmph.mXavgTime, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_consistentEmph.mXavgTime(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[-0.1, obj.mXavgTime_ylim], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';
            
            % random emphasis
            pnl(1,2,1,3).select(); 
            stem(obj.freq, obj.eeg_randomEmph.mXavgTime, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_randomEmph.mXavgTime(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[-0.1, obj.mXavgTime_ylim], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';
                        
            
            % ----------------------------------------------------
            % magnitude spectra (first take FFT for each trial, then average magnitudes)
            
            % no emphasis
            pnl(1,2,2,1).select(); 
            stem(obj.freq, obj.eeg_noEmph.mXavgFreq, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_noEmph.mXavgFreq(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[-0.1, obj.mXavgFreq_ylim], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';

            % consitent emphasis
            pnl(1,2,2,2).select(); 
            stem(obj.freq, obj.eeg_consistentEmph.mXavgFreq, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_consistentEmph.mXavgFreq(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[-0.1, obj.mXavgFreq_ylim], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';
            
            % random emphasis
            pnl(1,2,2,3).select(); 
            stem(obj.freq, obj.eeg_randomEmph.mXavgFreq, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_randomEmph.mXavgFreq(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[-0.1, obj.mXavgFreq_ylim], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';
                        

            % ----------------------------------------------------
            % impulse response
            
            pnl(2,1,2).select(); 
            plot(obj.irT,obj.ir,'color',obj.irCol,'linew',obj.linew)
            pnl(2,1,2).title('impulse response')
            box off
            set(gca,'xtick',[],'ytick',[],'xlim',[0,obj.irT(end)], 'fontsize',obj.fontsize)


            % ----------------------------------------------------
            % phase at beat frequency and the mean vector (polar plot)

            % NO emphasis            
            axCart = axes('Color', 'none'); % dummy cartesian axes
            axCart.XAxis.Visible = 'off'; 
            axCart.YAxis.Visible = 'off'; 
            axPolar = polaraxes; 
            clear h
            for fi=1:length(obj.meterFrex)
                polarplot(obj.eeg_noEmph.aXmeter(:,fi), repmat(1,1,obj.nTrials)','o','color',obj.meterCols(fi,:))
                hold on
                h(fi) = polarplot([obj.eeg_noEmph.meanAngle(obj.meterFrexIdx(fi)), ...
                                   obj.eeg_noEmph.meanAngle(obj.meterFrexIdx(fi))]', ...
                                  [0,obj.eeg_noEmph.itpc(obj.meterFrexIdx(fi))]', ...
                                  'color',obj.meterCols(fi,:),'linew',obj.linew);                    
            end            
            set(gca,'thetaAxisUnits','radians','thetatick',[0,pi/2,pi,3*pi/2],'thetaticklabel',{},'rtick',[],'rlim',[0,1.1], 'fontsize',obj.fontsize)
                        
            legendTxt = cellfun(@(x)sprintf('%.2f Hz',x), num2cell(obj.meterFrex),'uni',0); 
            lPol = legend(h,legendTxt); 
            lPol.Box = 'off'; 

            pnl(2,2,1,1).select([axCart,axPolar]); 
            

            % consitent emphasis            
            axCart = axes('Color', 'none'); % dummy cartesian axes
            axCart.XAxis.Visible = 'off'; 
            axCart.YAxis.Visible = 'off'; 
            axPolar = polaraxes; 
            clear h
            for fi=1:length(obj.meterFrex)
                polarplot(obj.eeg_consistentEmph.aXmeter(:,fi), repmat(1,1,obj.nTrials)','o','color',obj.meterCols(fi,:))
                hold on
                h(fi) = polarplot([obj.eeg_consistentEmph.meanAngle(obj.meterFrexIdx(fi)), ...
                                   obj.eeg_consistentEmph.meanAngle(obj.meterFrexIdx(fi))]', ...
                                  [0,obj.eeg_consistentEmph.itpc(obj.meterFrexIdx(fi))]', ...
                                  'color',obj.meterCols(fi,:),'linew',obj.linew);                    
            end            
            set(gca,'thetaAxisUnits','radians','thetatick',[0,pi/2,pi,3*pi/2],'thetaticklabel',{},'rtick',[],'rlim',[0,1.1], 'fontsize',obj.fontsize)
            pnl(2,2,2,1).select([axCart,axPolar]); 
            
            
            % random emphasis
            axCart = axes('Color', 'none'); % dummy cartesian axes
            axCart.XAxis.Visible = 'off'; 
            axCart.YAxis.Visible = 'off'; 
            axPolar = polaraxes; 
            clear h
            for fi=1:length(obj.meterFrex)
                polarplot(obj.eeg_randomEmph.aXmeter(:,fi), repmat(1,1,obj.nTrials)','o','color',obj.meterCols(fi,:))
                hold on
                h(fi) = polarplot([obj.eeg_randomEmph.meanAngle(obj.meterFrexIdx(fi)), ...
                                   obj.eeg_randomEmph.meanAngle(obj.meterFrexIdx(fi))]', ...
                                  [0,obj.eeg_randomEmph.itpc(obj.meterFrexIdx(fi))]', ...
                                  'color',obj.meterCols(fi,:),'linew',obj.linew);                    
            end            
            set(gca,'thetaAxisUnits','radians','thetatick',[0,pi/2,pi,3*pi/2],'thetaticklabel',{},'rtick',[],'rlim',[0,1.1], 'fontsize',obj.fontsize)
            pnl(2,2,3,1).select([axCart,axPolar]); 
            
            
            % ----------------------------------------------------
            %  ITPC 
            
            % no emphasis
            pnl(2,2,1,2).select(); 
            stem(obj.freq, obj.eeg_noEmph.itpc, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_noEmph.itpc(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[0,1], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';
            
            % consitent emphasis
            pnl(2,2,2,2).select(); 
            stem(obj.freq, obj.eeg_consistentEmph.itpc, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_consistentEmph.itpc(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[0,1], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';
            
            % random emphasis
            pnl(2,2,3,2).select(); 
            stem(obj.freq, obj.eeg_randomEmph.itpc, 'k','marker','none','linew',obj.linew)
            hold on
            for fi=1:length(obj.meterFrex)
                stem(obj.freq(obj.meterFrexIdx(fi)),obj.eeg_randomEmph.itpc(obj.meterFrexIdx(fi)),'color',obj.meterCols(fi,:),'marker','none','linew',obj.linew)
            end
            box off
            ax = gca; 
            set(ax,'xcolor','none','xtick',[],'ytick',[], 'xlim',[0,obj.maxFreqLim],'ylim',[0,1], 'fontsize',obj.fontsize)
            ax.XAxis.Label.Color=[0 0 0];
            ax.XAxis.Label.Visible='on';

            
            % ----------------------------------------------------
            %  labels and titles
            pnl(1,1,1).title('no emphasis'); 
            pnl(1,1,2).title('consistent emphasis'); 
            pnl(1,1,3).title('random emphasis'); 
            pnl(1,1).xlabel('time (s)'); 
            pnl(1,1).ylabel('amplitude (a.u.)'); 

            pnl(2,1).xlabel('time (s)'); 
            
            pnl(1,2,1).title(sprintf('time-domain avg')); 

            pnl(1,2,2).title(sprintf('freq-domain avg')); 
            
            pnl(1,2).xlabel('frequency (Hz)')
            pnl(1,2).ylabel('magnitude (a.u.)')
  

            pnl(2,2,1).title('no emphasis')
            pnl(2,2,2).title('consistent emphasis')
            pnl(2,2,3).title('random emphasis')

            pnl(2,2).xlabel('frequency (Hz)')
            pnl(2,2,1,2).ylabel('ITPC')


            % ----------------------------------------------------
            %  set margins (from smallest to largest)
            
            pnl.de.margin = 5; 
            
            pnl(1).de.marginbottom = 10; 

            
            pnl(2,2,1).marginright = 15; 
            pnl(2,2,2).marginright = 15; 
            pnl(2,2,3).marginright = 15; 
            
            pnl(1,1).marginbottom = 20;             
            pnl(2,1).marginbottom = 20; 
            
            pnl(2).marginleft = 15;
            
            pnl.margin = [13 10 2 15];

            
            
            
        end
        

        
        
        
        
        
    end
    
    
    
    
end
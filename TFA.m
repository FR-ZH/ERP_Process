addpath('/Users/Ray/Documents/MATLAB_ToolBoxes/Hu_L')
clc
for sub = 1:73
    tempD1 = FL.ERP_incm{sub};
    tempD2 = FL.ERP_cm{sub};
    for ch = 1:size(tempD1,1)
        disp(['Processing sub: ',num2str(sub), ',  ch: ',num2str(ch)]);  
        Fs = 500;
        x1 = squeeze(tempD1(ch,:,:)); % [Time × Trials]
        x2 = squeeze(tempD2(ch,:,:)); % [Time × Trials]
        xtimes = TimeStamps/1000; % -0.5:0.002:0.998
        t_stft = xtimes; 
        f_stft = 0 : 0.5 : 48; 
        winsize = 0.2; % 200 ms window
        [~,P1] = sub_stft(x1, xtimes, t_stft, f_stft, Fs, winsize);
        [~,P2] = sub_stft(x2, xtimes, t_stft, f_stft, Fs, winsize);       
        baseline_idx = (t_stft >= -0.4 & t_stft <= -0.05);
        baseline_power1 = mean(P1(:, baseline_idx, :), 2);
        baseline_power2 = mean(P2(:, baseline_idx, :), 2);
        P_dB1 = 10*log10(P1 ./ baseline_power1);
        P_dB2 = 10*log10(P2 ./ baseline_power2);
        P_dB_All_incm(sub,ch,:,:)=mean(P_dB1,3);
        P_dB_All_cm(sub,ch,:,:)=mean(P_dB2,3);
        clear trial P_dB1 P1 baseline_power1 x1 P_dB2 P2 baseline_power2 x2
    end
    clear tempD1 tempD2
end
%%
addpath('/Users/Ray/Documents/MATLAB_ToolBoxes/Hu_L')
clc
for sub = 1:72
    tempD1 = VS.ERP{sub};
    for ch = 1:size(tempD1,1)
        disp(['Processing sub: ',num2str(sub), ',  ch: ',num2str(ch)]);  
        Fs = 500;
        x1 = squeeze(tempD1(ch,:,:)); % [Time × Trials]
        xtimes = TimeStamps/1000; % -0.5:0.002:0.998
        t_stft = xtimes; 
        f_stft = 0 : 0.5 : 48; 
        winsize = 0.2; % 200 ms window
        [~,P1] = sub_stft(x1, xtimes, t_stft, f_stft, Fs, winsize);
        baseline_idx = (t_stft >= -0.4 & t_stft <= -0.05);
        baseline_power1 = mean(P1(:, baseline_idx, :), 2);
        P_dB1 = 10*log10(P1 ./ baseline_power1);
        P_dB_All_VS(sub,ch,:,:)=mean(P_dB1,3);
        clear trial P_dB1 P1 baseline_power1 x1 P_dB2 P2 baseline_power2 x2
    end
    clear tempD1 tempD2
end

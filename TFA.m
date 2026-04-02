addpath('/Users/Re-Re/Documents/MATLAB_ToolBoxes/Hu_L')
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
        winsize = 0.3; % 300 ms window
        [~,P1] = sub_stft(x1, xtimes, t_stft, f_stft, Fs, winsize);
        [~,P2] = sub_stft(x2, xtimes, t_stft, f_stft, Fs, winsize);       
        baseline_idx = (t_stft >= -0.45 & t_stft <= -0.05);
        
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
addpath('/Users/Re-Re/Documents/MATLAB_ToolBoxes/Hu_L')
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
        winsize = 0.3; % 300 ms window
        [~,P1] = sub_stft(x1, xtimes, t_stft, f_stft, Fs, winsize);
        
        baseline_idx = (t_stft >= -0.45 & t_stft <= -0.05);
        
        baseline_power1 = mean(P1(:, baseline_idx, :), 2);

        P_dB1 = 10*log10(P1 ./ baseline_power1);
 
        P_dB_All_VS(sub,ch,:,:)=mean(P_dB1,3);
        clear trial P_dB1 P1 baseline_power1 x1 P_dB2 P2 baseline_power2 x2
    end
    clear tempD1 tempD2
end

%%
figure
for i = 1:26
subplot(5,6,i)
imagesc(t_stft,f_stft,squeeze(mean(P_dB_All_cm(:,i,:,:),1)))
colormap jet
axis xy
caxis([-3 3])
ylim([0.5 45])
xlim([-0.3 0.805])
title(ChanLocs(i).labels)
end
%%
% FL: F3 (3), Fz (2), and F4 (24)   P3 (13), Pz (11), and P4 (17)
figure
subplot(2,1,1)
imagesc(t_stft,f_stft,squeeze(mean(P_dB_All_cm(:,[2 3 11 13 17 24],:,:),[1 2])))
colormap jet
axis xy
caxis([-4 3])
ylim([0.5 45])
xlim([-0.2 0.8])
colorbar
title('Cm')
subplot(2,1,2)
imagesc(t_stft,f_stft,squeeze(mean(P_dB_All_incm(:,[2 3 11 13 17 24],:,:),[1 2])))
colormap jet
axis xy
caxis([-4 3])
ylim([0.5 45])
xlim([-0.2 0.8])
colorbar
title('Incm')
% VS: P7 13, P3 12, Pz 11, P4 17, and P8 18; O1 14, Oz 15, and O2 16;
%F3 3, Fz 2, F4 24, FC1 6, and FC2 23; 
figure
imagesc(t_stft,f_stft,squeeze(mean(VS_dB_All(:,[13 12 11 17 18 14 15 16 3 2 24 6 23],:,:),[1 2])))
colormap jet
axis xy
caxis([-4 3])
ylim([0.5 45])
xlim([-0.2 0.8])
colorbar
%%
ROI_Value(:,2) = squeeze(mean(P_dB_All_cm(:,[2 3 11 13 17 24],f_stft>=8 & f_stft< 30,t_stft>=0.2 & t_stft<=0.5),[2 3 4]));
ROI_Value(:,1) = squeeze(mean(VS_dB_All(:,[13 12 11 17 18 14 15 16 3 2 24 6 23],f_stft>=8 & f_stft< 30,t_stft>=0.2 & t_stft<=0.6),[2 3 4]));
[~,p]=ttest(ROI_Value(:,1),ROI_Value(:,2))
%%
clc
P = mafdr(A,'BHFDR',true);
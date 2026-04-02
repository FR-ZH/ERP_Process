eeglab
close all
clear
clc
All_Sub_Names = arrayfun(@(x) sprintf('%03d', x), 1:127, 'UniformOutput', false);
Folder = '/Desktop/Researches/SES/ds006018-download';
for sub = 1:numel(All_Sub_Names)
    clc
    disp(All_Sub_Names{sub})
    subFolfer = [Folder,'/sub-',All_Sub_Names{sub},'/eeg'];
    vhdrFiles = dir(fullfile(subFolfer, '*.vhdr'));
    for c = 1:length(vhdrFiles)
        work_file_name =  vhdrFiles(c).name;
        EEG = pop_loadbv([subFolfer,'/'],work_file_name);
        EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',45);
        [P,Frq]=spectopo(EEG.data, 0, EEG.srate, 'freqrange', [1 45], 'plot', 'off');
        P_1=mean(P(:,Frq > 1 & Frq<10),2);
        P_2=mean(P(:,Frq>= 10 & Frq<20),2);
        P_3=mean(P(:,Frq>= 20 & Frq<45),2);
        Aver(1)=mean(P_1);
        Aver(2)=mean(P_2);
        Aver(2)=mean(P_3);
        STD(1)=std(P_1);
        STD(2)=std(P_2);
        STD(3)=std(P_3);
        Bad_CH =[];
        Bad_CH1=[];
        Bad_CH2=[];
        Bad_CH3=[];
        for ch = 1:EEG.nbchan
            if (P_1(ch) <= Aver(1) - 3*STD(1)) || (P_1(ch) >= Aver(1) + 3*STD(1))
                Bad_CH1=[Bad_CH1 ch];
            end
            if (P_2(ch) <= Aver(2) - 3*STD(2)) || (P_3(ch) >= Aver(2) + 3*STD(2))
                Bad_CH2=[Bad_CH2 ch];
            end
            if (P_3(ch) <= Aver(3) - 3*STD(3)) || (P_3(ch) >= Aver(2) + 3*STD(3))
                Bad_CH2=[Bad_CH2 ch];
            end
        end
        Bad_CH = union(Bad_CH1, Bad_CH2,Bad_CH3);
        clc
        if ~isempty(Bad_CH)
            EEG = pop_interp(EEG, Bad_CH, 'spherical');
        end
        EEG = pop_reref( EEG, [10 21] );
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'rndreset','yes');
        switch length(work_file_name)
            case 29
                task_name = 'FL';
            case 34
                task_name = 'VS';
        end
        EEG.setname=['Sub',All_Sub_Names{sub},'_',task_name,'_1to45Hz_reRef'];
        pop_saveset( EEG, 'filename',['Sub',All_Sub_Names{sub},'_',task_name,'_ICA.set'],...
            'filepath',['/Users/Re-Re/Desktop/Researches/SES/ICA/',task_name,'/']);
        clear  work_file_name task_name P Frq Bad_CH Bad_CH1 Bad_CH2 Aver STD P_1 P_2 P_3 EEG ch
        pause(0.1)
    end
    clear c vhdrFiles subFolfer
end
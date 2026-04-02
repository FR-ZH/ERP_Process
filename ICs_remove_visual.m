eeglab
close all
clear
clc
Folder = '/Users/Re-Re/Desktop/Researches/SES/ICA/MMN';
SaveFolder = '/Users/Re-Re/Desktop/Researches/SES/ICA/MMN/rmv_ICs';
setFiles = dir(fullfile(Folder, '*.set'));
screenSize = get(0, 'ScreenSize'); 
for i = 26:length(setFiles)
    clc
    workFile = setFiles(i).name;
    EEG = pop_loadset('filename',workFile,'filepath',[Folder,'/']);
    disp([num2str(i),'/',num2str(length(setFiles))])
    pop_topoplot(EEG, 0, [1:15] ,workFile(1:end-4),[4 5] ,0,'electrodes','on');
    hFig = gcf;
    set(hFig, 'Units', 'pixels', 'Position', [0 800 650 600]);
    Answer0 = inputdlg('Suspect ICs (no more than 8)', 'Suspect', [1 50], {''});
    Answer1 = str2num(Answer0{1});
    close all
    for n = 1:length(Answer1)
        pop_prop( EEG, 0, Answer1(n), NaN, {'freqrange',[2 50] });
        hFig = gcf;
        figWidth = floor((screenSize(3)-10)/4);
        figHeight = 430;
        if n<=4
            newPosition = [(n-1)*(figWidth+1), screenSize(4) - figHeight, figWidth, figHeight];
            set(hFig, 'Units', 'pixels', 'Position', newPosition);
        else
            newPosition = [(n-5)*(figWidth+1), screenSize(4)-550 - figHeight, figWidth, figHeight];
            set(hFig, 'Units', 'pixels', 'Position', newPosition);
        end
    end
    clc
    Answer2 = inputdlg('Input ICs to delete', 'Remove ICs', [1 50], {Answer0{1}});
    EEG = pop_subcomp( EEG, str2num(Answer2{1}), 0);
    close all
    EEG.setname=[workFile(1:end-4),'_Remove_ICs'];
    pop_saveset( EEG, 'filename',[workFile(1:end-4),'_rmvICs.set'],'filepath',[SaveFolder,'/']);
    clear EEG Answer0 Answer1 Answer2 n hFig newPosition figWidth figHeight ans workFile
end
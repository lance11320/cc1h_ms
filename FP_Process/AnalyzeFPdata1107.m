clear
close all

home = 'V:\FP_data\usedata\';
load('V:\FP_data\usedata\20231019-female 628 diestrus-data.mat')


%%

ori_trial = alldata{1};
crf_trial = alldata{2};
glu_trial = alldata{3};
oxy_trial = alldata{4};

offset_dura = [5000:8001];
onset_dura = [2000:5000];
peak1ori = max(ori_trial(:,onset_dura),[],2);
peak1crf = max(crf_trial(:,onset_dura),[],2);
peak1glu = max(glu_trial(:,onset_dura),[],2);
peak1oxy = max(oxy_trial(:,onset_dura),[],2);
checkpeak1ori = findmaxpeak(ori_trial(:,onset_dura));
checkpeak1crf = findmaxpeak(crf_trial(:,onset_dura));
checkpeak1glu = findmaxpeak(glu_trial(:,onset_dura));
checkpeak1oxy = findmaxpeak(oxy_trial(:,onset_dura));
peak2ori = max(ori_trial(:,offset_dura),[],2);
peak2crf = max(crf_trial(:,offset_dura),[],2);
peak2glu = max(glu_trial(:,offset_dura),[],2);
peak2oxy = max(oxy_trial(:,offset_dura),[],2);
checkpeak2ori = findmaxpeak(ori_trial(:,offset_dura));
checkpeak2crf = findmaxpeak(crf_trial(:,offset_dura));
checkpeak2glu = findmaxpeak(glu_trial(:,offset_dura));
checkpeak2oxy = findmaxpeak(oxy_trial(:,offset_dura));

for i = 1:length(peak2crf)
    if peak2ori(i) ~= checkpeak2ori(i)
        disp('Detect small pks in ori')
        checkpeak2ori = findminlocpeak(ori_trial(i,offset_dura));
    end
    if peak2crf(i) ~= checkpeak2crf(i)
        disp('Detect small pks in crf')
        checkpeak2crf(i) = findminlocpeak(crf_trial(i,offset_dura));
    end
    if peak2glu(i) ~= checkpeak2glu(i)
        disp('Detect small pks in glu')
        checkpeak2glu(i) = findminlocpeak(glu_trial(i,offset_dura));
    end
    if peak2oxy(i) ~= checkpeak2oxy(i)
        disp('Detect small pks in oxy')
        checkpeak2oxy(i) = findminlocpeak(oxy_trial(i,offset_dura));        
    end
end

for i = 1:length(peak1crf)
    if peak1ori(i) ~= checkpeak1ori(i)
        disp('Detect small pks in ori')
        checkpeak1ori = findminlocpeak(ori_trial(i,onset_dura));
    end
    if peak1crf(i) ~= checkpeak1crf(i)
        disp('Detect small pks in crf')
        checkpeak1crf(i) = findminlocpeak(crf_trial(i,onset_dura));
    end
    if peak1glu(i) ~= checkpeak1glu(i)
        disp('Detect small pks in glu')
        checkpeak1glu(i) = findminlocpeak(glu_trial(i,onset_dura));
    end
    if peak1oxy(i) ~= checkpeak1oxy(i)
        disp('Detect small pks in oxy')
        checkpeak1oxy(i) = findminlocpeak(oxy_trial(i,onset_dura));        
    end
end
%%
function [maxpeak] = findmaxpeak(datamat)
maxpeak = [];
for i = 1:size(datamat,1)
    pks = findpeaks(datamat(i,:));
    maxpeak(i) = max(pks);
end
end
function [minlocpk] = findminlocpeak(datamat)
minlocpk = [];
for i = 1:size(datamat,1)
    [pks,loc] = findpeaks(datamat(i,:));
    minloc = min(loc);
    minlocpk = pks(1);
end
end
%%
% compatible with computers with lower memory; 64GB memory could be
% recommended; depending on your video size.
clear

animal = {'86'};%
sess = {{'1'},{'1','2'},{'2'},{'2'},{'1','2'},{'1','2'},{'1','2'}};
state = {'ESTRUS'};

home = 'I:\CXM_Miniscope T1\ESTRUS\AHN';
writevid = 1;
for ii = 1:length(animal)
    for j = 1:length(state)
        session = sess{ii};
        for k = 1:length(session)
            dpath = [home,animal{ii},'\',state{j},'\Sess',session{k},'Res'];
            if exist(dpath,'dir')
                disp(['Now Processing ',dpath])
                % write new video (motion corrected) and use the video to
                % extract
                SaveVidPath = [dpath,'\minian_all.avi'];
                if exist(SaveVidPath,'file') ~= 0
                    writevid = 0;
                end
                if writevid
                    file1 = load([dpath,'\varr1']);
                    f1len = size(file1.array,1);
                    vid = VideoWriter(SaveVidPath,'Grayscale AVI');
                    open(vid)
                    for i = 1:f1len
                        writeVideo(vid,squeeze(file1.array(i,:,:)))
                    end
                    disp('file1 finished')
                    for idx = 2:100

                        if exist([dpath,'\varr',num2str(idx),'.mat'],'file')==2
                            disp(['file',num2str(idx),'detected'])
                            file3 = load([dpath,'\varr',num2str(idx)]);
                            f3len = size(file3.array,1);
                            for i = 1:f3len
                                writeVideo(vid,squeeze(file3.array(i,:,:)))
                            end
                        end
                    end
                    close(vid)
                    disp('done')
                end
                %%
                try
                    Apath = [dpath,'\',animal{ii},'A_add.mat'];
                    LocMat = load(Apath); LocMat = LocMat.array;
                    disp('Using Manually Selected Footprint to Extract')
                catch
                    Apath = [dpath,'\',animal{ii},'A.mat'];
                    LocMat = load(Apath); LocMat = LocMat.array;
                    disp('Using minian Selected Footprint to Extract')
                end

                %extract trace fast; only suit for work station with large
                %memory
                NeuTraceMat = func_GetTraceFast(LocMat,SaveVidPath);

                save([dpath,'\','NeuTrace.mat'],'NeuTraceMat')
                disp(['Saved to ',dpath,'\','NeuTrace.mat'])

            end
        end
    end
end
%send_email('Done in converting miniscope data into NeuTraceMat')
% Modified code for ReadAviAndndRegister
% Only suit for working station
function [NeuTraceMat] = func_GetTraceFast(LocMat,filename)

%% Load
tic
UnitLen = size(LocMat,1);
mov = VideoReader(filename);
FrameNum = mov.Duration*mov.FrameRate;
FrHeight = mov.Height;
FrWidth = mov.Width;
frame = read(mov);
Vdata = squeeze(frame(:,:,1,:));
clear frame
allV = Vdata;

toc
%% Calc
tic
disp('Finished loading! Calculating')
[k1,k2,k3] = size(allV);
[t1,t2,t3] = size(LocMat);
allvv = reshape(allV,k1*k2,k3);%turn H.W matrix to vector
lm = reshape(LocMat,t1,t2*t3);%turn H.W matrix to vector

seg_n = 5;
[starts,ends] = seg_seq(k3,seg_n);

for ix = 1:length(starts)
    partv = allvv(:,starts(ix):ends(ix));
    res = lm*im2double(partv);%0-1 matrix multiply
    var = sum(lm,2);
    var = repmat(var,1,size(res,2));
    if ix == 1
        NeuTraceMat = res./var;
    else
        NeuTraceMat = [NeuTraceMat res./var];
    end
end
%NeuTraceMat = im2uint8(NeuTraceMat);
NeuTraceMat = 255*im2double(NeuTraceMat);
toc
return
end

function [starts,ends]=seg_seq(total_len,seg_n)
if length(total_len)>1
    total_len = length(total_len);
end
starts = [1];
ends = ceil(total_len/seg_n);
for i = 1:seg_n-2
    starts = [starts i*ceil(total_len/seg_n)+1];
    ends = [ends (i+1)*ceil(total_len/seg_n)];
end
starts = [starts (seg_n-1)*ceil(total_len/seg_n)+1];
ends = [ends total_len];
end
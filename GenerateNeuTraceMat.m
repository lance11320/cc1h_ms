%%
clear


animal = {'408'};%
sess = {{'1'},{'1','2'},{'2'},{'2'},{'1','2'},{'1','2'},{'1','2'}};
state = {'Estrus','Diestrus'};

home = 'H:\NewFree\SortMS\M';
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
                if writevid
                    file1 = load([dpath,'\varr1']);
                    f1len = size(file1.array,1);
                    vid = VideoWriter(SaveVidPath,'Grayscale AVI');
                    open(vid)
                    for i = 1:f1len
                        writeVideo(vid,squeeze(file1.array(i,:,:)))
                    end
                    disp('file1 finished')
                    %                 file2 = load([dpath,'\varr2']);
                    %                 f2len = size(file2.array,1);
                    %                 for i = 1:f2len
                    %                     writeVideo(vid,squeeze(file2.array(i,:,:)))
                    %                 end
                    %                 disp('file2 finished')
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
                    Apath = [dpath,'\',animal{ii},'Chose_A.mat'];
                    LocMat = load(Apath); LocMat = LocMat.newarea;
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
Vdata = squeeze(frame(:,:,1,:));% squeeze to Height.Width.Frame
allV = Vdata;%(:,753:1504,:);
%NeuTraceMat = zeros(UnitLen,FrameNum);
toc
%% Calc
tic
disp('Finished loading! Calculating')
%LocMat(LocMat~=0)=1;
%LocMat(LocMat~=1)=0;
[k1,k2,k3] = size(allV);
[t1,t2,t3] = size(LocMat);
allvv = reshape(allV,k1*k2,k3);%turn H.W matrix to vector
lm = reshape(LocMat,t1,t2*t3);%turn H.W matrix to vector

partv1 = allvv(:,1:ceil(1/3*k3));
res = lm*im2double(partv1);%0-1 matrix multiply
var = sum(lm,2);
var = repmat(var,1,ceil(1/3*k3));
NeuTraceMat = res./var;

partv2 = allvv(:,(ceil(1/3*k3)+1):2*ceil(1/3*k3));
res = lm*im2double(partv2);%0-1 matrix multiply
var = sum(lm,2);
var = repmat(var,1,ceil(1/3*k3));
NeuTraceMat = [NeuTraceMat res./var];

partv3 = allvv(:,(2*ceil(1/3*k3)+1):end);
res = lm*im2double(partv3);%0-1 matrix multiply
var = sum(lm,2);
var = repmat(var,1,k3-2*ceil(1/3*k3));
NeuTraceMat = [NeuTraceMat res./var];

%NeuTraceMat = im2uint8(NeuTraceMat);
NeuTraceMat = 255*im2double(NeuTraceMat);
toc
return
end
% Manually extracting neurons

clear
concatidx = 0;

if concatidx
    animal = {'402'};%
    sess = {{'1','2'},{'1','2'},{'1','2'},{'1','2'},{'1','2'},{'1','2'},{'1','2'}};
    state = {'Estrus'};
    home = 'M:\MJH\NewFree\SortMS\M';
    for ii = 1:length(animal)
        for j = 1:length(state)
            session = sess{ii};
            for k = 1:length(session)
                dpath = [home,animal{ii},'\',state{j},'\Sess',session{k},'\My_V4_Miniscope\'];
                if exist(dpath,'dir')
                    SavePath = [dpath,'ConcatedAvi.avi'];
                    delete(SavePath)
                    %%
                    allf = dir([dpath,'*.avi']);

                    behloc = zeros(size(allf));
                    for i = 1:length(allf)
                        if contains(allf(i).name,'.avi')
                            behloc(i) = 1;
                        end
                    end

                    lenbehloc = length(find(behloc==1));
                    disp([num2str(lenbehloc),' Behavior Record is Found'])
                    %%
                    avi = VideoWriter(SavePath,'Motion JPEG AVI');
                    open(avi)
                    disp('Start Writing AVI Movie to the Same Direction')
                    for idx = 1:lenbehloc
                        disp(['Now Writing ',num2str(idx),'th Video'])
                        try
                            vidPath = [dpath,num2str(idx),'.avi'];
                            mov = VideoReader(vidPath);
                        catch
                            vidPath = [dpath,num2str(idx-1),'.avi'];
                            mov = VideoReader(vidPath);
                        end
                        while hasFrame(mov)
                            frame = readFrame(mov);
                            writeVideo(avi,frame)
                        end
                    end
                    close(avi)
                end
            end
        end
    end
end
%%

animalid = '391'; usemp = 1;
if usemp
    frame = load(['G:\NewData\SortMS\M',animalid,'\Diestrus\Sess1Res\',animalid,'max_proj.mat']);
    frame = frame.array;
else
    mov = VideoReader(['G:\NewData\SortMS\M',animalid,'\Diestrus\Sess2\My_V4_Miniscope\4.avi']);
    frame = readFrame(mov);
end
im = frame;
stopflag = 0;
Array = [];

originfootprint = load(['G:\NewData\SortMS\M',animalid,'\Diestrus\Sess1Res\',animalid,'A.mat']);
orfp = squeeze(sum(originfootprint.array,1)); orfp(orfp>0) = 1; orfp = logical(orfp);
im(orfp) = 0;

while stopflag==0
    imshow(im/100)
    roi = drawfreehand('color','r');
    mask = createMask(roi);
    %mask = roipoly(im);
    %pixel = regionprops(mask ,'all');
    im(mask) = 0;
    close gcf
    if size(mask,1)~=size(im,1)
        fpmat = zeros(size(im));
        fpmat(mask) = 1;
        mask = fpmat;
    end
    

    try
        if isempty(mask) || max(max(mask))==0
            error('No mask assigned')
        end
        Array(end+1,:,:) = mask;
        disp('Adding One Neuron')

    catch
        dostop = input('Please make sure that you have ended the extracting (y/n)','s');
        if strcmp(dostop,'y')
            stopflag = 1;
        else
            continue
        end
    end
end
array = cat(1,Array,originfootprint.array);
disp('Ending, Please save the Area matrix')
%%
% save('G:\NewData\SortMS\M391\Diestrus\Sess1Res\391A.mat','array')
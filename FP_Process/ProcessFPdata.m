close all
clear
home = 'V:\FP_data\';
date = '20231116';
animal = 'female 628 estrus';
%animal = 'male 709';
%state = {'OxtrA','GluRA', 'Before CRFR1A', 'CRFR1A', 'PTX', 'Before OxtrA'};
state = {'Before CRFR1A', 'CRFR1A', 'GluRA', 'OxtrA', 'mibefradil'};
%state = {'Before CRFR1A', 'CRFR1A', 'GluRA', 'OxtrA'};
%state = {'Before OxtrA', 'OxtrA'};
%state = {'GluRA','OxtrA'};

legendname = []; alldata = {};
addpath('colorplus\')
%% configs
%colors = {addcolorplus(1),addcolorplus(17),addcolorplus(115),addcolorplus(120)}; % each color for a state
colors = {'k','r','b','g','y'};
usePulse = false; % use the trials of pulse or not, if not, cw only
photoStimLength = 3; % the time duration of your photostimulation, in second.
preos = 2000; % preonset time
afterofs = 3000; % afteroffset time
TrialTrig = 'Photo'; % specify your trial trigger: could be 'Photo' or 'Behavior'
save2csv = false; % whether save  your data as the format compatible with old pipeline
dispinfo = false; % display your data groups in command line
%%
figure
title(animal)
for stidx = 1:length(state)

    filepath = [home,date,' test\',animal,'\',state{stidx},'\'];
    if exist(filepath,'dir')==0
        continue
    end
    allfiles = dir(filepath);
    alltrials = {};
    for i = 1:length(allfiles)
        if ~contains(allfiles(i).name,'DFF') && contains(allfiles(i).name,'.doric')
            fname = allfiles(i).name;
        else
            continue
        end
        DFFfilename = [fname(1:end-6),'_DFF.doric'];
        if exist([filepath,DFFfilename],'file')==0
            continue
        end
        if usePulse
            cond = ~(contains(fname,'hz') || contains(fname,'Hz') || contains(fname,'HZ'));
        else
            cond = contains(fname,'hz') || contains(fname,'Hz') || contains(fname,'HZ');
        end
        if cond
            continue
        end

        profp = FPdata();
        % config of the object
        profp.save2csv = save2csv; profp.dispinfo = dispinfo;
        profp.options.TrialTrig = TrialTrig;
        profp.options.preonset = preos; profp.options.afteroffset = afterofs;
        profp.filepath = filepath;


        profp = profp.makeobj(fname);
        frameRate = round(length(profp.dFFtime)/max(profp.dFFtime));
        disp(['The Approximate Framerate is ',num2str(frameRate),' Frames Per Second'])
        %profp.plotData();
        profp = profp.alignPhotoTime();
        %profp.plotintData();
        profp = profp.get_stimtrial();
        if isempty(profp.stimtrials)
            continue
        end
        for tlen = length(profp.stimtrials)
            alltrials{end+1} = profp.stimtrials{tlen};
        end

        %profp.plot_trial();
    end
    alllen = [];
    for i = 1:length(alltrials)
        itrial = alltrials{i}; triallen = size(itrial,2);
        alllen = [alllen triallen];
    end
    [sellen] = unique(alllen); 
    sellen = photoStimLength*1000+1+profp.options.preonset+profp.options.afteroffset;
    for i = 1:length(sellen)
        seltrials = [];
        for j = 1:length(alltrials)
            if size(alltrials{j},2)==sellen(i)
                seltrials = [seltrials; alltrials{j}];
            end
        end
        
        trials2plot = GetNorm(seltrials,1,profp.options.preonset); 
        %trials2plot = seltrials;
        obj = profp;
        if size(trials2plot,1)<2
            continue
        end
        
        alldata{end+1} = trials2plot;
        hi = plot(mean(trials2plot),'Color',colors{stidx},'LineWidth',1); hi.DisplayName = state{stidx};
        boundedline(1:size(trials2plot,2),mean(trials2plot),std(trials2plot)/sqrt(size(trials2plot,1)),'transparency',0.2,'alpha','cmap',colors{stidx});
        hold on
        
        line([obj.options.preonset obj.options.preonset],[-10 50],'Linestyle',':','Linewidth',2,'color','k')
        line([size(trials2plot,2)-obj.options.afteroffset size(trials2plot,2)-obj.options.afteroffset],[-10 50],'Linestyle',':','Linewidth',2,'color','k')
        xticks([0 obj.options.preonset size(trials2plot,2)-obj.options.afteroffset])
        xticklabels({'0s',['Light on(',num2str(obj.options.preonset/1000),'s)'],'Light off'})
        ylim([-3,28])
        set(gcf,'Color','w')
        hold on
        legendname = [legendname hi];
        legend(legendname)
    end

end
%save([home,date,'-',animal,'-data.mat',],'alldata')
%exportgraphics(gcf,'H:\OneDrive\SortScripts\Results\FP_Male_trace.emf','ContentType','vector')
function Normalized = GetNorm(NeuTraceMat,varargin)
if isempty(varargin)
    MEAN = mean(NeuTraceMat,2); Sigma = std(NeuTraceMat,0,2);
    Normalized = (NeuTraceMat - repmat(MEAN,1,size(NeuTraceMat,2)))./repmat(Sigma,1,size(NeuTraceMat,2));
else
    TimeS = varargin{1}; TimeE = varargin{2};
    BaselineTime = [TimeS:TimeE];
    Baseline = NeuTraceMat(:,BaselineTime); MEAN = mean(Baseline,2);Sigma = std(Baseline,0,2);
    Normalized = (NeuTraceMat - repmat(MEAN,1,size(NeuTraceMat,2)))./repmat(Sigma,1,size(NeuTraceMat,2));
end
end
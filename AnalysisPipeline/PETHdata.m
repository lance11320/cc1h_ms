classdef PETHdata
    properties
        dpath;PETHpath;
        NeuTraceMat; elab;
        cellAUC; cellID;
        trials; timeslides;
        params;
        Intermediate;
        saves;
        colorlist;
        beh;
        StartPoint;EndPoint;
        pseudoPETH;
        checkpoints;
        transparency
        ms_align2beh_time
        denoiseTrace; spikes;
        BehFile; BehtsFile; mstsFile;
        omitBeh;
        UseC;

    end

    methods
        function obj = PETHdata(obj,varargin)
            disp('对象已创建')
        end
        function obj = makeobj(obj,PETHpath,varargin)
            dosave = 0;
            if ~isempty(varargin)
                home = varargin{1}; root = varargin{2}; animal = varargin{3};
                state = varargin{4}; session = varargin{5};
                obj.dpath = [home,root,animal,'\',state,'\Sess',session];
                obj.saves.savepath = [obj.dpath, 'Res\'];
                obj.saves.savefig = [obj.dpath, 'Res\SmooUnitTrace\'];
                Arr = load([obj.dpath,'Res\','NeuTrace.mat']);
                obj.NeuTraceMat = Arr.NeuTraceMat;
                if obj.UseC
                    Arr = load([obj.dpath,'Res\',animal,'C.mat']);
                    obj.NeuTraceMat = Arr.array;
                end
                try
                    [aligned_elab,obj] = func_alignBeh(obj,[home,root],animal,state,session);
                    obj.elab = aligned_elab;
                    obj.elab = obj.elab(:,1:length(obj.NeuTraceMat));
                catch
                    if obj.omitBeh
                        obj.pseudoPETH = 'y';
                        disp('现在以无行为模式继续，AUC和trial不会被分析')
                    else
                        obj.pseudoPETH = input('你没有行为打标的结果或者打标文件为空，要使用无标记的结果继续吗？(输入y/n)','s');
                        disp('现在以无行为模式继续，AUC和trial不会被分析')
                    end
                    if strcmp(obj.pseudoPETH,'y')
                        obj.elab = zeros(1,length(obj.NeuTraceMat));
                    else
                        error('请进行打标或导出打标结果后继续')
                    end
                end

                PETHsavepath = [home,'MSRes\',root,animal,'\',state,'\Sess',session,'Res\'];
                mkdir(PETHsavepath)
                save([PETHsavepath,'PETH'],'obj');
                dosave = 1;
            end
            obj.PETHpath = PETHpath;
            currentdpath = obj.dpath;
            if dosave
                obj.PETHpath = [PETHsavepath,'PETH.mat'];
            end
            load(obj.PETHpath)
            obj.dpath = currentdpath;
            obj.saves.savepath = [obj.dpath, 'Res\'];
            obj.saves.savefig = [obj.dpath, 'Res\SmooUnitTrace\'];
        end

        function obj = matsmooth(obj,varargin)
            TraceMat = obj.NeuTraceMat;
            obj.Intermediate.NeuTracemat = obj.NeuTraceMat;
            if isempty(varargin)
                for i = 1:size(TraceMat,1)
                    Trace = smooth(TraceMat(i,:)',10);
                    TraceMat(i,:) = Trace';
                end
            else
                for i = 1:size(TraceMat,1)
                    Trace = smooth(TraceMat(i,:)',varargin{1});
                    TraceMat(i,:) = Trace';
                end
            end
            obj.NeuTraceMat = TraceMat;
            obj.checkpoints{end+1} = 'smoothed';
        end

        function [obj] = func_calcAUC(obj)
            behlist = unique(obj.elab);
            obj.cellAUC = cell(size(obj.elab,1),1);
            obj.cellID = cell(size(obj.elab,1),1);
            for idx = behlist(2):behlist(end)
                if length(obj.elab)>length(obj.NeuTraceMat)
                    obj.elab = obj.elab(:,1:size(obj.NeuTraceMat,2));
                end
                labs = obj.elab(idx,obj.elab(idx,:)==0 | obj.elab(idx,:)==idx);
                labs(obj.elab(idx,:)~=idx) = 0;
                labs(obj.elab(idx,:)==idx) = 1;
                if length(find(labs==1)) < 1
                    disp(['虽然',num2str(idx),'号行为被指定，但是并没有在过程中发生'])
                    continue
                end
                AUC = zeros(1,size(obj.NeuTraceMat,1));
                cellact = zeros(1,size(obj.NeuTraceMat,1));

                for id = 1:size(obj.NeuTraceMat,1)
                    tic
                    [~,~,~,auc] = perfcurve(labs,obj.NeuTraceMat(id,obj.elab(idx,:)==0 | obj.elab(idx,:)==idx),1);
                    AUC(1,id) = auc;
                    cons = zeros(1,1000);
                    % parfor JitterTime = 1:1000
                    %     I = randperm(length(obj.NeuTraceMat));
                    %     neumat = obj.NeuTraceMat(id,I);
                    %     [~,~,~,conauc] = perfcurve(labs,neumat,1);
                    %     cons(1,JitterTime) = conauc;
                    % end
                    % sortcons = sort(cons,'descend');
                    % uplim = sortcons(25);
                    % lowlim = sortcons(975);
                    % if auc>uplim %&& auc>0.7
                    %     cellact(1,id) = 1; 
                    % elseif auc<lowlim %&& auc<0.3
                    %     cellact(1,id) = -1;
                    % end
                    toc
                end
                obj.cellAUC{idx} = AUC;
                obj.cellID{idx} = cellact;
            end
            auc2save = obj.cellAUC; 
            cellid2save = obj.cellID;
            save([obj.saves.savepath,'AUC.mat'],'auc2save')
            save([obj.saves.savepath,'cellidentity.mat'],'cellid2save')
            obj.checkpoints{end+1} = 'AUC calc';
        end

        function [sorted,I] = func_sortresp(obj,varargin)
            behlist = unique(obj.elab);
            for idx = behlist(2):behlist(end)
                labs(obj.elab(idx,:)~=idx) = 0;
                labs(obj.elab(idx,:)==idx) = 1;
                if length(find(labs==1)) < 1
                    disp(['虽然',num2str(idx),'号行为被指定，但是并没有在过程中发生'])
                    continue
                end
                mtrial = func_getTrialData(obj,idx,obj.params.preonset,obj.params.afteronset);
                %% rearrange by neuron
                AverPSTH = [];
                for ineuron = 1:size(mtrial{1},1)
                    nresA = [];
                    for tid = 1:length(mtrial)
                        nres = mtrial{tid}(ineuron,:);
                        nresA = [nresA; nres];
                    end
                    if length(mtrial)==1
                        Aver_res = nresA;
                    else
                        Aver_res = mean(nresA);
                    end
                    AverPSTH = [AverPSTH;Aver_res];
                end
                %% sort

                if ~isempty(varargin)
                    pre = varargin{1};
                    post = varargin{2};
                else
                    pre = obj.params.preonset;
                    post = obj.params.afteronset;
                end
                AverPSTH = obj.Norms(AverPSTH,obj.params.baseline(1),obj.params.baseline(2));
                [~,I] = sort(mean(AverPSTH(:,pre:post),2));
                sorted = (AverPSTH(I,:));

                %% plot
                figure
                imagesc(sorted)
                colorbar
                clim([-5 5])
                line([pre pre],[0 size(AverPSTH,1)+1],'Linestyle',':','Linewidth',2,'color','k')
                yticks(1:1:size(AverPSTH,1)+1)
                yticklabels(I)
                title(['Sort Beh according to',obj.beh{idx}])
            end
        end

        function [mtrial,tslidesM,obj] = func_getTrialData(obj,behidx,varargin)
            if isempty(varargin)
                preonset = 480; afteronset = 360;
            else
                preonset = varargin{1};
                afteronset = varargin{2};
            end
            binsize = 1;
            [SniMStart,SniMEnd,obj] = func_getStartEnd(obj,behidx);
%             diff = SniMStart-SniMEnd; locs = find(diff<-30);
%             SniMStart = SniMStart(locs); SniMEnd = SniMEnd(locs);
            OnsetTraceA = {};
            tnum = 0;
            Mlen = length(SniMStart);
            try
                isempty(obj.params.trialgap);
            catch
                obj.params.trialgap = 90;
            end
            for tid = 1:Mlen
                if tid == 1 || SniMStart(tid)>SniMEnd(tid-1)+obj.params.trialgap
                    try
                        OnsetAct = obj.NeuTraceMat(:,SniMStart(tid)-preonset:SniMStart(tid)+afteronset);
                        OnsetTraceA{end+1} = OnsetAct;
                        tnum = tnum + 1;
                    catch
                        continue
                    end
                end
            end
            mtrial = OnsetTraceA;
            time_xtrainM = {};
            for idx = 1:binsize:preonset+afteronset
                xtrain = [];
                for ii = 1:length(OnsetTraceA)
                    xii = OnsetTraceA{ii}(:,idx:idx+binsize);
                    xtrain = [xtrain mean(xii,2)];
                end
                time_xtrainM{end+1} = xtrain;
            end
            tslidesM = time_xtrainM;

        end

        function obj = get_trial(obj)
            behlist = 0:max(unique(obj.elab));
            for idx = behlist(2):behlist(end)
                [mtrial,tslide,obj] = func_getTrialData(obj,idx,obj.params.preonset,obj.params.afteronset);
                obj.trials{idx} = mtrial;
                obj.timeslides{idx} = tslide;
            end
            obj.checkpoints{end+1} = 'Trial Got';
        end
        function Normalized = Norms(obj,neumat,varargin)
            if isempty(varargin)
                %% Baseline deltaf/f
                MEAN = mean(neumat,2); Sigma = std(neumat,0,2);
                Normalized = (neumat - repmat(MEAN,1,size(neumat,2)))./repmat(Sigma,1,size(neumat,2));
            else
                TimeS = varargin{1}; TimeE = varargin{2};
                BaselineTime = [TimeS:TimeE];
                Baseline = neumat(:,BaselineTime); MEAN = mean(Baseline,2);Sigma = std(Baseline,0,2);
                Normalized = (neumat - repmat(MEAN,1,size(neumat,2)))./repmat(Sigma,1,size(neumat,2));
            end
        end

        function obj = func_CalcDeltaf(obj,varargin)
            % zscore neutrace
            if isempty(varargin)
                %% Baseline deltaf/f
                MEAN = mean(obj.NeuTraceMat,2); Sigma = std(obj.NeuTraceMat,0,2);
                Normalized = (obj.NeuTraceMat - repmat(MEAN,1,size(obj.NeuTraceMat,2)))./repmat(Sigma,1,size(obj.NeuTraceMat,2));
            else
                TimeS = varargin{1}; TimeE = varargin{2};
                BaselineTime = [TimeS:TimeE];
                Baseline = obj.NeuTraceMat(:,BaselineTime); MEAN = mean(Baseline,2);Sigma = std(Baseline,0,2);
                Normalized = (obj.NeuTraceMat - repmat(MEAN,1,size(obj.NeuTraceMat,2)))./repmat(Sigma,1,size(obj.NeuTraceMat,2));
            end
            obj.Intermediate.NeuTracemat = obj.NeuTraceMat;
            obj.NeuTraceMat = Normalized;
            obj.checkpoints{end+1} = 'zscored';
        end

        function obj = func_percentDeltaf(obj)
            maxresp = max(obj.NeuTraceMat,[],2);
            Normalized = obj.NeuTraceMat./(repmat(maxresp,1,size(obj.NeuTraceMat,2)));
            obj.Intermediate.NeuTracemat = obj.NeuTraceMat;
            obj.NeuTraceMat = Normalized;

        end

        function  PlotAllTrace(obj,varargin)
            % Plot the neuron traces and return the epochs of sniffs;
            xlen = [1:size(obj.NeuTraceMat,2)]/30;
            for N = 1:size(obj.NeuTraceMat,1)
                ylim(1) = max(obj.NeuTraceMat(N,:));
                ylim(2) = min(obj.NeuTraceMat(N,:));

                scrsz = get(0,'ScreenSize');
                figure1 = figure('Position',[0 30 scrsz(3) scrsz(4)-105]);
                plot(xlen,obj.NeuTraceMat(N,:),'linewidth',1)
                behlist = unique(obj.elab);
                hold on
                all_handle = [];
                for ix = 2:length(behlist)
                    [SniMStart,SniMEnd] = func_getStartEnd(obj,behlist(ix));
                    SniMStart = SniMStart/30; SniMEnd = SniMEnd/30;
                    for j = 1:length(SniMStart)

                        hj = fill([SniMStart(j) SniMStart(j) SniMEnd(j) SniMEnd(j)],[ylim(2) ylim(1) ylim(1) ylim(2)],obj.colorlist{behlist(ix)});
                        set(hj,'edgealpha',0,'facealpha',obj.transparency)
                        if j==1
                            hj.DisplayName = obj.beh{ix-1};
                            all_handle = [all_handle hj];
                        end
                    end
                    hold on
                end
                if ~isempty(all_handle)
                    legend(all_handle)
                end
                xticks(1:120:length(obj.NeuTraceMat)/30)
                if obj.saves.saves
                    savepath = obj.saves.savefig;
                    if isempty(dir(savepath))
                        mkdir(savepath)
                    end
                    savefig(gcf,[savepath,'\UnitCalcium',num2str(N),'.fig'])
                    print(gcf,'-r600','-dpng',[savepath,'\UnitCalcium',num2str(N),'.png']);

                end
                close gcf

            end
        end

        function func_PlotCombineTrace(obj,varargin)
            Neu = obj.NeuTraceMat;
            scrsz = get(0,'ScreenSize');
            figure1 = figure('Position',[0 30 scrsz(3) scrsz(4)-105]);
            xlen = [1:size(Neu,2)]/30;
            ylabs = [];
            for ii = 1:size(Neu,1)
                if ii > 1
                    w = (ii-1)*5;
                    plot(xlen,Neu(ii,:)+w)
                    ylabs = [ylabs,w];
                else
                    plot(xlen,Neu(ii,:))
                    ylabs = [ylabs 0];
                end
                hold on
            end
            behlist = 0:length(obj.beh);
            ylim(1) = max(Neu(ii,:)+w);
            ylim(2) = min(Neu(1,:));
            all_handle = [];
            name2disp = {};

            for ix = 2:length(behlist)
                [SniMStart,SniMEnd] = func_getStartEnd(obj,behlist(ix));
                SniMStart = SniMStart/30; SniMEnd = SniMEnd/30;
                for j = 1:length(SniMStart)
                    hj = fill([SniMStart(j) SniMStart(j) SniMEnd(j) SniMEnd(j)],[ylim(2) ylim(1) ylim(1) ylim(2)],obj.colorlist{behlist(ix)});
                    set(hj,'edgealpha',0,'facealpha',obj.transparency)
                    if j==1
                        hj.DisplayName = obj.beh{ix-1};
                        name2disp{end+1} = obj.beh{ix-1};
                        all_handle = [all_handle hj];
                    end
                end
                hold on
            end
            if ~isempty(all_handle)
                legend(all_handle,name2disp)
            end
            if obj.params.align2behtime
                xticks(obj.ms_align2beh_time(obj.ms_align2beh_time~=0)/30)
            else
                xticks(1:120:length(Neu)/30)
            end
            yticks(ylabs);
            yticklabels(1:size(Neu,1));
            title(obj.dpath)
            

            if obj.saves.saves
                savepath = obj.saves.savefig;
                if isempty(dir(savepath))
                    mkdir(savepath)
                end
                savefig(gcf,[savepath,'\CombinedTrace.fig'])
                print(gcf,'-r600','-dpng',[savepath,'\CombinedTrace.png']);
            end
            figure
            ylim(1) = 3;
            ylim(2) = -0.5;
            imagesc(Neu)
            colorbar
            hold on
            all_handle = [];
            for ix = 2:length(behlist)
                [SniMStart,SniMEnd] = func_getStartEnd(obj,behlist(ix));
                for j = 1:length(SniMStart)
                    hj = fill([SniMStart(j) SniMStart(j) SniMEnd(j) SniMEnd(j)],-1*[ylim(2) ylim(1) ylim(1) ylim(2)],obj.colorlist{behlist(ix)});
                    set(hj,'edgealpha',0,'facealpha',obj.transparency)
                    if j==1
                        hj.DisplayName = obj.beh{ix-1};
                        all_handle = [all_handle hj];
                    end
                end
                
                hold on
            end
            if ~isempty(all_handle)
                legend(all_handle,obj.beh)
            end
            close gcf

        end

        function [Start,End,obj] = func_getStartEnd(obj,idx)
            % Get the timepoints of behavior onset and offset
            elabx = obj.elab; elabx = elabx(idx,:);
            elabx(elabx~=idx) = 0;
            elab1 = [0 elabx];
            EventP = [elabx 0] - elab1;
            Start = find(EventP==idx);
            End = find(EventP==-1*idx)-1;
            obj.StartPoint{end+1} = Start;
            obj.EndPoint{end+1} = End;
            obj.checkpoints{end+1} = 'StartEndPoint Got';

        end

        function [elab1,obj] = func_getBeh(obj,BehaviorFile)
            % load behavior to generate elab
            try
                opt = detectImportOptions(BehaviorFile);
                Beh_f = readtable(BehaviorFile,opt,'ReadVariableNames',true);
                len = length(Beh_f.RecordingTime);
            catch
                opt = detectImportOptions(BehaviorFile);
                opt.VariableNames = {'TrialTime' 'RecordingTime' 'Subject' 'Behavior' 'Event'};
                Beh_f = readtable(BehaviorFile,opt,'ReadVariableNames',true);
                for rows = 1:size(Beh_f,1)
                    if strcmp(Beh_f.Event(rows),'state start') || strcmp(Beh_f.Event(rows),'状态开始')
                        startrow = rows;
                        break
                    end
                end
                Beh_f = Beh_f(startrow:size(Beh_f,1),:);
                if strcmp(Beh_f.Event(1),'state stop')
                    Beh_f = readtable(BehaviorFile,opt,'ReadVariableNames',true);
                    Beh_f = Beh_f(42:size(Beh_f,1),:);
                end
            end
            if isempty(obj.beh)
                Beh_list = unique(Beh_f.Behavior);
                obj.beh = Beh_list;
            else
                Beh_list = obj.beh;
            end
            Behnum = length(Beh_list);
            for i = 1:Behnum
                disp(['The ',num2str(i),'th',' Behavior is ',Beh_list{i}])
            end
            elab1 = zeros(Behnum,18000);
            for idx = 1:Behnum
                sub_Behid = strcmp(Beh_f.Behavior,Beh_list{idx});
                sub_bf = Beh_f(sub_Behid,:);
                len = length(sub_bf.RecordingTime);
                Start = [1:2:len-1]';
                
                RecTime = str2double(sub_bf.RecordingTime);
                if isnan(RecTime)
                    RecTime = sub_bf.RecordingTime;
                end
                start_time = RecTime(Start);
                Stop = [2:2:len]';
                stop_time = RecTime(Stop);

                for i = 1:length(Start)
                    timark = ceil(30*start_time(i)):floor(30*stop_time(i)); timark = timark';
                    if isempty(timark)
                        continue
                    end
                    if timark(1)==0
                        timark = timark(2:end);
                    end
                    elab1(idx,timark) = idx;
                end
            end
            obj.elab = elab1;

        end

        function [aligned_elab,obj] = func_alignBeh(obj,dpath,animal,state,session)
            % align the behavior time stamp to miniscope time
            BehaviorFile = obj.BehFile; timestampFile = obj.BehtsFile; msstampFile = obj.mstsFile;
            if isempty(obj.BehFile)
                BehaviorFile = [dpath,animal,'\',state,'\Behavior',session,'\Export Files\','Beh.csv'];
            end
            if isempty(obj.BehtsFile)
                timestampFile = [obj.dpath,'\behcam\timeStamps.csv'];
            end
            if isempty(obj.mstsFile)
                msstampFile = [obj.dpath,'\timeStamps.csv'];
            end
            if  exist(BehaviorFile,'file')==2 && exist(timestampFile,'file')==2
                [~,obj] = func_getBeh(obj,BehaviorFile);
                elab_origin = obj.elab;
                behidx = unique(elab_origin);
                timestamp_beh = readtable(timestampFile);
                beh_tsmat = table2array(timestamp_beh);
                timestamp_ms = readtable(msstampFile);
                ms_tsmat = table2array(timestamp_ms);
                elabx = zeros(size(elab_origin,1),size(ms_tsmat(:,1),1));

                obj.ms_align2beh_time = zeros(size(1:450:length(ms_tsmat)));
                for points = 1:300:length(ms_tsmat)
                    [~,alignloc] = min(abs(beh_tsmat(:,2) - ms_tsmat(points,2)));
                    obj.ms_align2beh_time(points) = alignloc;
                end
                for id = 2:length(behidx)
                    [Start,End] = func_getStartEnd(obj,behidx(id));
                    [msloc,~] = find(beh_tsmat(:,1)==Start);
                    [meloc,~] = find(beh_tsmat(:,1)==End);

                    m_start = beh_tsmat(msloc,2);
                    m_end = beh_tsmat(meloc,2);

                    x_m_start = zeros(size(m_start)); x_m_end = x_m_start;
                    try
                        for iter = 1:length(m_start)
                            [x_m,~] = find(abs(ms_tsmat(:,2)-m_start(iter))<40); % a frame is around 30ms
                            x_m_start(iter) = x_m(1);
                            [x_m,~] = find(abs(ms_tsmat(:,2)-m_end(iter))<40);
                            x_m_end(iter) = x_m(1);
                        end
                    catch
                        disp('WARRNING! Align Failed!! Due to missing frame!!!')
                        for iter = 1:length(m_start)
                            [minval,x_m] = min(abs(ms_tsmat(:,2)-m_start(iter))); % a frame is around 30ms
                            x_m_start(iter) = x_m(1);
                            if minval>30
                                disp(['The frame gap is ',num2str(minval),'ms, Please check.'])
                            end
                            [minval,x_m] = min(abs(ms_tsmat(:,2)-m_end(iter)));
                            x_m_end(iter) = x_m(1);
                            if minval>30
                                disp(['The frame gap is ',num2str(minval),'ms, Please check.'])
                            end

                        end
                    end


                    for iter = 1:length(x_m_start)
                        elabx(behidx(id),x_m_start(iter):x_m_end(iter)) = behidx(id);
                    end
                end
            elseif exist([dpath,animal,'\',state,'\Sess',session,'\timestamp.dat'],'file')==2
                [elab_ori,obj] = func_getBeh(obj,BehaviorFile);
                behidx = unique(elab_ori);
                tsdat = importdata([dpath,animal,'\',state,'\Sess',session,'\timestamp.dat']);
                timestamp_all = tsdat.data(3:end,:);
                camera_series = sort(unique(timestamp_all(:,1)));
                beh_tsmat = timestamp_all(find(timestamp_all(:,1)==camera_series(1)),2:3);
                ms_tsmat = timestamp_all(find(timestamp_all(:,1)==camera_series(2)),2:3);
                elabx = zeros(size(elab_ori,1),size(ms_tsmat(:,1),1));
                for id = 2:length(behidx)
                    [Start,End,obj] = func_getStartEnd(obj,behidx(id));
                    [msloc,~] = find(beh_tsmat(:,1)==Start);
                    [meloc,~] = find(beh_tsmat(:,1)==End);

                    m_start = beh_tsmat(msloc,2);
                    m_end = beh_tsmat(meloc,2);

                    x_m_start = zeros(size(m_start)); x_m_end = x_m_start;

                    for iter = 1:length(m_start)
                        [x_m,~] = find(abs(ms_tsmat(:,2)-m_start(iter))<30); % a frame is around 30ms
                        x_m_start(iter) = x_m(1);
                        [x_m,~] = find(abs(ms_tsmat(:,2)-m_end(iter))<30);
                        x_m_end(iter) = x_m(1);
                    end


                    for iter = 1:length(x_m_start)
                        elabx(x_m_start(iter):x_m_end(iter)) = behidx(id);
                    end
                end

            end
            try
                aligned_elab = elabx;
            catch
                disp('Check your behavior path as the program cannot detect related file')
            end

        end

        function [svmres] = dosvmtrain(obj,svmparams)
            % this function has a realization of the method used in 'Social behaviour shapes hypothalamic
            % neural ensemble representations of conspecific sex' by David Anderson.
            Neu = obj.NeuTraceMat; repiter = 32; 
            if svmparams.pca > 0
                options.ReducedDim = svmparams.pca;
                eig = PCA(Neu,options); Neu = eig'*Neu;
            end
            svmres = zeros(repiter,size(obj.elab,1)); 
            parfor reps = 1:repiter
                behf1 = [];
                for idx = 1:size(obj.elab,1)
                    [starts,ends] = obj.func_getStartEnd(idx);
                    Accus = [];
                    randchosenum = 10;
                    if length(starts)<randchosenum
                        testtrial = 1:length(starts);
                        randchosenum = length(starts);
                    else
                        testtrial = randsample(1:length(starts),randchosenum);
                    end
                    for testiter = 1:randchosenum
                        testidx = testtrial(testiter);
                        trainidx = setdiff(1:length(starts),testidx);
                        trainset = [];
                        for ix = 1:length(trainidx)
                            trialid = trainidx(ix);
                            trainset = [trainset Neu(:,starts(trialid):ends(trialid))];
                        end
                        if isempty(trainset)
                            continue
                        end
                        numset = size(trainset,2); nullidx = randsample(find(obj.elab(idx,:)==0),numset);
                        trainset = [trainset Neu(:,nullidx)];
                        labels = [ones(1,numset)*idx,zeros(1,numset)];
                        model = svmtrain(labels',trainset','-t 0 -c 2 -q');

                        testset = Neu(:,starts(testidx):ends(testidx));
                        numtestset = size(testset,2); nulltest = randsample(find(obj.elab(idx,:)==0),numtestset);
                        testset = [testset Neu(:,nulltest)];
                        labeltest = [ones(1,numtestset)*idx,zeros(1,numtestset)];
                        pred = svmpredict(labeltest',testset',model);
                        if strcmp(svmparams.method,'F1')
                            accu = obj.f1_score(labeltest',pred);
                        elseif strcmp(svmparams.method,'auROC')
                            [~,~,~,accu] = perfcurve(labeltest',pred,idx);
                        else
                            accu = mean(pred==labeltest');
                        end
                        Accus = [Accus accu];
                    end
                    behf1 = [behf1 mean(Accus)];
                end
                svmres(reps,:) = behf1;
            end
            svmres = mean(svmres);

        end

        function [score, TPR, TNR] = f1_score(obj,label, predict)
            % F1 score calculation
            M = confusionmat(label, predict);
            TPR = M(2,2) / (M(2,1) + M(2,2)); 
            TNR = M(1,1) / (M(1,1) + M(1,2)); 
            M = M';
            precision = diag(M)./(sum(M,2) + 0.0001);
            recall = diag(M)./(sum(M,1)+0.0001)';
            precision = mean(precision);
            recall = mean(recall);
            score = 2*precision*recall/(precision + recall);
        end

        function [w,b] = get_w_b(obj,model,x,y)
            % calculate the weights and bias for SVM model
            SVs_idx = model.sv_indices;
            x_SVs = x(SVs_idx,:);
            y_SVs = y(SVs_idx);
            alpha_SVs = model.sv_coef;
            w = sum(diag(alpha_SVs)*x_SVs)';
            SVs_on = (abs(alpha_SVs)<1);
            y_SVs_on = y_SVs(SVs_on,:);
            x_SVs_on = x_SVs(SVs_on,:);
            b_temp = zeros(1,sum(SVs_on));
            for idx=1:sum(SVs_on)
                b_temp(idx) = 1/y_SVs_on(idx)-x_SVs_on(idx,:)*w;
            end
            b = mean(b_temp);
        end

        function [obj] = selfupdate(obj)

        end

        function [obj] = Gethelp(obj)
            disp('类有以下的方法：')
            methods(obj)
            hint = input('请输入你想获得帮助的函数','s');
            eval(['help ', hint])
        end

        function [obj] = getDeconvolve(obj)
            % use oasis package to get putative spikes
            oasis_setup
            spikex = zeros(size(obj.NeuTraceMat));
            denoisex = zeros(size(obj.NeuTraceMat));
            for i = 1:size(obj.NeuTraceMat,1)
                [denoised, spiketrain] = deconvolveCa(obj.NeuTraceMat(i,:),'ar1','foopsi', 'lambda', 1, 'optimize_pars');
                spikex(i,:) = spiketrain';
                denoisex(i,:) = denoised';
            end
            obj.denoiseTrace = denoisex;
            obj.spikes = spikex;
        end
        
    end

end


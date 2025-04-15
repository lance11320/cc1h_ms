classdef SVMdata < PETHdata
    properties
        w
        b
        x
        y
        model
        modelshu
        doshu
        svmplotflag
        svmres
        projected
        beh_projected
        samples
        svmparam
        donorm = 0;

    end

    methods
        function obj = copyProperties(obj,parentObj)
            props = properties(parentObj);
            for i = 1:length(props)
                prop = props{i};
                obj.(prop) = parentObj.(prop);
            end
        end
        function obj = SVMdata(inobj)
            obj = copyProperties(obj,inobj);
        end
        function redfea = PCAvisual(obj,fea,label,reducedim)
            options.ReducedDim = reducedim;
            [eigvec] = PCA(fea,options);
            redfea = fea*eigvec;
            scatter3(redfea(label==1,1),redfea(label==1,2),redfea(label==1,3))
            hold on
            scatter3(redfea(label==0,1),redfea(label==0,2),redfea(label==0,3))
        end

        function [obj] = trainSVM(obj,kind1,kind2)
            ploted = 0; 

            if obj.doshu
                maxiter = 50;
            else
                maxiter = 100;
            end

            auROC = [];
            OriNeu = obj.NeuTraceMat;
            svmobj = obj;
            Neu = OriNeu;
            % options.ReducedDim = 20;
            % eigvec = PCA(Neu',options);
            % redmat = eigvec'*Neu;
            if obj.samples == 0
                downsp = size(Neu,1);
            else
                downsp = obj.samples;
            end
            for iter = 1:maxiter
                if downsp == size(Neu,1)
                    randid = 1:size(Neu,1);
                else
                    randid = randsample(1:size(Neu,1),downsp);
                end
                redmat = Neu(randid,:);
                svmobj.NeuTraceMat = redmat;
                svmobj = svmobj.get_trial();
                femaletrial = svmobj.trials{kind1};
                maletrial = svmobj.trials{kind2};
                [feseqstart,feseqend] = obj.segseq(length(femaletrial),5);
                [maseqstart,maseqend] = obj.segseq(length(maletrial),5);
                for ix = 1:5
                    fetestidx = feseqstart(ix):feseqend(ix);
                    matestidx = maseqstart(ix):maseqend(ix);
                    fetrainidx = setdiff(1:length(femaletrial),fetestidx);
                    matrainidx = setdiff(1:length(maletrial),matestidx);
                    sniffmale = []; snifffemale = [];
                    for iters = matrainidx
                        sniffmale = [sniffmale obj.binseq(maletrial{iters}(:,120:181))];
                    end
                    for iters = fetrainidx
                        snifffemale = [snifffemale obj.binseq(femaletrial{iters}(:,120:181))];
                    end

                    dataset = [sniffmale snifffemale];
                    label = [zeros(size(sniffmale,2),1);ones(size(snifffemale,2),1)];

                    if obj.doshu
                        randidx = randperm(length(label));
                        obj.modelshu = svmtrain(label(randidx),dataset','-t 0 -c 2 -q');
                    else
                        obj.model = svmtrain(label,dataset','-t 0 -c 2 -q');
                        obj.x = dataset';
                        obj.y = label;
                    end

                    maletest = []; fetest = [];

                    for iters = matestidx
                        maletest = [maletest obj.binseq(maletrial{iters}(:,120:181))];
                    end
                    for iters = fetestidx
                        fetest = [fetest obj.binseq(femaletrial{iters}(:,120:181))];
                    end
                    testset = [maletest,fetest];
                    truelabel = [zeros(size(maletest,2),1);ones(size(fetest,2),1)];

                    if obj.doshu
                        pred = svmpredict(truelabel,testset',obj.modelshu);
                    else
                        pred = svmpredict(truelabel,testset',obj.model);
                    end
                    [~,~,~,clsAUC] = perfcurve(truelabel,pred,1);
                    auROC = [auROC clsAUC];
                end
                if clsAUC>0.7 && ploted==0 && obj.svmplotflag
                    ploted = 1;
                    figure
                    obj.PCAvisual(dataset',label,size(dataset,1));
                    %tSNEvisual(dataset,label)
                end
            end
            auROC = mean(auROC);
            obj.svmres = auROC;
        end

        function redmat = easyPCA(fea,dim)
            options.ReducedDim = dim;
            eigvec = PCA(fea',options);
            redmat = eigvec'*fea;
        end

        function obj = calc_hyperplane(obj)
            % calculate the weights and bias for SVM model
            svmmodel = obj.model;
            SVs_idx = svmmodel.sv_indices;
            x_SVs = obj.x(SVs_idx,:);
            y_SVs = obj.y(SVs_idx);
            alpha_SVs = svmmodel.sv_coef;
            obj.w = sum(diag(alpha_SVs)*x_SVs)';
            SVs_on = (abs(alpha_SVs)<1);
            y_SVs_on = y_SVs(SVs_on,:);
            x_SVs_on = x_SVs(SVs_on,:);
            b_temp = zeros(1,sum(SVs_on));
            for idx=1:sum(SVs_on)
                b_temp(idx) = 1/y_SVs_on(idx)-x_SVs_on(idx,:)*obj.w;
            end
            obj.b = mean(b_temp);
        end

        function obj = project_svm_direct(obj)
            obj.projected = obj.w' * obj.NeuTraceMat;
            resp = {};
            for ix = 1:size(obj.elab,1)
                resp{ix} = obj.projected(obj.elab(ix,:)==ix);
            end
            obj.beh_projected = resp;
        end

        function [seqstart,seqend] = segseq(obj,wholelength,numseq)
            seqstart = zeros(numseq,1); seqend = zeros(numseq,1);
            seqstart(1) = 1;
            seqend(numseq) = wholelength;
            for i = 1:numseq-1
                seqend(i) = i*floor(1/numseq*wholelength);
                seqstart(i+1) = i*floor(1/numseq*wholelength)+1;
            end

        end

        function [newseq] = binseq(obj,seq,varargin)
            if isempty(varargin)
                numbins = 12;
            else
                numbins = varargin{1};
            end
            seqstart = zeros(numbins,1); seqend = zeros(numbins,1);
            seqstart(1) = 1;
            seqend(numbins) = size(seq,2);
            newseq = zeros(size(seq,1),numbins);
            for i = 1:numbins-1
                seqend(i) = i*floor(1/numbins*size(seq,2));
                seqstart(i+1) = i*floor(1/numbins*size(seq,2))+1;
            end
            for i = 1:numbins
                newseq(:,i) = mean(seq(:,seqstart(i):seqend(i)),2);
            end 
        end

        function [tw_accu,tw_con_accu] = timewise_svm(obj,kind1,kind2)
            svmobj = obj;
            maxiter = obj.svmparam.maxiter; stepsize = obj.svmparam.stepsize;
            binsize = obj.svmparam.binsize;
            
            OriNeu = obj.NeuTraceMat;
            svmobj.params.preonset = 150; svmobj.params.afteronset = 150;
            Neu = OriNeu;

            if obj.samples == 0
                downsp = size(Neu,1);
            else
                downsp = obj.samples;
            end

            tw_accu = []; tw_con_accu = [];
            for iter = 1:maxiter
                if obj.svmparam.doPCA == 1
                    options.ReducedDim = obj.svmparam.pcadim;
                    eigvec = PCA(Neu',options);
                    redmat = eigvec'*Neu;
                else
                    if downsp == size(Neu,1)
                        randid = 1:size(Neu,1);
                    else
                        randid = randsample(1:size(Neu,1),downsp);
                    end
                    redmat = Neu(randid,:);
                end
                svmobj.NeuTraceMat = redmat;
                svmobj = svmobj.get_trial();
                femaletrial = svmobj.trials{kind1};
                maletrial = svmobj.trials{kind2};

                if svmobj.donorm == 1
                    for iters = 1:length(femaletrial)
                        femaletrial{iters} = svmobj.GetNorm(femaletrial{iters},1,svmobj.params.preonset);
                    end
                    for iters = 1:length(maletrial)
                        maletrial{iters} = svmobj.GetNorm(maletrial{iters},1,svmobj.params.preonset);
                    end
                end
                Maleslice = {}; Femaleslice = {};
                for idx = 1:stepsize:svmobj.params.preonset+svmobj.params.afteronset-svmobj.svmparam.binsize
                    xtrain = [];
                    for ii = randsample(1:length(femaletrial),svmobj.params.randtrialnum)
                        xii = femaletrial{ii}(:,idx:idx+binsize-1);
                        xtrain = [xtrain mean(xii,2)];
                    end
                    Maleslice{end+1} = xtrain;
                end
                for idx = 1:stepsize:svmobj.params.preonset+svmobj.params.afteronset-svmobj.svmparam.binsize
                    xtrain = [];
                    for ii = randsample(1:length(maletrial),svmobj.params.randtrialnum)
                        xii = maletrial{ii}(:,idx:idx+binsize-1);
                        xtrain = [xtrain mean(xii,2)];
                    end
                    Femaleslice{end+1} = xtrain;
                end

                accu_time = []; shu_accu_time = [];
                for time_i = 1:length(Maleslice)
                    mtrain = Maleslice{time_i};
                    ftrain = Femaleslice{time_i};
                    [segst,segend] = obj.segseq(size(mtrain,2),5);
                    ac_time = []; shu_ac_time = [];
                    for cviter = 1:5
                        testid = segst(cviter):segend(cviter);
                        trainid = setdiff(1:size(mtrain,2),testid);
                    % stid = 1:3:6*(binsize+1); edid = 4:3:6*(binsize+1);
                    % ac_time = []; shu_ac_time = [];
                    % for cviter = 1:5
                    %     testid = stid(cviter):edid(cviter);
                    %     trainid = setdiff(1:2*binsize,testid);
                        x_train = [ftrain(:,trainid) mtrain(:,trainid)]; y_train = [zeros(1,length(trainid)),ones(1,length(trainid))];
                        y_shu_id = randperm(length(y_train)); y_rd = y_train(y_shu_id);
                        x_test = [ftrain(:,testid) mtrain(:,testid)]; y_test = [zeros(1,length(testid)),ones(1,length(testid))];
                        model = svmtrain(y_train',x_train','-t 0 -c 2 -q');
                        model_shu = svmtrain(y_rd',x_train','-t 0 -c 2 -q');

                        pred = svmpredict(y_test',x_test',model);
                        accu = mean(y_test' == pred);
                        %[~,~,~,accu] = perfcurve(y_test',pred,1);

                        ac_time = [ac_time accu];

                        pred_shu = svmpredict(y_test',x_test',model_shu);
                        con_accu = mean(y_test' == pred_shu);
                        shu_ac_time = [shu_ac_time con_accu];
                    end
                    accu_time = [accu_time mean(ac_time)];
                    shu_accu_time = [shu_accu_time mean(shu_ac_time)];
                end
                tw_accu = [tw_accu ;accu_time];
                tw_con_accu = [tw_con_accu ; shu_accu_time];
            end
            tw_accu = mean(tw_accu); tw_con_accu = mean(tw_con_accu);

            accu_pre = tw_accu(:,1:1/2*size(tw_accu,2));
            accu_post = tw_accu(:,1/2*size(tw_accu,2)+1:end);
            % [lox,loy] = find(accu_pre>0.56);
            % for ii = 1:length(lox)
            %     tw_accu(lox(ii),loy(ii)) = tw_accu(lox(ii),loy(ii))-0.06;
            % end
        end

        function [Normalized] = GetNorm(obj,NeuTraceMat,varargin)
            if isempty(NeuTraceMat)
                Normalized = [];
                return
            end

            % zscore neutrace
            if isempty(varargin)
                %% Baseline deltaf/f
                MEAN = mean(NeuTraceMat,2); Sigma = std(NeuTraceMat,0,2);
                Normalized = (NeuTraceMat - repmat(MEAN,1,size(NeuTraceMat,2)))./repmat(Sigma,1,size(NeuTraceMat,2));
            else
                TimeS = varargin{1}; TimeE = varargin{2};
                BaselineTime = [TimeS:TimeE];
                Baseline = NeuTraceMat(:,BaselineTime); MEAN = mean(Baseline,2);Sigma = std(Baseline,0,2);
                Normalized = (NeuTraceMat - repmat(MEAN,1,size(NeuTraceMat,2)))./repmat(Sigma,1,size(NeuTraceMat,2));
            end
        end

    end

end
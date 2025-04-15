classdef FPdata
    properties
        filepath; fname
        dFFdata
        dFFtime
        oriref
        oridata
        photostim
        orireftime
        oridatatime
        photostimtime
        save2csv
        dispinfo
        intDFF; intStim; intDFFtime; intStimtime
        options
        stimtrials
    end
    methods
        function obj = FPdata()
            disp('---Set Object---')
        end
        function obj = makeobj(obj,filename)
            obj = obj.loadOri(filename);
            DFFfilename = [filename(1:end-6),'_DFF.doric'];
            obj = obj.loadDFF(DFFfilename);
        end
        function obj = loadDFF(obj,filename)

            obj.fname = filename;
            filename = [obj.filepath, filename];
            if obj.dispinfo
                disp('--- Please Check The Data Struct Below !!! ---')
                h5disp(filename,'/DataProcessed','min')
                disp('--- Please Focus on the Datasets, Ending. ---')
            end
            obj.dFFdata = h5read(filename,'/DataProcessed/FPConsole/DFFSignals/Series0001/AIN01xAOUT02-LockIn/Values');
            obj.dFFtime = h5read(filename,'/DataProcessed/FPConsole/DFFSignals/Series0001/AIN01xAOUT02-LockIn/Time');
            %     plot(dFFtime,dFFdata)
            if obj.save2csv
                csv2write = [obj.dFFtime, obj.dFFdata];
                writematrix(csv2write,[obj.filepath,'#',filename(1:end-6),' ','.csv'])
            end
            %     table2write = table(dFFtime,dFFdata);
            %     writetable(table2write,[filepath,'#',fname(1:end-6),' ','.csv'])
        end

        function obj = loadOri(obj,filename)

            filename = [obj.filepath, filename];
            if obj.dispinfo
                disp('--- Please Check The Data Struct Below !!! ---')
                h5disp(filename,'/DataAcquisition','min')
                disp('--- Please Focus on the Datasets, Ending. ---')
            end
            obj.oriref = h5read(filename,'/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Values');
            obj.oridata = h5read(filename,'/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT02-LockIn/Values');
            obj.orireftime = h5read(filename,'/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Time');
            obj.oridatatime = h5read(filename,'/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT02-LockIn/Time');
            obj.photostim = h5read(filename,'/DataAcquisition/FPConsole/Signals/Series0001/AnalogOut/AOUT04');
            obj.photostimtime = h5read(filename,'/DataAcquisition/FPConsole/Signals/Series0001/AnalogOut/Time');
            if obj.save2csv
                csv2write = [obj.oritime, obj.oridata];
                writematrix(csv2write,[obj.filepath,'#',obj.fname(1:end-6),' ','.csv'])
            end
            %     table2write = table(dFFtime,dFFdata);
            %     writetable(table2write,[filepath,'#',fname(1:end-6),' ','.csv'])
        end

        function obj = alignPhotoTime(obj)
            maxTimepoint = round(1e3*obj.dFFtime(end));
            maxStimTimePoint = round(1e3*obj.photostimtime(end));
            obj.intDFF = interp1(round(1e3*obj.dFFtime),obj.dFFdata,1:maxTimepoint);
            obj.intDFFtime = 1:maxTimepoint;
            obj.intStim = interp1(round(1e3*obj.photostimtime),obj.photostim,1:maxStimTimePoint); obj.intStim(obj.intStim<1)=0;
            obj.intStimtime = 1:maxStimTimePoint;

        end
        function obj = plotData(obj)
            tbldff = table(obj.dFFtime, obj.dFFdata);
            tbldff = renamevars(tbldff,["Var1","Var2"],["Time","DFFdata"]);
            tblphoto = table(obj.photostimtime, obj.photostim);
            tblphoto = renamevars(tblphoto,["Var1","Var2"],["Time","Stimdata"]);
            figure
            stackedplot(tblphoto,tbldff);
        end
        function obj = plotintData(obj)
            figure
            plot(obj.intDFFtime, obj.intDFF);
            hold on
            plot(obj.intStimtime, obj.intStim);
        end
        function obj = get_stimtrial(obj)
            elabx = obj.intStim;
            elab1 = [0 elabx];
            EventP = [elabx 0] - elab1;
            Starts = find(EventP==1);
            Ends = find(EventP==-1)-1;
            locs = [];
            if length(Starts) >= 2
                for i = 2:length(Starts)
                    if Starts(i)<Ends(i-1)+2000
                        continue
                    end
                    locs = [locs i];
                end
                Starts = Starts([1 locs]); Ends = [Ends(locs-1) Ends(end)];
            end

            obj.stimtrials = {};
            dura = [];
            for i = 1:length(Starts)
                dura = [dura Ends(i)-Starts(i)];
            end
            try
                if dura(1)~=mean(dura) && max(dura)-min(dura)>240
                    timetrial = ceil(max(dura)/1000);
                    if timetrial > 0
                        lentrial = timetrial*1000;
                        disp(['Trial length not matched, will be using fixed trial length, which is: ', num2str(timetrial),'s'])
                    else
                        lentrial = max(dura);
                        disp(['Trial length not matched, will be using fixed trial length, which is: ', num2str(lentrial),'ms'])
                    end
                elseif dura(1)~=mean(dura) && max(dura)-min(dura)<240
                    timetrial = ceil(max(dura)/1000);
                    if timetrial > 0
                        lentrial = timetrial*1000;
                        disp(['Trial length is: ', num2str(timetrial),'s'])
                    else
                        lentrial = max(dura);
                        disp(['Trial length is: ', num2str(lentrial),'ms'])
                    end
                elseif dura(1) == mean(dura) && max(dura)-min(dura) == 0
                    timetrial = ceil(max(dura)/1000);
                    if timetrial > 0
                        lentrial = timetrial*1000;
                        disp(['Perfectly matched! Trial length is: ', num2str(timetrial),'s'])
                    else
                        lentrial = max(dura);
                        disp(['Perfectly matched! Trial length is: ', num2str(lentrial),'ms'])
                    end
                end
                for i = 1:length(Starts)
                    obj.stimtrials{end+1} = obj.intDFF(Starts(i)-obj.options.preonset:Starts(i)+lentrial+obj.options.afteroffset);
                end
            catch
                disp('--- Warning!!! Please Check the Data, No Photostim Detected!!! ---')
            end
        end
        function obj = plot_trial(obj)
            trials2plot = [];

            for i = 1:length(obj.stimtrials)
                trials2plot = [trials2plot; obj.stimtrials{i}];
            end

            figure
            if size(trials2plot,1)>=2
                boundedline(1:size(trials2plot,2),mean(trials2plot),std(trials2plot)/sqrt(size(trials2plot,1)))
            else
                plot(trials2plot,'LineWidth',1)
            end
            line([obj.options.preonset obj.options.preonset],[-10 10],'Linestyle',':','Linewidth',2,'color','k')
            line([size(trials2plot,2)-obj.options.afteroffset size(trials2plot,2)-obj.options.afteroffset],[-10 10],'Linestyle',':','Linewidth',2,'color','k')
            xticks([0 obj.options.preonset size(trials2plot,2)-obj.options.afteroffset])
            xticklabels({'0s',['Light on(',num2str(obj.options.preonset/1000),'s)'],'Light off'})
            ylim([-3,10])
            set(gcf,'Color','w')
            mkdir([obj.filepath,'save\'])
            print(gcf,'-r600','-dpng',[obj.filepath,'save\',obj.fname,'Trial.png'])

        end

    end
end
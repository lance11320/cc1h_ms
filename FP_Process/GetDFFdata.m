% processing fiber photometry

clear 
close all

filepath = 'I:\20240130 FP\0202 VSDS\';
allfiles = dir(filepath);
for i = 1:length(allfiles)
    if contains(allfiles(i).name,'DFF') && ~contains(allfiles(i).name,'.csv')
        fname = allfiles(i).name;
    else
        continue
    end
    filename = [filepath, fname];
    disp('--- Please Check The Data Struct Below !!! ---')
    h5disp(filename,'/DataProcessed','min')
    disp('--- Please Focus on the Datasets, Ending. ---')
    dFFdata = h5read(filename,'/DataProcessed/FPConsole/DFFSignals/Series0001/AIN01xAOUT02-LockIn/Values');
    dFFtime = h5read(filename,'/DataProcessed/FPConsole/DFFSignals/Series0001/AIN01xAOUT02-LockIn/Time');
%     plot(dFFtime,dFFdata)

    csv2write = [dFFtime, dFFdata];
    writematrix(csv2write,[filepath,'#',fname(1:end-6),' ','.csv'])
%     table2write = table(dFFtime,dFFdata);
%     writetable(table2write,[filepath,'#',fname(1:end-6),' ','.csv'])
end
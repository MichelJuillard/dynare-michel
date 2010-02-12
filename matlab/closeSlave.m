function closeSlave(Parallel),
% In a parallelc context, this utility closes all remote matlab instances
% called by masteParallelMan (which leaves open remote matlab instances)

delete( 'slaveParallel_input*.mat');
for indPC=1:length(Parallel),
    if Parallel(indPC).Local==0,
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' rm -fr ',Parallel(indPC).RemoteFolder,'/slaveParallel_input*.mat']);
        else
            mydelete('slaveParallel_input*.mat',['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
        end
    end
end

function GAM = iteration(model,GAM,exp_data)

    fitting = ones(size(GAM))*1000;
    
    for i = 1:length(GAM)
        %Modify GAM:
        model_i = changeGAM(model,GAM(i),2.81);%NGAM=2.81 (Glucose)
        disp(['changed GAM:' num2str(GAM(i))])
        %Simulate model and calculate fitting:
        mod_data   = simulateChemostat(model_i,exp_data);
        R          = (mod_data - exp_data)./exp_data;
        fitting(i) = sqrt(sum(sum(R.^2)));
        disp(['GAM = ' num2str(GAM(i)) ' -> Error = ' num2str(fitting(i))])
    end
    
    %Choose best:
    [~,best] = min(fitting);
    
    if best == 1 || best == length(GAM)
        error('GAM found is sub-optimal: please expand GAM search bounds.')
    else
        GAM = GAM(best);
    end
end
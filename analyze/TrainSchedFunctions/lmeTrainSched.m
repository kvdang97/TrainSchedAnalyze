function lme = lmeTrainSched(data, var, func)
% data - data in a table format
% var - enter the dependent variable to plot as a string e.g. 'RT'
% func = 0 for mean data, 1 for variance data
if func == 0
    %%% AVERAGE DATA for every subject and every abs(tgt rotation)
    subj_mean = varfun(@nanmean, data, 'GroupingVariables',{'SN','Group','abs_tgt_rot'},'OutputFormat','table'); 
    
    subj_mean.Group = categorical(subj_mean.Group);

    %%% Linear Mixed Effect
    nanmeanVar = {strcat('nanmean_',var,'~ Group + abs_tgt_rot + (1|SN)')};
    
    lme = fitlme(subj_mean,nanmeanVar{1});
elseif func == 1
    %%% VARIANCE DATA for every subject and every abs(tgt rotation)
    subj_std = varfun(@nanstd, data, 'GroupingVariables',{'SN','Group','abs_tgt_rot'},'OutputFormat','table'); 
    
    subj_std.Group = categorical(subj_std.Group);

    %%% Linear Mixed Effect
    nanstdVar = {strcat('nanstd_',var,'~ Group + abs_tgt_rot + (1|SN)')};
    
    lme = fitlme(subj_std,nanstdVar{1});
end
    
end
%% Load data %% 
clear all

addpath(genpath(pwd))
load('TrainSched_data.mat'); % load data file

%%% REAL DATA
    K_subj_names = {'RBA_K1','RBA_K15','RBC_K6','RBC_K19','RBD_K8','RBD_K22',...
    'RBH_K12','RBH_K24','RBF_K13','RBG_K3','RBG_K17','RBF_K26',... %Regular-Blocked Group
    'RRA_K2','RRA_K16','RRC_K5','RRC_K21','RRD_K7','RRD_K23',...
    'RRH_K10','RRH_K25','RRF_K11','RRF_K28','RRG_K14','RRG_K20',... %Regular-Random Group
    'ERA_K1','ERC_K3','ERG_K2','ERD_K4','ERH_K5','ERF_K6','ERA_K7',...
    'ERG_K8','ERC_K9','ERD_K10','ERH_K12','ERF_K13'...%ErrorClamp-Random Group
    'EBA_K11','EBC_K15','EBD_K16','EBG_K14','EBH_K17','EBF_K18','EBA_K19'...
    'EBC_K21','EBD_K22','EBF_K24','EBG_K20','EBH_K23'}; %ErrorClamp-Blocked Group

    % Check that # of subject names == # of subjects in datafile
    if length(K_subj_names) ~= length(unique(T.SN))
        error_please_copy_in_subject_names
    end

    K_Group_names = {'Regular-Blocked','Regular-Random','ErrorClamp-Blocked','ErrorClamp-Random'};
    K_col = hsv(length(unique(T.Group)));

%%% Set LINES for figures
    lines_all_trials = [15.5 45.5 135.5 225.5 315.5 345.5 435.5 465.5 495.5]; %% current schedule
    lines_rot = [-60 -45 -30 0 30 45 60];

%%% Remove OUTLIERS
    bad_trial = 0;
    bad_handtheta = 0;
    bad_rt = 0;

    for trial_num = 1:length(T.TN)
        if abs(T.hand_theta(trial_num)) > 90 %remove trials with reach angle greater than 90
            T.hand_theta(trial_num) = nan;
            T.RT(trial_num) = nan;
            bad_handtheta = bad_handtheta + 1;
        end
    end
    for trial_num = 1:length(T.TN)
        if abs(T.RT(trial_num)) > 5 %remove trials with RT greater than 5 sec
            T.RT(trial_num) = nan;
            T.hand_theta(trial_num) = nan;
            bad_rt = bad_rt + 1;
        end
    end

    %Show how many outliers of each were exclued
    bad_handtheta
    bad_rt
    
%%% Set recurring VARIABLES
    group = [1 2 3 4];
    rot = [30 45 60];
    target = [30 150 270];
    color = {'b','r','g','m'};
    
%%% Calculate BASELINE SUBTRACTION (D1)
    D1=T;
    baseCN = 7:15; % (6-15) are all the feedback baseline cycles per tgt
    base_idx = D1.CN >= min(baseCN) & D1.CN <= max(baseCN); % baseline index
    
    % Average baseline data for every subject and every tgt
    base_mean = varfun(@nanmean,D1(base_idx ,:),'GroupingVariables',{'SN','ti'},'OutputFormat','table'); 

    for SN = unique(D1.SN)'
        for ti = unique(D1.ti)' % subtract baseline for each target
            trial_idx = (D1.SN==SN & D1.ti==ti);
            base_idx = (base_mean.SN==SN & base_mean.ti==ti);
            
            D1.hand_theta(trial_idx) = D1.hand_theta(trial_idx) - base_mean.nanmean_hand_theta(base_idx);
        end
    end

%%% FLIP baseline subtracted HAND THETA from CCW to CW (D2)
    %Also flip cursor error
    D2 = D1;
    rot_sizes = [45 60 -30; 30 -45 -60; -30 60 45; -60 30 -45; 60 45 -30; -45 -60 30];

    for trial_num = 1:length(D2.TN)
        if D2.tgt_rot(trial_num) > 0
            D2.hand_theta(trial_num) = D2.hand_theta(trial_num)*(-1);
            D2.hand_theta_maxv(trial_num) = D2.hand_theta_maxv(trial_num)*(-1);
            D2.hand_theta_maxradv(trial_num) = D2.hand_theta_maxradv(trial_num)*(-1);
            D2.handMaxRadExt(trial_num) = D2.handMaxRadExt(trial_num)*(-1);
            D2.hand_theta_50(trial_num) = D2.hand_theta_50(trial_num)*(-1);
            D2.CE(trial_num) = D2.CE(trial_num)*(-1);
        end
    end

%% Example Blocked Rotation Schedule %% 
% Use baseline-subtracted, flipped data
% And set which rotatio condition to use
E = D2(D2.rot_cond == 1,:); 
subjs = unique(E.SN)';
i = 1;

%%% Create figure
    Every_sub = figure(); set(gcf,'units','centimeters','pos',[5 5 25 10]); hold on;
    set(gcf,'PaperPositionMode','auto'); set(gca,'FontSize',9);

%%% Select data for subj
    x_data = E.TN (E.SN == subjs(i));
    y_data = E.ri (E.SN == subjs(i));

    % Plot by target
    plot(x_data, y_data,'linewidth',1.5,'color','k');
    ylim([-70 90]);
    xlim([-2 530]);

%%% Shade the no feedback trials
    %Baseline nf
    no_fb_base =patch([0 15.5 15.5 0],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
    set(no_fb_base,'facecolor',.1*[1 1 1]); set(no_fb_base,'edgecolor',.7*[1 1 1]);
    % Aftereffect nf
    no_fb_post =patch([315.5 345.5 345.5 315.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
    set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
    alpha(.1)
    % Retest nf
    no_fb_post =patch([435.5 465.5 465.5 435.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
    set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
    alpha(.1)

%%% Line at 0 (baseline)
    drawline(0, 'dir', 'horz', 'linestyle', '-')

%%% Lines between block
    drawline(lines_all_trials, 'dir', 'vert', 'linestyle', ':');

%%% Title and Axis labels
    str = sprintf('Example Blocked Rotation Schedule');
    title(str,'interpreter','none');
    xlabel('Trial Number') % x-axis label
    ylabel('Target Rotation (deg)') % y-axis label
    set(gca,'FontSize',13);
    
%% Example Random Rotation Schedule %% 
% Use baseline-subtracted, flipped data
% And set which rotatio condition to use
E = D2(D2.rot_cond == 1,:); 
subjs = unique(E.SN)';
i = 3;

%%% Create figure
    Every_sub = figure(); set(gcf,'units','centimeters','pos',[5 5 25 10]); hold on;
    set(gcf,'PaperPositionMode','auto'); set(gca,'FontSize',9);

%%% Select data for subj
    x_data = E.TN (E.SN == subjs(i));
    y_data = E.ri (E.SN == subjs(i));

    % Plot by target
    plot(x_data, y_data,'linewidth',1.5,'color','k');
    ylim([-70 90]);
    xlim([-2 530]);

%%% Shade the no feedback trials
    %Baseline nf
    no_fb_base =patch([0 15.5 15.5 0],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
    set(no_fb_base,'facecolor',.1*[1 1 1]); set(no_fb_base,'edgecolor',.7*[1 1 1]);
    % Aftereffect nf
    no_fb_post =patch([315.5 345.5 345.5 315.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
    set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
    alpha(.1)
    % Retest nf
    no_fb_post =patch([435.5 465.5 465.5 435.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
    set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
    alpha(.1)

%%% Line at 0 (baseline)
    drawline(0, 'dir', 'horz', 'linestyle', '-')

%%% Lines between block
    drawline(lines_all_trials, 'dir', 'vert', 'linestyle', ':');

%%% Title and Axis labels
    str = sprintf('Example Random Rotation Schedule');
    title(str,'interpreter','none');
    xlabel('Trial Number') % x-axis label
    ylabel('Target Rotation (deg)') % y-axis label
    set(gca,'FontSize',13);
    
%% Example Blocked Subject (#11) - Every trial (Hand Angle) split by ABS_TGT_ROT %% 
E=D2; %use baseline-subtracted, flipped data

%%% Every subject
 subjs = unique(E.SN)';
 i=11;

%%% Create figure
 figure(); set(gcf,'units','centimeters','pos',[5 5 25 10]); hold on;
 set(gcf,'PaperPositionMode','auto'); set(gca,'FontSize',9); hold on

%%% Select data for subj
 x_data = E.TN (E.SN == subjs(i));
 y_data = E.hand_theta (E.SN == subjs(i));

%%% Plot by rotation
 scatterplot(x_data, y_data, 'leg','auto','split', E.abs_tgt_rot(E.SN == subjs(i)));
 ylim([-70 90]);
 xlim([-2 530]);

%%% Shade the no feedback trials
 %Baseline nf
 no_fb_base =patch([0 15.5 15.5 0],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_base,'facecolor',.1*[1 1 1]); set(no_fb_base,'edgecolor',.7*[1 1 1]);
 % Aftereffect nf
 no_fb_post =patch([315.5 345.5 345.5 315.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)
 % Retest nf
 no_fb_post =patch([435.5 465.5 465.5 435.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)

%%% Line at hand solution
 %Training block
 plot([45.5 315.5],[30 30],'k','LineWidth',1);
 plot([45.5 315.5],[45 45],'k','LineWidth',1);
 plot([45.5 315.5],[60 60],'k','LineWidth',1);
 %Retest nf block
 plot([435.5 465.5],[30 30],'k','LineWidth',1);
 plot([435.5 465.5],[45 45],'k','LineWidth',1);
 plot([435.5 465.5],[60 60],'k','LineWidth',1);
 %Retest w/f block
 plot([495.5 525.5],[30 30],'k','LineWidth',1);
 plot([495.5 525.5],[45 45],'k','LineWidth',1);
 plot([495.5 525.5],[60 60],'k','LineWidth',1);

%%% Line at 0 (baseline)
 drawline(0, 'dir', 'horz', 'linestyle', '-')

%%% Lines between block
 drawline(lines_all_trials, 'dir', 'vert', 'linestyle', ':');

%%% Title and Axis labels
 str = sprintf('Example Blocked Subject');
 title(str,'interpreter','none');
 xlabel('Trial Number') % x-axis label
 ylabel('Hand Angle (deg)') % y-axis label
 set(gca,'FontSize',13);
    
%% Example Random Subject (#20) - Every trial (Hand Angle) split by ABS_TGT_ROT %% 
 
E=D2; %use baseline-subtracted, flipped data

%%% Every subject
 subjs = unique(E.SN)';
 i=20;

%%% Create figure
 figure(); set(gcf,'units','centimeters','pos',[5 5 25 10]); hold on;
 set(gcf,'PaperPositionMode','auto'); set(gca,'FontSize',9); hold on

%%% Select data for subj
 x_data = E.TN (E.SN == subjs(i));
 y_data = E.hand_theta (E.SN == subjs(i));

%%% Plot by rotation
 scatterplot(x_data, y_data, 'leg','auto','split', E.abs_tgt_rot(E.SN == subjs(i)));
 ylim([-70 90]);
 xlim([-2 530]);

%%% Shade the no feedback trials
 %Baseline nf
 no_fb_base =patch([0 15.5 15.5 0],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_base,'facecolor',.1*[1 1 1]); set(no_fb_base,'edgecolor',.7*[1 1 1]);
 % Aftereffect nf
 no_fb_post =patch([315.5 345.5 345.5 315.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)
 % Retest nf
 no_fb_post =patch([435.5 465.5 465.5 435.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)

%%% Line at hand solution
 %Training block
 plot([45.5 315.5],[30 30],'k','LineWidth',1);
 plot([45.5 315.5],[45 45],'k','LineWidth',1);
 plot([45.5 315.5],[60 60],'k','LineWidth',1);
 %Retest nf block
 plot([435.5 465.5],[30 30],'k','LineWidth',1);
 plot([435.5 465.5],[45 45],'k','LineWidth',1);
 plot([435.5 465.5],[60 60],'k','LineWidth',1);
 %Retest w/f block
 plot([495.5 525.5],[30 30],'k','LineWidth',1);
 plot([495.5 525.5],[45 45],'k','LineWidth',1);
 plot([495.5 525.5],[60 60],'k','LineWidth',1);

%%% Line at 0 (baseline)
 drawline(0, 'dir', 'horz', 'linestyle', '-')

%%% Lines between block
 drawline(lines_all_trials, 'dir', 'vert', 'linestyle', ':');

%%% Title and Axis labels
 str = sprintf('Example Random Subject');
 title(str,'interpreter','none');
 xlabel('Trial Number') % x-axis label
 ylabel('Hand Angle (deg)') % y-axis label
 set(gca,'FontSize',13);
     
%% RT - REG - GROUP averaged EVERY trial split by ROTATION %% 
D2.ReactionTime = D2.RT;

E = D2(D2.Group < 3, :); % Plot only for groups 1 and 2

% Plot Data
dpPlotTrainSched(E, 'ReactionTime', [-1 176 0 1.5])

%% RT - REG - Bar Graph between Groups %% 
% Use baseline-subtracted, flipped data
% And only grab data for Training (3), RNF (6), and RF (8) Blocks
    E=D2(((D2.compBlock == 3 | D2.compBlock == 6 | D2.compBlock == 8)),:); 

%%% AVERAGE DATA %%%
    %for every subject and every abs(tgt rotation)
    subj_mean = varfun(@nanmean,E,'GroupingVariables',{'SN','Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
    %Calculate group MEAN and SEM (by Experiment, Group, abs(tgt rotation))
    grp_mean = varfun(@nanmean,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
%%% RT BAR GRAPH %%%
    xpos = [1 2 3];
    space = [0 .25];

    groups = [1 2];
    cond = [3 6 8];
    cond_name = {'Training', 'Retest No Feedback', 'Retest Feedback'};

    % Loop through groups, conditions, and rotations to plot graph
    for ci = 1:length(cond)
        figure;
        for ri = 1:length(rot)
            for gi = 1:length(groups)
                
                % Set index
                indx = grp_mean.Group == groups(gi) & grp_mean.compBlock == cond(ci)...
                    & grp_mean.abs_tgt_rot == rot(ri);
                
                % Plot group bars
                bar(xpos(ri)+space(gi),grp_mean.nanmean_nanmean_RT(indx),.25,color{gi}); hold on;
                
                % Plot SEM
                plot([xpos(ri)+space(gi) xpos(ri)+space(gi)],...
                    [grp_mean.nanmean_nanmean_RT(indx) - grp_sem.sem_nanmean_RT(indx), ...
                    grp_mean.nanmean_nanmean_RT(indx) + grp_sem.sem_nanmean_RT(indx)], ...
                    'linewidth',2,'color','k');
                
                hold on
                % Plot individual data points
                subj_indx = subj_mean.Group == groups(gi) & subj_mean.compBlock == cond(ci)...
                    & subj_mean.abs_tgt_rot == rot(ri);
                xcoords = repmat(xpos(ri)+space(gi), length(subj_mean.nanmean_RT(subj_indx)),1) +...
                    (randn(length(subj_mean.nanmean_RT(subj_indx)),1)/40) - 0.07; %add jitter
                ycoords = subj_mean.nanmean_RT(subj_indx);
                scatterplot(xcoords,ycoords,'markercolor', 'k','markersize',5, 'markerfill','w'); 

            end
        end
        
        % Set axis limits
        ylims = [0 1.25];
        xlims = [.75 3.5];
        
        % Title and axes labels
        str = sprintf('Reaction Time - %s Block', cond_name{ci});
        title(str,'interpreter','none');
        set(gca,'xtick',[1 1.25 2 2.25 3 3.25],'xticklabel',...
            {'Blocked 30 deg','Random 30 deg','Blocked 45 deg','Random 45 deg','Blocked 60 deg','Random 60 deg'}...
            ,'ylim',ylims,'xlim',xlims,'xticklabelrotation',45)
        ylabel('Reaction Time (sec)') % y-axis label
        set(gca,'FontSize',13);
    end

%%% SIGNIFICANCE TEST %%%
    %Create matrix
    pvalue = nan(3,3); %columns=train, rnf, rf
                       %rows = 30, 45, 60

   for ci = 1:length(cond)
       for ri = 1:length(rot)
           rand_indx = subj_mean.Group == 2 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);
           block_indx = subj_mean.Group == 1 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);

           %Loop through blocks to calculate ttest
           [t,p] = ttest(subj_mean.nanmean_RT(rand_indx),...
               subj_mean.nanmean_RT(block_indx),2,'independent');

           %Place in appropriate cell for matrix
           pvalue(ri,ci) = p;
       end
       
       E=D2( (D2.Group < 3 & D2.compBlock == cond(ci) ) , : );
       lme = lmeTrainSched(E,'RT',0);  
       if ci == 1
           rtTrain = lme.Coefficients;
       elseif ci == 2
           rtRNF = lme.Coefficients;
       elseif ci == 3
           rtRF = lme.Coefficients;
       end
   end

%% LME on RT - REG - for first two cycles of RETEST blocks

% Relevant CN
rnfCN = D2.CN == 146 | D2.CN == 147;
rfCN = D2.CN == 166 | D2.CN == 167;

E = D2((D2.Group < 3 & rnfCN),:);
lme = lmeTrainSched(E,'RT',0);
rtEarlyRNF = lme.Coefficients;

E = D2((D2.Group < 3 & rfCN),:);
lme = lmeTrainSched(E,'RT',0);
rtEarlyRF = lme.Coefficients;

%% RT - REG - Bar Graph for EARLY RNF %% 
% Use baseline-subtracted, flipped data
% And only grab data for Training (3), RNF (6), and RF (8) Blocks
    rnfCN = D2.CN == 146 | D2.CN == 147;
    E = D2((D2.Group < 3 & rnfCN),:);

%%% AVERAGE DATA %%%
    %for every subject and every abs(tgt rotation)
    subj_mean = varfun(@nanmean,E,'GroupingVariables',{'SN','Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
    %Calculate group MEAN and SEM (by Experiment, Group, abs(tgt rotation))
    grp_mean = varfun(@nanmean,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
%%% RT BAR GRAPH %%%
    xpos = [1 2 3];
    space = [0 .25];

    groups = [1 2];

    % Loop through groups, conditions, and rotations to plot graph
        figure;
        for ri = 1:length(rot)
            for gi = 1:length(groups)
                
                % Set index
                indx = grp_mean.Group == groups(gi) & grp_mean.abs_tgt_rot == rot(ri);
                
                % Plot group bars
                bar(xpos(ri)+space(gi),grp_mean.nanmean_nanmean_RT(indx),.25,color{gi}); hold on;
                
                % Plot SEM
                plot([xpos(ri)+space(gi) xpos(ri)+space(gi)],...
                    [grp_mean.nanmean_nanmean_RT(indx) - grp_sem.sem_nanmean_RT(indx), ...
                    grp_mean.nanmean_nanmean_RT(indx) + grp_sem.sem_nanmean_RT(indx)], ...
                    'linewidth',2,'color','k');
                
                hold on
                % Plot individual data points
                subj_indx = subj_mean.Group == groups(gi) & subj_mean.abs_tgt_rot == rot(ri);
                xcoords = repmat(xpos(ri)+space(gi), length(subj_mean.nanmean_RT(subj_indx)),1) +...
                    (randn(length(subj_mean.nanmean_RT(subj_indx)),1)/40) - 0.07; %add jitter
                ycoords = subj_mean.nanmean_RT(subj_indx);
                scatterplot(xcoords,ycoords,'markercolor', 'k','markersize',5, 'markerfill','w'); 

            end
        end
        
        % Set axis limits
        ylims = [0 1.25];
        xlims = [.75 3.5];
        
        % Title and axes labels
        str = sprintf('Reaction Time - Early Retest No Feedback');
        title(str,'interpreter','none');
        set(gca,'xtick',[1 1.25 2 2.25 3 3.25],'xticklabel',...
            {'Blocked 30 deg','Random 30 deg','Blocked 45 deg','Random 45 deg','Blocked 60 deg','Random 60 deg'}...
            ,'ylim',ylims,'xlim',xlims,'xticklabelrotation',45)
        ylabel('Reaction Time (sec)') % y-axis label
        set(gca,'FontSize',13);

%% RT - REG - Bar Graph for EARLY RNF split by TARGET %% 
% Use baseline-subtracted, flipped data
% And only grab data for Training (3), RNF (6), and RF (8) Blocks
    rnfCN = D2.CN == 146;
    E = D2((D2.Group < 3 & rnfCN),:);

%%% AVERAGE DATA %%%
    %for every subject and every abs(tgt rotation)
    subj_mean = varfun(@nanmean,E,'GroupingVariables',{'SN','Group','ti','compBlock'},'OutputFormat','table'); 
    
    %Calculate group MEAN and SEM (by Experiment, Group, abs(tgt rotation))
    grp_mean = varfun(@nanmean,subj_mean,'GroupingVariables',{'Group','ti','compBlock'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,subj_mean,'GroupingVariables',{'Group','ti','compBlock'},'OutputFormat','table'); 
    
%%% RT BAR GRAPH %%%
    xpos = [1 2 3];
    space = [0 .25];

    groups = [1 2];

    % Loop through groups, conditions, and rotations to plot graph
        figure;
        for ti = 1:length(target)
            for gi = 1:length(groups)
                
                % Set index
                indx = grp_mean.Group == groups(gi) & grp_mean.ti == target(ti);
                
                % Plot group bars
                bar(xpos(ti)+space(gi),grp_mean.nanmean_nanmean_RT(indx),.25,color{gi}); hold on;
                
                % Plot SEM
                plot([xpos(ti)+space(gi) xpos(ti)+space(gi)],...
                    [grp_mean.nanmean_nanmean_RT(indx) - grp_sem.sem_nanmean_RT(indx), ...
                    grp_mean.nanmean_nanmean_RT(indx) + grp_sem.sem_nanmean_RT(indx)], ...
                    'linewidth',2,'color','k');
                
                hold on
                % Plot individual data points
                subj_indx = subj_mean.Group == groups(gi) & subj_mean.ti == target(ti);
                xcoords = repmat(xpos(ti)+space(gi), length(subj_mean.nanmean_RT(subj_indx)),1) +...
                    (randn(length(subj_mean.nanmean_RT(subj_indx)),1)/40) - 0.07; %add jitter
                ycoords = subj_mean.nanmean_RT(subj_indx);
                scatterplot(xcoords,ycoords,'markercolor', 'k','markersize',5, 'markerfill','w'); 

            end
        end
        
        % Set axis limits
        ylims = [0 1.25];
        xlims = [.75 3.5];
        
        % Title and axes labels
        str = sprintf('Reaction Time - Early Retest No Feedback');
        title(str,'interpreter','none');
        set(gca,'xtick',[1 1.25 2 2.25 3 3.25],'xticklabel',...
            {'Blocked 30 deg','Random 30 deg','Blocked 150 deg','Random 150 deg','Blocked 270 deg','Random 270 deg'}...
            ,'ylim',ylims,'xlim',xlims,'xticklabelrotation',45)
        ylabel('Reaction Time (sec)') % y-axis label
        set(gca,'FontSize',13);

%%% Linear Mixed Effect %%%
    subj_mean.Group = categorical(subj_mean.Group);
    subj_mean.ti = categorical(subj_mean.ti);
    lme = fitlme(subj_mean,'nanmean_RT ~ Group + ti + (1|SN)');
    RTgtEarlyRNF = lme.Coefficients;

%% Hand Angle - REG - GROUP averaged EVERY trial split by ROTATION %% 
D2.HandAngle = D2.hand_theta;

E = D2(D2.Group < 3, :); % Plot only for groups 1 and 2

% Plot Data
dpPlotTrainSched(E, 'HandAngle', [-1 176 -10 75])

%% Hand Angle - REG - Bar Graph between Groups %% 
% Use baseline-subtracted, flipped data
% Only grab data for Training (3), RNF (6), and RF (8) Blocks
E=D2(D2.compBlock == 3 | D2.compBlock == 6 | D2.compBlock == 8,:); 

%%% AVERAGE DATA %%% 
    % For every subject and every abs(tgt rotation)
    subj_mean = varfun(@nanmean,E,'GroupingVariables',{'SN','Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
    % Calculate group MEAN and SEM (by Experiment, Group, abs(tgt rotation))
    grp_mean = varfun(@nanmean,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
%%% RT BAR GRAPH %%%
    xpos = [1 2 3];
    space = [0 .25];

    groups = [1 2]; %only data for experiment 1
    cond = [3 6 8];
    cond_name = {'Training', 'Retest No Feedback', 'Retest Feedback'};

    % Loop through groups, conditions, and rotations to plot graph
    for ci = 1:length(cond)
        figure;
        for ri = 1:length(rot)
            
            for gi = 1:length(groups)
                
                % Set index
                indx = grp_mean.Group == groups(gi) & grp_mean.compBlock == cond(ci)...
                    & grp_mean.abs_tgt_rot == rot(ri);
                
                % Plot group bars
                RT_bar = bar(xpos(ri)+space(gi),grp_mean.nanmean_nanmean_hand_theta(indx),.25,color{gi}); hold on;
                
                % Plot SEM
                plot([xpos(ri)+space(gi) xpos(ri)+space(gi)],...
                    [grp_mean.nanmean_nanmean_hand_theta(indx) - grp_sem.sem_nanmean_hand_theta(indx), ...
                    grp_mean.nanmean_nanmean_hand_theta(indx) + grp_sem.sem_nanmean_hand_theta(indx)], ...
                    'linewidth',2,'color','k');
                
                hold on
                % Plot individual data points
                subj_indx = subj_mean.Group == groups(gi) & subj_mean.compBlock == cond(ci)...
                    & subj_mean.abs_tgt_rot == rot(ri);
                xcoords = repmat(xpos(ri)+space(gi), length(subj_mean.nanmean_hand_theta(subj_indx)),1) +...
                    (randn(length(subj_mean.nanmean_hand_theta(subj_indx)),1)/40) - 0.07; %add jitter
                ycoords = subj_mean.nanmean_hand_theta(subj_indx);
                scatterplot(xcoords,ycoords,'markercolor', 'k','markersize',5, 'markerfill','w'); 
                                
            end
        end
        
        % Set axis limits
        ylims = [0 80];
        xlims = [.75 3.5];
        
        % Title and axes labels
        str = sprintf('Hand Angle - %s Block', cond_name{ci});
        title(str,'interpreter','none');
        set(gca,'xtick',[1 1.25 2 2.25 3 3.25],'xticklabel',...
            {'Blocked 30 deg','Random 30 deg','Blocked 45 deg','Random 45 deg','Blocked 60 deg','Random 60 deg'}...
            ,'ylim',ylims,'xlim',xlims,'xticklabelrotation',45)
        ylabel('Hand Angle (deg)') % y-axis label
        set(gca,'FontSize',13);
    end

%%% SIGNIFICANCE TEST %%%
    %Create matrix
    pvalue = nan(3,3); %columns=train, rnf, rf
                       %rows = 30, 45, 60

   for ci = 1:length(cond)
       for ri = 1:length(rot)
           rand_indx = subj_mean.Group == 2 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);
           block_indx = subj_mean.Group == 1 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);

           %Loop through blocks to calculate ttest
           [t,p] = ttest(subj_mean.nanmean_hand_theta(rand_indx),...
               subj_mean.nanmean_hand_theta(block_indx),2,'independent');

           %Place in appropriate cell for matrix
           pvalue(ri,ci) = p;
       end
       
       E=D2( (D2.Group < 3 & D2.compBlock == cond(ci) ) , : );
       lme = lmeTrainSched(E,'hand_theta',0);  
       if ci == 1
           handTrain = lme.Coefficients;
       elseif ci == 2
           handRNF = lme.Coefficients;
       elseif ci == 3
           handRF = lme.Coefficients;
       end
   end

%% LME on Hand Angle - REG - for first two cycles of RETEST blocks

% Relevant CN
rnfCN = D2.CN == 146 | D2.CN == 147;
rfCN = D2.CN == 166 | D2.CN == 167;

E = D2((D2.Group < 3 & rnfCN),:);
lme = lmeTrainSched(E,'hand_theta',0);
handEarlyRNF = lme.Coefficients;

E = D2((D2.Group < 3 & rfCN),:);
lme = lmeTrainSched(E,'hand_theta',0);
handEarlyRF = lme.Coefficients;   

%% Train Variance - REG - GROUP averaged EVERY trial split by ROTATION %% 
%use baseline-subtracted, flipped data
    E =D2; 
    E.hand = E.hand_theta; % Set which hand variable to look at
    E = E(E.Group < 3, :); % Plot only for groups 1 and 2

%%% CREATE TRAINING DATA %%%
    varianceData = table; % define table
    rotations = unique(E.abs_tgt_rot)';
    train_indx = E.TN_tgt >= 16 & E.TN_tgt <= 105; % training index (Trial number per target)

    for sn = unique(E.SN)'
        for rotAng = rotations

            % empty variables for each loop
            tgt_movmean=[];SN=[];Group=[];TN_tgt=[];tgtRot=[];

            % Calculate moving STD for each target
            idx = (E.SN == sn & E.abs_tgt_rot == rotAng & train_indx) ; % training trials idx for specific subj & target

            handAng = E.hand(idx);
            residuals = handAng - movmean(handAng,10,'omitnan'); %calculate residuals on moving mean (window=10)

            movStd = movstd(residuals,10,'omitnan');


            % Save dependent variables for table
            SN = repmat(sn, length(movStd), 1);
            Group = repmat( unique(E.Group(idx)), length(movStd), 1);
            TN_tgt = E.TN_tgt(idx);
            tgtRot = repmat(rotAng, length(movStd), 1);

            % Append data to overall table
            varianceData = [varianceData; table(SN, Group, TN_tgt, tgtRot, movStd)];
        end
    end

%%% CREATE GROUP MEAN AND SEM %%%
    meanData = varfun(@nanmean, varianceData, 'GroupingVariables', {'Group','tgtRot','TN_tgt'},'OutputFormat','table');
    semData = varfun(@sem, varianceData, 'GroupingVariables', {'Group','tgtRot','TN_tgt'},'OutputFormat','table');

%%% PLOT DATA %%%
    figure;
    for ri = 1:length(rotations)
        for gi = unique(E.Group)'
            subplot(3,1,ri); hold on;

            idx = meanData.Group == gi & meanData.tgtRot == rotations(ri);

            % Plot data
            x = meanData.TN_tgt(idx);
            y = meanData.nanmean_movStd(idx);
            err = semData.sem_movStd(idx);
            shadedErrorBar(x, y, err, color{gi}, 1);

            % Title and axes labels
            axis([min(x)-1 max(x)+1 0 45])
            str = sprintf('%iº Rotation, Moving StD (Window=10)', rotations(ri));
            title(str, 'interpreter','none')
            xlabel('Cycle') % x-axis label
            ylabel('Standard Deviation') % y-axis label
            set(gca,'FontSize',11);

        end
    end

%% Variance - REG - Bar Graph between Groups %% 
E=D2; %use baseline-subtracted, flipped data

%%% AVERAGE DATA %%%
    % STD for every subject and every abs(tgt rotation)
    subj_std = varfun(@nanstd,E,'GroupingVariables',{'SN','Group','abs_tgt_rot','blockNum'},'OutputFormat','table'); 
    
    % Calculate group MEAN and SEM (by Experiment, Group, abs(tgt rotation))
    grp_mean = varfun(@nanmean,subj_std,'GroupingVariables',{'Group','abs_tgt_rot','blockNum'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,subj_std,'GroupingVariables',{'Group','abs_tgt_rot','blockNum'},'OutputFormat','table'); 
    
 %%% BAR GRAPH %%%
    xpos = [0 1 2];
    space = [0 .25];

    groups = [1 2]; %experiment 1 only
    cond = [3 6 8]; %train & rnf & rf only
    cond_names = {'Training','Retest No Feedback','Retest w/Feedback'};
    rot = [30 45 60];

    % Loop through groups, conditions, and rotations to plot graph
    for ci = 1:length(cond)
        figure;
        
        for gi = 1:length(groups)
            for ri = 1:length(rot)
                
                % Set index
                indx = grp_mean.Group == groups(gi) & grp_mean.blockNum == cond(ci)...
                    & grp_mean.abs_tgt_rot == rot(ri);
                
                
                % Plot group bars
                bar(xpos(ri)+space(gi),grp_mean.nanmean_nanstd_hand_theta(indx),.25,color{gi}); 
                hold on;
                
                % Plot SEM
                plot([xpos(ri)+space(gi) xpos(ri)+space(gi)],...
                    [grp_mean.nanmean_nanstd_hand_theta(indx) - grp_sem.sem_nanstd_hand_theta(indx), ...
                    grp_mean.nanmean_nanstd_hand_theta(indx) + grp_sem.sem_nanstd_hand_theta(indx)], ...
                    'linewidth',2,'color','k'); 
                
                hold on;
                % Plot individual data points
                subj_indx = subj_std.Group == groups(gi) & subj_std.blockNum == cond(ci)...
                    & subj_std.abs_tgt_rot == rot(ri);
                xcoords = repmat(xpos(ri)+space(gi), length(subj_std.nanstd_hand_theta(subj_indx)),1) +...
                    (randn(length(subj_std.nanstd_hand_theta(subj_indx)),1)/40) - 0.07; %add jitter
                ycoords = subj_std.nanstd_hand_theta(subj_indx);
                scatterplot(xcoords,ycoords,'markercolor', 'k','markersize',5, 'markerfill','w');
                
                
            end
        end
        
        % Set axis limits
        ylims = [0 27];
        xlims = [-.5 3];
        
        % Title and axes labels
        str = sprintf('Standard Deviation - %s Block', cond_names{ci});
        title(str, 'interpreter','none')        
        set(gca,'xtick',[-.05 .25 .95 1.25 1.95 2.25],'xticklabel',...
            {'Blocked 30 deg','Random 30 deg','Blocked 45 deg','Random 45 deg','Blocked 60 deg','Random 60 deg'}...
            ,'ylim',ylims,'xlim',xlims,'xticklabelrotation',45)
        ylabel('Standard Deviation') % y-axis label
        set(gca,'FontSize',13);
    end

%%% SIGNIFICANCE TEST %%%
    %Create matrix
    pvalue = nan(3,3); %columns= train, rnf, rf
                       %rows = 30, 45, 60

    for ci = 1:length(cond)
       for ri = 1:length(rot)
           rand_indx = subj_std.Group == 2 & subj_std.blockNum == cond(ci) ...
               & subj_std.abs_tgt_rot == rot(ri);
           block_indx = subj_std.Group == 1 & subj_std.blockNum == cond(ci) ...
               & subj_std.abs_tgt_rot == rot(ri);

           %Loop through blocks to calculate ttest
           [t,p] = ttest(subj_std.nanstd_hand_theta(rand_indx),...
               subj_std.nanstd_hand_theta(block_indx),2,'independent');

           %Place in appropriate cell for matrix
           pvalue(ri,ci) = p;
       end
       
       E=D2( (D2.Group < 3 & D2.compBlock == cond(ci) ) , : );
       lme = lmeTrainSched(E,'hand_theta',1);  
       if ci == 1
           varTrain = lme.Coefficients;
       elseif ci == 2
           varRNF = lme.Coefficients;
       elseif ci == 3
           varRF = lme.Coefficients;
       end
    end

%% LME on Variance - REG - for first two cycles of RETEST blocks

% Relevant CN
rnfCN = D2.CN == 146 | D2.CN == 147;
rfCN = D2.CN == 166 | D2.CN == 167;

E = D2((D2.Group < 3 & rnfCN),:);
lme = lmeTrainSched(E,'hand_theta',1);
varEarlyRNF = lme.Coefficients;

E = D2((D2.Group < 3 & rfCN),:);
lme = lmeTrainSched(E,'hand_theta',1);
varEarlyRF = lme.Coefficients;

%% Aftereffect - REG %% 
E=D2; %use baseline-subtracted, flipped data
    
%%% Set variables %%%
    ae_idx = E.CN >= 106 & E.CN <= 115; % aftereffect index
    group = [1 2]; %1 = reg-blocked, 2 = reg-random
    color = {'b','r'};
    
%%% Calculate group MEAN and SEM (by Group, CN) %%%
    grp_mean = varfun(@nanmean,E(ae_idx ,:),'GroupingVariables',...
        {'Group','CN'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,E(ae_idx ,:),'GroupingVariables',...
        {'Group','CN'},'OutputFormat','table');
 
%%% PLOT FIGURE %%%
figure;
    for gi = 1:length(group)
        x = unique(grp_mean.CN);
        y = grp_mean.nanmean_hand_theta(grp_mean.Group == group(gi)); 
        err = grp_sem.sem_hand_theta(grp_sem.Group == group(gi));
        
        plot(x,y,'-','markersize',15,'color',color{gi}); hold all;
        shadedErrorBar(x, y, err, color{gi}, 1); hold all
    end
    
    % Set axes
    xlim([105 116]);
    ylim([0 25]);

    % Draw lines
    drawline(x, 'dir', 'vert', 'linestyle', ':'); %target angle
    drawline(0, 'dir', 'hor', 'linestyle', '-'); %0º hand angle

    % Title and axes labels
    title('Average Aftereffect');
    xlabel('Cycle Number') % x-axis label
    ylabel('Hand Angle (deg)') % y-axis label
    set(gca,'FontSize',13);

%%% MIXED EFFECTS MODEL %%%
    E=D2( (D2.Group < 3 & D2.CN >= 106 & D2.CN <= 110 ) , : ); %select only first few AE cycles
    lme = lmeTrainSched(E,'hand_theta',0);
    ae = lme.Coefficients;

%% Aftereffect - REG - Split by TGT %%
E=D2; %use baseline-subtracted, flipped data

%%% Set variables %%%
    ae_idx = E.CN >= 106 & E.CN <= 115; % aftereffect index
    group = [1 2]; %1 = reg-blocked, 2 = reg-random

%Calculate group MEAN and SEM by Group, CN, target
grp_mean = varfun(@nanmean,E(ae_idx ,:),'GroupingVariables',...
    {'Group','CN','ti'},'OutputFormat','table');
grp_sem = varfun(@sem,E(ae_idx ,:),'GroupingVariables',...
    {'Group','CN','ti'},'OutputFormat','table');

color = {'b','r','g'};

figure;
for gi = 1:length(group)
    for ti = 1:length(target)

        x = unique(grp_mean.CN);
        y = grp_mean.nanmean_hand_theta(grp_mean.Group == group(gi)...
            & grp_mean.ti == target(ti));
        err = grp_sem.sem_hand_theta(grp_sem.Group == group(gi) &...
            grp_sem.ti == target(ti));

        plot(x,y,'-','markersize',15,'color',color{ti}); hold all;
        shadedErrorBar(x, y, err, color{ti},1); hold all
    end
end

% Set axes
xlim([105 116]);
ylim([5 25]);

% Draw lines
drawline(x, 'dir', 'vert', 'linestyle', ':'); %target angle
drawline(0, 'dir', 'hor', 'linestyle', '-'); %0º hand angle

% Title and axes labels
title('Average Aftereffect Split By Target');
xlabel('Cycle Number') % x-axis label
ylabel('Hand Angle (deg)') % y-axis label
set(gca,'FontSize',13);

%% Example EC Blocked Subject (#39) - Every trial (Hand Angle) split by ABS_TGT_ROT %% 
E=D2; %use baseline-subtracted, flipped data

%%% Every subject
 subjs = unique(E.SN)';
 i=39;

%%% Create figure
 figure(); set(gcf,'units','centimeters','pos',[5 5 25 10]); hold on;
 set(gcf,'PaperPositionMode','auto'); set(gca,'FontSize',9); hold on

%%% Select data for subj
 x_data = E.TN (E.SN == subjs(i));
 y_data = E.hand_theta (E.SN == subjs(i));

%%% Plot by rotation
 scatterplot(x_data, y_data, 'leg','auto','split', E.abs_tgt_rot(E.SN == subjs(i)));
 ylim([-70 90]);
 xlim([-2 530]);

%%% Shade the no feedback trials
 %Baseline nf
 no_fb_base =patch([0 15.5 15.5 0],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_base,'facecolor',.1*[1 1 1]); set(no_fb_base,'edgecolor',.7*[1 1 1]);
 % Aftereffect nf
 no_fb_post =patch([315.5 345.5 345.5 315.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)
 % Retest nf
 no_fb_post =patch([435.5 465.5 465.5 435.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)

%%% Line at hand solution
 %Training block
 plot([45.5 315.5],[30 30],'k','LineWidth',1);
 plot([45.5 315.5],[45 45],'k','LineWidth',1);
 plot([45.5 315.5],[60 60],'k','LineWidth',1);
 %Retest nf block
 plot([435.5 465.5],[30 30],'k','LineWidth',1);
 plot([435.5 465.5],[45 45],'k','LineWidth',1);
 plot([435.5 465.5],[60 60],'k','LineWidth',1);
 %Retest w/f block
 plot([495.5 525.5],[30 30],'k','LineWidth',1);
 plot([495.5 525.5],[45 45],'k','LineWidth',1);
 plot([495.5 525.5],[60 60],'k','LineWidth',1);

%%% Line at 0 (baseline)
 drawline(0, 'dir', 'horz', 'linestyle', '-')

%%% Lines between block
 drawline(lines_all_trials, 'dir', 'vert', 'linestyle', ':');

%%% Title and Axis labels
 str = sprintf('Example Error-Clamp Blocked Subject');
 title(str,'interpreter','none');
 xlabel('Trial number') % x-axis label
 ylabel('Hand angle (deg)') % y-axis label
 set(gca,'FontSize',13);

%% Example EC Random Subject (#33) - Every trial (Hand Angle) split by ABS_TGT_ROT %% 
E=D2; %use baseline-subtracted, flipped data

%%% Every subject
 subjs = unique(E.SN)';
 i=33;

%%% Create figure
 figure(); set(gcf,'units','centimeters','pos',[5 5 25 10]); hold on;
 set(gcf,'PaperPositionMode','auto'); set(gca,'FontSize',9); hold on

%%% Select data for subj
 x_data = E.TN (E.SN == subjs(i));
 y_data = E.hand_theta (E.SN == subjs(i));

%%% Plot by rotation
 scatterplot(x_data, y_data, 'leg','auto','split', E.abs_tgt_rot(E.SN == subjs(i)));
 ylim([-70 90]);
 xlim([-2 530]);

%%% Shade the no feedback trials
 %Baseline nf
 no_fb_base =patch([0 15.5 15.5 0],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_base,'facecolor',.1*[1 1 1]); set(no_fb_base,'edgecolor',.7*[1 1 1]);
 % Aftereffect nf
 no_fb_post =patch([315.5 345.5 345.5 315.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)
 % Retest nf
 no_fb_post =patch([435.5 465.5 465.5 435.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
 set(no_fb_post,'facecolor',.1*[1 1 1]); set(no_fb_post,'edgecolor',.7*[1 1 1]);
 alpha(.1)

%%% Line at hand solution
 %Training block
 plot([45.5 315.5],[30 30],'k','LineWidth',1);
 plot([45.5 315.5],[45 45],'k','LineWidth',1);
 plot([45.5 315.5],[60 60],'k','LineWidth',1);
 %Retest nf block
 plot([435.5 465.5],[30 30],'k','LineWidth',1);
 plot([435.5 465.5],[45 45],'k','LineWidth',1);
 plot([435.5 465.5],[60 60],'k','LineWidth',1);
 %Retest w/f block
 plot([495.5 525.5],[30 30],'k','LineWidth',1);
 plot([495.5 525.5],[45 45],'k','LineWidth',1);
 plot([495.5 525.5],[60 60],'k','LineWidth',1);

%%% Line at 0 (baseline)
 drawline(0, 'dir', 'horz', 'linestyle', '-')

%%% Lines between block
 drawline(lines_all_trials, 'dir', 'vert', 'linestyle', ':');

%%% Title and Axis labels
 str = sprintf('Example Error-Clamp Random Subject');
 title(str,'interpreter','none');
 xlabel('Trial number') % x-axis label
 ylabel('Hand angle (deg)') % y-axis label
 set(gca,'FontSize',13);
     
%% RT - EC - GROUP averaged EVERY trial split by ROTATION %% 
D2.ReactionTime = D2.RT;

E = D2(D2.Group > 2, :); % Plot only for groups 1 and 2

% Plot Data
dpPlotTrainSched(E, 'ReactionTime', [-1 176 0 1.5])

%% RT - EC - Bar Graph between Groups %% 
% Use baseline-subtracted, flipped data
% Only grab data for Training (3), RNF (6), and RF (8) Blocks
    E=D2(((D2.compBlock == 3 | D2.compBlock == 6 | D2.compBlock == 8)),:); 

%%% AVERAGE DATA %%% 
    % For every subject and every abs(tgt rotation)
    subj_mean = varfun(@nanmean,E,'GroupingVariables',{'SN','Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
    % Calculate group MEAN and SEM (by Experiment, Group, abs(tgt rotation))
    grp_mean = varfun(@nanmean,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
%%% RT BAR GRAPH %%%
    xpos = [1 2 3];
    space = [0 .25];

    groups = [3 4];
    cond = [3 6 8];
    color = {'g','m'};
    cond_name = {'Training', 'Retest No Feedback', 'Retest Feedback'};

    % Loop through groups, conditions, and rotations to plot graph
    for ci = 1:length(cond)
        figure;
        for ri = 1:length(rot)
            for gi = 1:length(groups)
                
                % Set index
                indx = grp_mean.Group == groups(gi) & grp_mean.compBlock == cond(ci)...
                    & grp_mean.abs_tgt_rot == rot(ri);
                
                % Plot group bars
                bar(xpos(ri)+space(gi),grp_mean.nanmean_nanmean_RT(indx),.25,color{gi}); hold on;
                
                % Plot SEM
                plot([xpos(ri)+space(gi) xpos(ri)+space(gi)],...
                    [grp_mean.nanmean_nanmean_RT(indx) - grp_sem.sem_nanmean_RT(indx), ...
                    grp_mean.nanmean_nanmean_RT(indx) + grp_sem.sem_nanmean_RT(indx)], ...
                    'linewidth',2,'color','k');
                
                hold on
                % Plot individual data points
                subj_indx = subj_mean.Group == groups(gi) & subj_mean.compBlock == cond(ci)...
                    & subj_mean.abs_tgt_rot == rot(ri);
                xcoords = repmat(xpos(ri)+space(gi), length(subj_mean.nanmean_RT(subj_indx)),1) +...
                    (randn(length(subj_mean.nanmean_RT(subj_indx)),1)/40) - 0.07; %add jitter
                ycoords = subj_mean.nanmean_RT(subj_indx);
                scatterplot(xcoords,ycoords,'markercolor', 'k','markersize',5, 'markerfill','w'); 

            end
        end
        
        % Set axis limits
        ylims = [0 1.25];
        xlims = [.75 3.5];
        
        % Title and axes labels
        str = sprintf('Reaction Time - %s Block', cond_name{ci});
        title(str,'interpreter','none');
        set(gca,'xtick',[1 1.25 2 2.25 3 3.25],'xticklabel',...
            {'Blocked 30 deg','Random 30 deg','Blocked 45 deg','Random 45 deg','Blocked 60 deg','Random 60 deg'}...
            ,'ylim',ylims,'xlim',xlims,'xticklabelrotation',45)
        ylabel('Reaction Time (sec)') % y-axis label
        set(gca,'FontSize',13);
    end

%%% SIGNIFICANCE TEST %%%
    %Create matrix
    pvalue = nan(3,3); %columns=train, rnf, rf
                       %rows = 30, 45, 60

   for ci = 1:length(cond)
       for ri = 1:length(rot)
           rand_indx = subj_mean.Group == 2 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);
           block_indx = subj_mean.Group == 1 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);

           %Loop through blocks to calculate ttest
           [t,p] = ttest(subj_mean.nanmean_RT(rand_indx),...
               subj_mean.nanmean_RT(block_indx),2,'independent');

           %Place in appropriate cell for matrix
           pvalue(ri,ci) = p;
       end
       
       E=D2( (D2.Group > 2 & D2.compBlock == cond(ci) ) , : );
       lme = lmeTrainSched(E,'RT',0);  
       if ci == 1
           rtTrain = lme.Coefficients;
       elseif ci == 2
           rtRNF = lme.Coefficients;
       elseif ci == 3
           rtRF = lme.Coefficients;
       end
   end

%% Hand Angle - EC - GROUP averaged EVERY trial split by ROTATION %% 
D2.HandAngle = D2.hand_theta;

E = D2(D2.Group > 2, :); % Plot only for groups 1 and 2

% Plot Data
dpPlotTrainSched(E, 'HandAngle', [-1 176 -10 75])

%% Hand Angle - EC - Bar Graph between Groups %% 
% Use baseline-subtracted, flipped data
% Only grab data for Training (3), RNF (6), and RF (8) Blocks
E=D2(D2.compBlock == 3 | D2.compBlock == 6 | D2.compBlock == 8,:); 

%%% AVERAGE DATA %%% 
    % For every subject and every abs(tgt rotation)
    subj_mean = varfun(@nanmean,E,'GroupingVariables',{'SN','Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
    % Calculate group MEAN and SEM (by Experiment, Group, abs(tgt rotation))
    grp_mean = varfun(@nanmean,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,subj_mean,'GroupingVariables',{'Group','abs_tgt_rot','compBlock'},'OutputFormat','table'); 
    
%%% RT BAR GRAPH %%%
    xpos = [1 2 3];
    space = [0 .25];

    groups = [3 4];
    color = {'g', 'm'};
    cond = [3 6 8];
    cond_name = {'Training', 'Retest No Feedback', 'Retest Feedback'};

    % Loop through groups, conditions, and rotations to plot graph
    for ci = 1:length(cond)
        figure;
        for ri = 1:length(rot)
            
            for gi = 1:length(groups)
                
                % Set index
                indx = grp_mean.Group == groups(gi) & grp_mean.compBlock == cond(ci)...
                    & grp_mean.abs_tgt_rot == rot(ri);
                
                % Plot group bars
                RT_bar = bar(xpos(ri)+space(gi),grp_mean.nanmean_nanmean_hand_theta(indx),.25,color{gi}); hold on;
                
                % Plot SEM
                plot([xpos(ri)+space(gi) xpos(ri)+space(gi)],...
                    [grp_mean.nanmean_nanmean_hand_theta(indx) - grp_sem.sem_nanmean_hand_theta(indx), ...
                    grp_mean.nanmean_nanmean_hand_theta(indx) + grp_sem.sem_nanmean_hand_theta(indx)], ...
                    'linewidth',2,'color','k');
                
                hold on
                % Plot individual data points
                subj_indx = subj_mean.Group == groups(gi) & subj_mean.compBlock == cond(ci)...
                    & subj_mean.abs_tgt_rot == rot(ri);
                xcoords = repmat(xpos(ri)+space(gi), length(subj_mean.nanmean_hand_theta(subj_indx)),1) +...
                    (randn(length(subj_mean.nanmean_hand_theta(subj_indx)),1)/40) - 0.07; %add jitter
                ycoords = subj_mean.nanmean_hand_theta(subj_indx);
                scatterplot(xcoords,ycoords,'markercolor', 'k','markersize',5, 'markerfill','w'); 
                                
            end
        end
        
        % Set axis limits
        ylims = [0 80];
        xlims = [.75 3.5];
        
        % Title and axes labels
        str = sprintf('Hand Angle - %s Block', cond_name{ci});
        title(str,'interpreter','none');
        set(gca,'xtick',[1 1.25 2 2.25 3 3.25],'xticklabel',...
            {'Blocked 30 deg','Random 30 deg','Blocked 45 deg','Random 45 deg','Blocked 60 deg','Random 60 deg'}...
            ,'ylim',ylims,'xlim',xlims,'xticklabelrotation',45)
        ylabel('Hand Angle (deg)') % y-axis label
        set(gca,'FontSize',13);
    end

%%% SIGNIFICANCE TEST %%%
    %Create matrix
    pvalue = nan(3,3); %columns=train, rnf, rf
                       %rows = 30, 45, 60

   for ci = 1:length(cond)
       for ri = 1:length(rot)
           rand_indx = subj_mean.Group == 2 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);
           block_indx = subj_mean.Group == 1 & subj_mean.compBlock == cond(ci) ...
               & subj_mean.abs_tgt_rot == rot(ri);

           %Loop through blocks to calculate ttest
           [t,p] = ttest(subj_mean.nanmean_hand_theta(rand_indx),...
               subj_mean.nanmean_hand_theta(block_indx),2,'independent');

           %Place in appropriate cell for matrix
           pvalue(ri,ci) = p;
       end
       
       E=D2( (D2.Group > 2 & D2.compBlock == cond(ci) ) , : );
       lme = lmeTrainSched(E,'hand_theta',0);  
       if ci == 1
           handTrain = lme.Coefficients;
       elseif ci == 2
           handRNF = lme.Coefficients;
       elseif ci == 3
           handRF = lme.Coefficients;
       end
   end

%% Aftereffect - EC %% 

E=D2; %use baseline-subtracted, flipped data
    
%%% Set variables %%%
    ae_idx = E.CN >= 106 & E.CN <= 115; % aftereffect index
    group = [3 4]; %1 = reg-blocked, 2 = reg-random
    color = {'g','m'};
    
%%% Calculate group MEAN and SEM (by Group, CN) %%%
    grp_mean = varfun(@nanmean,E(ae_idx ,:),'GroupingVariables',...
        {'Group','CN'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,E(ae_idx ,:),'GroupingVariables',...
        {'Group','CN'},'OutputFormat','table');
 
%%% Plot figure %%%
figure;
    for gi = 1:length(group)
        x = unique(grp_mean.CN);
        y = grp_mean.nanmean_hand_theta(grp_mean.Group == group(gi)); 
        err = grp_sem.sem_hand_theta(grp_sem.Group == group(gi));
        
        plot(x,y,'-','markersize',15,'color',color{gi}); hold all;
        shadedErrorBar(x, y, err, color{gi}, 1); hold all
    end
    
    % Set axes
    xlim([105 116]);
    ylim([0 25]);

    % Draw lines
    drawline(x, 'dir', 'vert', 'linestyle', ':'); %target angle
    drawline(0, 'dir', 'hor', 'linestyle', '-'); %0º hand angle

    % Title and axes labels
    title('Average Aftereffect (Error-Clamp Experiment)');
    xlabel('Cycle Number') % x-axis label
    ylabel('Hand Angle (deg)') % y-axis label
    set(gca,'FontSize',13);
    
%%% MIXED EFFECTS MODEL %%%
    E=D2( (D2.Group > 2 & D2.CN >= 106 & D2.CN <= 110 ) , : ); %select only first few AE cycles
    lme = lmeTrainSched(E,'hand_theta',0);
    ae = lme.Coefficients;

%% Aftereffect - EC - Split by TGT %%

%% Aftereffect - REG - Split by TGT %%
E=D2; %use baseline-subtracted, flipped data

%%% Set variables %%%
    ae_idx = E.CN >= 106 & E.CN <= 115; % aftereffect index
    group = [3 4]; %3 = ec-blocked, 4 = ec-random

%Calculate group MEAN and SEM by Group, CN, target
grp_mean = varfun(@nanmean,E(ae_idx ,:),'GroupingVariables',...
    {'Group','CN','ti'},'OutputFormat','table');
grp_sem = varfun(@sem,E(ae_idx ,:),'GroupingVariables',...
    {'Group','CN','ti'},'OutputFormat','table');

color = {'b','r','g'};

figure;
for gi = 1
    for ti = 1:length(target)

        x = unique(grp_mean.CN);
        y = grp_mean.nanmean_hand_theta(grp_mean.Group == group(gi)...
            & grp_mean.ti == target(ti));
        err = grp_sem.sem_hand_theta(grp_sem.Group == group(gi) &...
            grp_sem.ti == target(ti));

        plot(x,y,'-','markersize',15,'color',color{ti}); hold all;
        shadedErrorBar(x, y, err, color{ti},1); hold all
    end
end

% Set axes
xlim([105 116]);
ylim([5 25]);

% Draw lines
drawline(x, 'dir', 'vert', 'linestyle', ':'); %target angle
drawline(0, 'dir', 'hor', 'linestyle', '-'); %0º hand angle

% Title and axes labels
title('Average Aftereffect Split By Target');
xlabel('Cycle Number') % x-axis label
ylabel('Hand Angle (deg)') % y-axis label
set(gca,'FontSize',13);

%% Aftereffect - All Groups %% 

E=D2; %use baseline-subtracted, flipped data
    
%%% Set variables %%%
    ae_idx = E.CN >= 106 & E.CN <= 115; % aftereffect index
    group = [1 2 3 4]; %1 = reg-blocked, 2 = reg-random
    color = {'b','r','g','m'};
    
%%% Calculate group MEAN and SEM (by Group, CN) %%%
    grp_mean = varfun(@nanmean,E(ae_idx ,:),'GroupingVariables',...
        {'Group','CN'},'OutputFormat','table'); 
    grp_sem = varfun(@sem,E(ae_idx ,:),'GroupingVariables',...
        {'Group','CN'},'OutputFormat','table');
 
%%% Plot figure %%%
figure;
    for gi = 1:length(group)
        x = unique(grp_mean.CN);
        y = grp_mean.nanmean_hand_theta(grp_mean.Group == group(gi)); 
        err = grp_sem.sem_hand_theta(grp_sem.Group == group(gi));
        
        plot(x,y,'-','markersize',15,'color',color{gi}); hold all;
        shadedErrorBar(x, y, err, color{gi}, 1); hold all
    end
    
    % Set axes
    xlim([105 116]);
    ylim([0 25]);

    % Draw lines
    drawline(x, 'dir', 'vert', 'linestyle', ':'); %target angle
    drawline(0, 'dir', 'hor', 'linestyle', '-'); %0º hand angle

    % Title and axes labels
    title('Average Aftereffect');
    xlabel('Cycle Number') % x-axis label
    ylabel('Hand Angle (deg)') % y-axis label
    set(gca,'FontSize',13) 

%%% MIXED EFFECTS MODEL %%%
    E=D2( (D2.CN >= 106 & D2.CN <= 110) , : ); %select only first few AE cycles
    lme = lmeTrainSched(E,'hand_theta',0);
    
%% Aiming in RNF vs. Hand Angle in AE - Split by GROUP & ROTATION %%
E=D2; %use baseline-subtracted, flipped data
%blockNum 4 = AE, blockNum 6 = RNF 

%%% Average data across subjects for every group %%%
    grpMean = varfun(@nanmean,E,'GroupingVariables',{'Experiment','Group','SN',...
        'blockNum','abs_tgt_rot'},'OutputFormat','table'); %blockNum 4 = AE, blockNum 6 = RNF 

    group = unique(E.Group);
%%% Graph figure %%%
    for gi = 1:length(group)
        sn = unique(E.SN(E.Group == group(gi)));

        subplot(2,2,gi)
        for si = 1:length(sn)
            for ri = 1:length(rot)

                indx = grpMean.Group == group(gi) & grpMean.SN == sn(si) &...
                    grpMean.abs_tgt_rot == rot(ri);

                x = grpMean.nanmean_hand_theta(indx & grpMean.blockNum == 4);
                y = grpMean.nanmean_hand_theta(indx & grpMean.blockNum == 6);
                plot(x,y,'.','markersize',15,'color',color{ri}); hold all;

            end

            % Title and axes labels
            str = sprintf('%s', K_Group_names{gi});
            title(str, 'interpreter','none')
            xlabel('Aftereffect (deg)') % x-axis label
            ylabel('Aiming in RNF (deg)') % y-axis label
            set(gca,'FontSize',13);

            legend('30','45','60')
        end
    end

%% Aiming in RNF vs. Hand Angle in AE - Split by GROUP & TARGET %%
E=D2; %use baseline-subtracted, flipped data
%blockNum 4 = AE, blockNum 6 = RNF 

%%% Average data across subjects for every group %%%
    grpMean = varfun(@nanmean,E,'GroupingVariables',{'Experiment','Group','SN',...
        'blockNum','ti'},'OutputFormat','table'); %blockNum 4 = AE, blockNum 6 = RNF 

    group = unique(E.Group);
%%% Graph figure %%%
    for gi = 1:length(group)
        sn = unique(E.SN(E.Group == group(gi)));

        subplot(2,2,gi)
        for si = 1:length(sn)
            for ti = 1:length(target)

                indx = grpMean.Group == group(gi) & grpMean.SN == sn(si) &...
                    grpMean.ti == target(ti);

                x = grpMean.nanmean_hand_theta(indx & grpMean.blockNum == 4);
                y = grpMean.nanmean_hand_theta(indx & grpMean.blockNum == 6);
                plot(x,y,'.','markersize',15,'color',color{ti}); hold all;

            end

            % Title and axes labels
            str = sprintf('%s', K_Group_names{gi});
            title(str, 'interpreter','none')
            xlabel('Aftereffect (deg)') % x-axis label
            ylabel('Aiming in RNF (deg)') % y-axis label
            set(gca,'FontSize',13);

            legend('30','150','270')
        end
    end
    
    
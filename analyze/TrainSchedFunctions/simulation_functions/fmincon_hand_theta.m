subj = unique(D2.SN);
target = unique(D2.ti);
si = 1;

%%% Create figure
Every_sub = figure(); set(gcf,'units','centimeters','pos',[5 5 30 10]); hold on;
set(gcf,'PaperPositionMode','auto'); set(gca,'FontSize',9);

while si > 0 && si <= length(subj)
    hold on
    
    for ti = 1:length(target)
        
        subplot (1,3,ti)
        
        y_data = D2(D2.SN == subj(si) & D2.ti == target(ti) & D2.CN >= 16 & D2.CN <= 105,:).hand_theta;
        
        rot = -D2(D2.SN == subj(si) & D2.ti == target(ti) & D2.CN >= 16 & D2.CN <= 105,:).abs_tgt_rot;
        ntrials = length(rot);
        
        % sq_err = squared_err(params, trials, rot, y_data)
        
        % Initial parameters, upper & lower bounds
        A0 = 0;
        B0 = 0;
        
        init = [A0, B0];
        LB = [0 0];
        UB = [1 1];
        
        options=optimset('MaxFunEvals',1e14,'TolFun',1e-14,'TolX',1e-14,'display','off');
        [soln, fval] = fmincon(@squared_err,init,[],[],[],[],LB,UB,[],options, ntrials, rot, y_data);
        
        A = soln(1);
        B = soln(2);
        
        y_pred = SSM_simulator(A, B, ntrials, rot);
        
        %     figure; hold on;
        plot(1:ntrials, y_data, '.b'); hold on;
        plot(1:ntrials, y_pred, '.r')
        
         %%% Title and Axis labels
        str = sprintf('SN %d, %uº, A = %.2f, B = %.2f', si, target(ti), A, B);
        title(str,'interpreter','none');
        set(gca,'FontSize',13);
        
    end
    
    %%% Scroll through subjects
    w = waitforbuttonpress;
    if w == 0
        si = si + 1;
    elseif w==1
        si = si - 1;
    end
    clf(Every_sub)
end


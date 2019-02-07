subj = unique(D2.SN(D2.Group > 2));
target = unique(D2.ti);

A = nan(length(subj)*length(target),1);
B = nan(length(subj)*length(target),1);
Group = nan(length(subj)*length(target),1);
SN = nan(length(subj)*length(target),1);
tgt = nan(length(subj)*length(target),1);
rotation = nan(length(subj)*length(target),1);

ind = 0;
for si = 1:length(subj)
    
    for ti = 1:length(target)
        
        y_data = D2(D2.SN == subj(si) & D2.ti == target(ti) & D2.CN >= 16 & D2.CN <= 105,:).hand_theta;
        
        rot = -D2(D2.SN == subj(si) & D2.ti == target(ti) & D2.CN >= 16 & D2.CN <= 105,:).abs_tgt_rot;
        ntrials = length(rot);
        
        % sq_err = squared_err(params, trials, rot, y_data)
        
        % Initial parameters, upper & lower bounds
        A0 = 0;
        B0 = 0;
        
        init = [A0, B0];
        LB = [0.5 0];
        UB = [1 10];
        
        options=optimset('MaxFunEvals',1e14,'TolFun',1e-14,'TolX',1e-14,'display','off');
        [soln, fval] = fmincon(@squared_err,init,[],[],[],[],LB,UB,[],options, ntrials, rot, y_data);
        
        ind = ind + 1; % index
        A(ind) = soln(1);
        B(ind) = soln(2);
        Group(ind) = unique(D2.Group(D2.SN == subj(si)));
        SN(ind) = subj(si);
        tgt(ind) = target(ti);
        rotation(ind) = abs(unique(rot));
        
        
    end
    
end

errorClampFits = table(Group,SN,tgt,rotation,A,B);

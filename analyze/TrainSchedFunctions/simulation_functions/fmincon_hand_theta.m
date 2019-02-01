y_data = D2(D2.SN == 1 & D2.ti == 30 & D2.CN >= 16 & D2.CN <= 105,:).hand_theta;

rot = -D2(D2.SN == 1 & D2.ti == 30 & D2.CN >= 16 & D2.CN <= 105,:).abs_tgt_rot;
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

A = soln(1)
B = soln(2)


% 
% y_pred = SSM_simulator(A, B, ntrials, rot)
% 
% figure; hold on;
% plot(1:ntrials, y_data, '.b')
% plot(1:ntrials, y_pred, '.r')




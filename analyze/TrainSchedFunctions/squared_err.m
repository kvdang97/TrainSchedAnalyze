function sq_err = squared_err(params, ntrials, rot, y_data)

A = params(1);
B = params(2);

y_pred = SSM_simulator(A, B, ntrials, rot);

% Calculate Squared Error - between data and simulation
sq_err = nansum( (y_data' - y_pred).^2);

end
% function sq_err = squared_err(params, trials, y_data , rot)


% A = params(1);
% B = params(2);

A = 0.95;
B = 0.2;

rot = -[repmat(0,10,1) ; repmat(5,40,1) ; repmat(30,50,1)];

e = nan(1,99);
z = nan(1,99);

% starting value
z(1) = 0;

for n = trials
    
    e(n) = z(n) + rot(n);
    
    z(n+1) = A*z(n) - B*(sign(e(n)) ) ;
    
end


figure; 
plot((1:100), z, '.')


y_pred = z;

% Calculate Squared Error - between data and simulation
sq_err = nansum( (y_data - y_pred).^2);


% end
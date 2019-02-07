function z = SSM_simulator(A, B, ntrials, rot)

e = nan(1,ntrials);
z = nan(1,ntrials);

% starting value
z(1) = 0;

for n = 1:ntrials-1
    
    e(n) = z(n) + rot(n); % for VMR
    z(n+1) = A*z(n) - B*e(n) ;
    
    
    % e(n) = rot(n); % for error clamp
    % z(n+1) = A*z(n) - B*(sign(e(n)) ) ;
    
end

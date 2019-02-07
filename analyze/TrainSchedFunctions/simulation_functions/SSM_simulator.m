function z = SSM_simulator(A, B, ntrials, rot)

e = nan(1,ntrials);
z = nan(1,ntrials);

% starting value
z(1) = 0;

for n = 1:ntrials-1    
% % for VMR
% e(n) = z(n) + rot(n);
% z(n+1) = A*z(n) - B*e(n) ;

% for error clamp
e(n) = rot(n);
z(n+1) = A*z(n) - B*(sign(e(n)) ) ;
end

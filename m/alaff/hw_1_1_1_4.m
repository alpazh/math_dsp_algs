format long
		
rng( 0 );      
n = 100;          
U = triu( rand( n,n ) );
x = rand( n,1 );

b = U * x;

xhat = U \ b;

% norm of residual is small, and xhat seems to be good estimate
norm_res = norm( b - U * xhat )
cond_U = cond(U)
% but cond(U) is large and norm of estimation error is very large
norm_err = norm( xhat - x )

% It means that residual norm is a bad indicator 
% of the solution quality
min( abs( diag( U ) ) )

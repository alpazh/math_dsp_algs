format long
		
rng( 0 );      
n = 3          
U = triu( rand( n,n ) )
x = rand( n,1 )

b = U * x
xhat = U \ b
% xhat is different from x and there is non-zero error e_x
e_x = xhat - x
% but residual is still zero
residual_value = U * xhat - b
% e_b = U * xhat - b
% e_b = U * xhat - U * x
% e_b = U * (xhat - x)
e_b2 = U * (xhat - x)

cond_U = cond(U)

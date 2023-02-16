% function integrate_on_S2

ggrav = 9.81;
d =  0.0;
k = 10.0;

tol = 1.0e-10;

t0 = 0; te = 10.0;

% q0 = [ 1/2, 1/2*sqrt(3), 0 ]';
rng ( 'shuffle' )
q0rand = -1 + 2 * rand ( 3, 1 );
q0 = q0rand / norm ( q0rand );

options = odeset ( 'AbsTol', tol, 'RelTol', tol, 'Events', @myEventsFcn );
sol = ode45 ( @(t,y)S2rhs(y,ggrav,d,k), [ t0 te ], [ q0; zeros(3,1) ], options );

nt = 100001;
tt = linspace ( sol.x(1), sol.x(end), nt );

qq   = deval ( sol, tt,    1:3  )';
omom = deval ( sol, tt, 3+(1:3) )';

figure ( 1 )

subplot ( 1, 2, 1 )
plot ( tt, qq )
xlabel ( 't' ),  ylabel ( 'q' ), set ( gca, 'XLim', [ t0 te ] )

subplot ( 1, 2, 2 )
plot ( tt, omom )
xlabel ( 't' ),  ylabel ( '\omega' ), set ( gca, 'XLim', [ t0 te ] )

tstdat = zeros ( nt, 3 );
for it=1:nt
    qqcur = qq(it,:)';
    omcur = omom(it,:)';
    qdcur = cross(omcur,qqcur);
    tstdat(it,1) = norm ( qdcur )^2 / 2 + ggrav * qqcur(3);
    tstdat(it,2) = norm ( qqcur ) - 1;
    tstdat(it,3) = omcur' * qqcur;
end

figure ( 2 )

subplot ( 1, 2, 1 )
plot ( tt, tstdat(:,2) )
xlabel ( 't' ),  ylabel ( '||q||_2 - 1' ), set ( gca, 'XLim', [ t0 te ] )

subplot ( 1, 2, 2 )
semilogy ( tt, abs(tstdat(:,3)) )
xlabel ( 't' ),  ylabel ( '|  \omega^T q  |' ), set ( gca, 'XLim', [ t0 te ] )

figure ( 3 )
plot ( tt, tstdat(:,1) )
xlabel ( 't' ),  ylabel ( 'energy' ), set ( gca, 'XLim', [ t0 te ] )

% end

%==============================================================================%

function yprime = S2rhs ( y, ggrav, d, k )

q = y(1:3);
om = y(3+(1:3));

yprime = zeros ( size ( y ) );

yprime(   1:3 ) = cross ( om, q );
yprime(3+(1:3)) = ggrav * cross ( [0,0,1]', q ) - d * om + k * (om'*q) * q;

end

%==============================================================================%

function [ value, isterminal, direction ] = myEventsFcn ( ~, y )

normq  = norm ( y(   1:3 ) );
normom = norm ( y(3+(1:3)) );

value = zeros ( 3, 1 );
value(1) = normq -  0.1;
value(2) = normq - 10.0;
value(3) = normom - 1.0e+6;

isterminal = ones ( size ( value ) );
direction = zeros ( size ( value ) );

end
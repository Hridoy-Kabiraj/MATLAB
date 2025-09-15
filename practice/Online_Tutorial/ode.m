syms y(t);

ode1 = diff(y, t ) ==t*y;

ySol(t ) = dsolve(ode1)


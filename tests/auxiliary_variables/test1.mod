var x, y;
varexo e;

model;
expectation(-2)(x(+4)*y(+1)) = 1;
y = exp(e);
end;

initval;
y = 1;
x = 1;
end;

steady;

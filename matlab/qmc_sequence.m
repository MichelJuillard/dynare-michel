%@info:
%! @deftypefn {Mex File} {[@var{a}, @var{s}, @var{info}] =} qmc_sequence (@var{d}, @var{s}, @var{t}, @var{n}, @var{lu})
%! @anchor{qmc_sequence}
%! @sp 1
%! Computes quasi Monte-Carlo sequence (Sobol numbers).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item d
%! Integer scalar, dimension.
%! @item s
%! Integer scalar (int64), seed.
%! @item t
%! Integer scalar, sequence type:
%!  @sp 1
%!  @table @ @samp
%!  @item @var{t}=0
%!  Uniform numbers in a n-dimensional (unit by default) hypercube 
%!  @item @var{t}=1
%!  Gaussian numbers
%!  @item @var{t}=2
%!  Uniform numbers on a n-dimensional (unit by default) hypersphere
%!  @end table
%! @item n
%! Integer scalar, number of elements in the sequence.
%! @item lu
%! Optional argument, the interpretation depends on its size:
%!  @sp 1
%!  @table @ @samp
%!  @item @var{d}x2 array of doubles
%!  Lower and upper bounds of the hypercube (default is 0-1 in all dimensions). @var{t} must be equal to zero.
%!  @item @var{d}x@var{d} array of doubles
%!  Covariance matrix of the Gaussian variates (default is the identity matrix). @var{t} must be equal to one.
%!  @item scalar double
%!  Radius of the hypershere (default is one). @var{t} must be equal to two.
%!  @end table
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item a
%! @var{n}x@var{d} matrix of doubles, the Sobol sequence.
%! @item s
%! Integer scalar (int64), current value of the seed.
%! @item info
%! Integer scalar, equal to 1 if mex routine fails, 0 otherwise.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

%@test:1
%$ t = ones(3,1);
%$
%$ d = 2;
%$ n = 100;
%$ s = int64(0);
%$
%$ try
%$   [draws, S] = qmc_sequence(d,s,0,n);
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ try
%$   [draws, S] = qmc_sequence(d,s,1,n);
%$ catch
%$   t(2) = 0;
%$ end
%$
%$ try
%$   [draws, S] = qmc_sequence(d,s,2,n);
%$ catch
%$   t(3) = 0;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = ones(3,1);
%$
%$ d = 2;
%$ n = 100;
%$ s = int64(0);
%$
%$ [draws1, S] = qmc_sequence(d,s,0,n);
%$ [draws2, Q] = qmc_sequence(d,S,0,n);
%$ [draws3, P] = qmc_sequence(d,s,0,2*n);
%$
%$ t(1) = dyn_assert(s,int64(0));
%$ t(2) = dyn_assert(P,Q);
%$ t(3) = dyn_assert([draws1,draws2],draws3);
%$ T = all(t);
%@eof:2

%@test:3
%$ t = ones(3,1);
%$
%$ d = 2;
%$ n = 100;
%$ s = int64(0);
%$
%$ [draws1, S] = qmc_sequence(d,s,0,n,[0 , 2; -1, 2]);
%$ [draws2, Q] = qmc_sequence(d,s,0,n);
%$ 
%$ draws3 = draws2;
%$ draws3(1,:) = 2*draws2(1,:);
%$ draws3(2,:) = 3*draws2(2,:)-1;
%$ t(1) = dyn_assert(S,Q);
%$ t(2) = dyn_assert(draws1,draws3);
%$ T = all(t);
%@eof:3

%@test:4
%$ t = ones(3,1);
%$
%$ d = 2;
%$ n = 100;
%$ s = int64(0);
%$ radius = pi;
%$
%$ [draws, S] = qmc_sequence(d,s,2,n,radius);
%$ 
%$ t(1) = dyn_assert(sqrt(draws(:,3)'*draws(:,3)),radius,1e-14);
%$ t(2) = dyn_assert(sqrt(draws(:,5)'*draws(:,5)),radius,1e-14);
%$ t(3) = dyn_assert(sqrt(draws(:,7)'*draws(:,7)),radius,1e-14);
%$ t(4) = dyn_assert(sqrt(draws(:,11)'*draws(:,11)),radius,1e-14);
%$ t(5) = dyn_assert(sqrt(draws(:,13)'*draws(:,13)),radius,1e-14);
%$ t(6) = dyn_assert(sqrt(draws(:,17)'*draws(:,17)),radius,1e-14);
%$ t(7) = dyn_assert(sqrt(draws(:,19)'*draws(:,19)),radius,1e-14);
%$ t(8) = dyn_assert(sqrt(draws(:,23)'*draws(:,23)),radius,1e-14);
%$ t(9) = dyn_assert(sqrt(draws(:,29)'*draws(:,29)),radius,1e-14);
%$ T = all(t);
%@eof:4

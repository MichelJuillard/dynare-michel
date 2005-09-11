function CreateBenchmark

% stephane.adjemian@cepremap.cnrs.fr [12-06-2004]

global oo_

eval([M_.fname '_oo_ = oo_;'])
eval(['save ' M_.fname '_benchmark_oo ' M_.fname '_oo_'])
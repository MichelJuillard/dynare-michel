function CreateBenchmark

% stephane.adjemian@cepremap.cnrs.fr [12-06-2004]

global fname_ oo_

eval([fname_ '_oo_ = oo_;'])
eval(['save ' fname_ '_benchmark_oo ' fname_ '_oo_'])
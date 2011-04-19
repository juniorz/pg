%Adds benchmakrs to Path
addpath('./benchmarks');

%Globals
global Verbosity;
global NumTrials NumDimensions MaxIterations MaxNoise ETA_p;

MaxIterations   = 1000;
NumDimensions   = 30;
NumTrials       = 30;
MaxNoise        = 1.0;
ETA_p           = -1;  % default eta percentage (see hs_jump.m)

Verbosity = 2;

benchmark = 'all';
jump = 'cj';
noise_list = 0:0.5:1;

%Constants
log_files_path = '~/academico/harmony_search_pg/graphs';
benchmarks = { 'sphere', 'schaffer', 'ackley', 'rosenbrock', 'rastrigin', 'griewank', 'schwefel' };

if strcmp(benchmark, 'all')
  for i = 1:length(benchmarks)
    plot_noise_convergence( 'hs_jump', 'IHS', char(benchmarks(i)), jump, noise_list, log_files_path )
  end
else
  plot_noise_convergence( 'hs_jump', 'IHS', benchmark, jump, noise_list, log_files_path )
end

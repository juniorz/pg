log_files_path = '/Users/juniorz/academico/pg';
iteration_size = {
  1500
%  1600
%  1700
%  1800
%  1900
%  2000
%  3000
%  4000
%  5000
%  6000
%  7000
%  8000
%  9000
%  10000
};

for i = 1:length(iteration_size)
  num_iterations = iteration_size{i,1};
  filename = sprintf('%s/%d.csv', log_files_path, num_iterations);
  
  disp(filename);
  gera_tabela_benchmark(filename, 'all', num_iterations);
end
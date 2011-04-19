log_files_path = '~/academico/harmony_search_pg/benchmarks';

methods_list = { 'IHS' };
jumps_list = { 'cj' };

iteration_size = {
  1000
};



for i = 1:length(iteration_size)
  num_iterations = iteration_size{i,1};
  filename = sprintf('%s/%d.csv', log_files_path, num_iterations);
  
  disp(filename);
  gera_tabela_benchmark(filename, 'all', num_iterations, methods_list, jumps_list);
end
% Gera as tabelas com os resultados do HS
function gera_tabelas(filename, benchmark, verbose)
   
   dimensions = 30;
   trials = 30;
   max_noise = 0.0;

   fid = fopen(filename,'wb');
   
   fprintf(fid,'Function;Method;Noise\n');
   fprintf(fid,';');

   for noise=0:0.2:max_noise
      fprintf(fid,';%g',noise);
   end
   fprintf(fid,'\n');

   if strcmp(benchmark,'all')
      tabela_stats(fid,'sphere',     dimensions, max_noise, verbose, trials);
      tabela_stats(fid,'schaffer',   dimensions, max_noise, verbose, trials);
      tabela_stats(fid,'ackley',     dimensions, max_noise, verbose, trials);
      tabela_stats(fid,'rosenbrock', dimensions, max_noise, verbose, trials);
      tabela_stats(fid,'rastrigin',  dimensions, max_noise, verbose, trials);
      tabela_stats(fid,'griewank',   dimensions, max_noise, verbose, trials);
      % falta a "generalized penalized function"
      tabela_stats(fid,'schwefel',   dimensions, max_noise, verbose, trials);
   else
      tabela_stats(fid, benchmark, dimensions, max_noise, verbose, trials);
   end


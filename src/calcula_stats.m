% Calcula as Estatisticas
% media, desvio padrao, percentual de sucessos
function [media, desvio, success] = calcula_stats(bench, dim, noise, method, flag_jump, eta_p, verbose, tries)

   success = 0;
   results = zeros(tries, 1);

   for i=1:tries
      [out, rate] = feval('hs_jump', bench, dim, noise, method, flag_jump, eta_p, verbose);
      results(i) = out;
      success = success + rate;

      fprintf('run (%d) = %g [%s %s jump=%s noise=%g]\n', i, out, bench, method, flag_jump, noise);
   end
  fprintf('%s\n', results);

   media = mean(results);
   desvio = std(results);
   
   success = success / tries;
end


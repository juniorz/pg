% Calcula as Estatisticas
% media, desvio padrao, percentual de sucessos
function [media, desvio, success] = calcula_stats(bench,dim,noise,method,flag_jump,eta_p,tries)

   success = 0;
   results = zeros(tries, 1);

   for i=1:tries
      [out, rate] = feval('hs_jump', bench, dim, noise, method, flag_jump, eta_p, 0);
      results(i) = out;
      success = success + rate;

      fprintf('run = %d - ans = %g [%s %s %s noise=%g]\n', i, out, bench, method, flag_jump, noise);
   end

   media = mean(results);
   desvio = std(results);
   
   success = success / tries;
end


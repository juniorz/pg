% Calcula as Estatisticas
% media, desvio padrao, percentual de sucessos
function [media, desvio, success] = calcula_stats(bench, method, dim, iteractions, noise, flag_jump, eta_p, tries, verbose)

  success = 0;
  results = zeros(tries, 1);

  for i=1:tries
    [out, rate] = feval('hs_jump', bench, method, dim, iteractions, noise, flag_jump, eta_p, verbose);
    results(i) = out;
    success = success + rate;

    if verbose > 0
      fprintf('run (%d) = %g [%s %s jump=%s noise=%g]\n', i, out, bench, method, flag_jump, noise);
    end
  end

  if verbose > 1
    fprintf('%s\n', results);
  end

  media = mean(results);
  desvio = std(results);

  success = success / tries;
end


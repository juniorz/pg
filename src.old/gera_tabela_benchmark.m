function gera_tabela_benchmark( filename, benchmark, max_iterations )
%GERA_TABELA_BENCHMARK Gera a tabela de um benchmark
%   Para facilitar, o numero de iterações é definido pelo usuario

  % Coloquei fixo
  list_methods = { 'HS' };

  dimensions = 30;
  trials = 30;
  max_noise = 1.0;
  eta_p = -1; % default eta percentage (see chaoticpso.m)
  
  verbose = 1;
  
  fid = fopen(filename, 'wb');
  gera_cabecalho_tabela();
  corpo_tabela();
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function gera_cabecalho_tabela()
    fprintf(fid,'Function;Method;Noise\n');
    fprintf(fid,';');

    for noise=0:0.2:max_noise
      fprintf(fid, ';%g;', noise);
    end
    fprintf(fid,'\n');
  end



  function corpo_tabela()
    if strcmp(benchmark,'all')
      tabela_benchmark('sphere');
      tabela_benchmark('schaffer');
      tabela_benchmark('ackley');
      tabela_benchmark('rosenbrock');
      tabela_benchmark('rastrigin');
      tabela_benchmark('griewank');
      % falta a "generalized penalized function"
      tabela_benchmark('schwefel');
    else
      tabela_benchmark(benchmark);
    end
  end



  % Gera uma tabela
  function tabela_benchmark(bench)
    fprintf(fid,'%s',bench);
    linha_tabela_noise(bench, 'no');
    linha_tabela_noise(bench, 'yes');
  end


  function linha_tabela_noise(bench, flag_jump)

    for met=1:length(list_methods)
      method = char(list_methods(met));

      name = strcat(';',method);

      if strcmp(flag_jump, 'yes')
        name = strcat(name,'+jump');
      end

      fprintf(fid,name);

      for noise=0:0.2:max_noise
        [media,desvio, success] = calcula_stats(bench, method, noise, flag_jump);

        fprintf(fid,';%g (%g)',media,desvio);
      end  

      fprintf(fid,'\n');
    end
  end


  % Calcula as Estatisticas
  % media, desvio padrao, percentual de sucessos
  function [media, desvio, success] = calcula_stats(bench, method, noise, flag_jump)

    success = 0;
    results = zeros(trials, 1);

    for i=1:trials
      [out, rate] = feval('hs_jump', bench, method, dimensions, max_iterations, noise, flag_jump, eta_p, verbose);
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

    success = success / trials;
  end


end


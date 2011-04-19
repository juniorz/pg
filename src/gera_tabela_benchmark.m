function gera_tabela_benchmark( filename, benchmark, max_iterations, methods_list, jumps_list )
%GERA_TABELA_BENCHMARK Gera a tabela de um benchmark
%   Para facilitar, o numero de iterações é definido pelo usuario

  dimensions  = 30;
  trials      = 30;
  max_noise   = 1.0;
  eta_p       = -1;  % default eta percentage (see chaoticpso.m)

  verbose = 1;
  
  fid = fopen(filename, 'wb');
  gera_cabecalho_tabela();
  corpo_tabela();


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function gera_cabecalho_tabela()
    fprintf(fid, 'Function;Method;Noise\n');
    fprintf(fid, ';');

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

    for jmp=1:length(jumps_list)
      jump = char(jumps_list(jmp));
      linha_tabela_noise(bench, jump);
    end
  end


  function linha_tabela_noise(bench, flag_jump)

    for met=1:length(methods_list)
      method = char(methods_list(met));

      name = strcat(';',method);

      if strcmp(flag_jump, 'no') == false
        name = strcat(name, '+', flag_jump);
      end

      fprintf(fid,name);

      for noise=0:0.2:max_noise
        [media, desvio, success] = calcula_stats(bench, method, noise, flag_jump);

        fprintf(fid, ';%g;%g', media, desvio);
      end  

      fprintf(fid, '\n');
    end
  end


  % Calcula as Estatisticas
  % media, desvio padrao, percentual de sucessos
  function [media, desvio, success, mean_progress] = calcula_stats(bench, method, noise, jump_name)

    success = 0;
    results = zeros(trials, 1);
    progress_history = NaN(trials, max_iterations);

    for i=1:trials
      [out, rate, progress] = feval('hs_jump', bench, method, dimensions, max_iterations, noise, jump_name, eta_p, verbose);
      results(i) = out;
      success = success + rate;
      progress_history(i) = progress;
      
      if verbose > 0
        fprintf('run (%d) = %g [%s %s jump=%s noise=%g]\n', i, out, bench, method, jump_name, noise);
      end
    end

    if verbose > 1
      fprintf('%s\n', results);
    end

    media = mean(results);
    desvio = std(results);
    mean_progress = mean(progress_history);
    
    success = success / trials;
  end


end


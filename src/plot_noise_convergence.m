function plot_noise_convergence( metaheuristic, method, bench, jump, noise_list, log_files_path )
%PLOT_CONVERGENCE Summary of this function goes here
%   Detailed explanation goes here

  global MaxIterations;
  series = NaN(length(noise_list), MaxIterations);
  
  filename = sprintf('%s/%s+%s_%s_%d', log_files_path, method, jump, bench, MaxIterations);
  disp(filename);

  for i = 1:length(noise_list)
    noise = noise_list(i);
    [~, ~, ~, progress] = calculate_statistics(metaheuristic, method, bench, noise, jump);
    series(i,:) = progress;
  end

  t = 1:MaxIterations;
  %t = 1:10;
  h = plot_series(series, t);
  saveas(h, filename, 'epsc');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function h=plot_series(series, t)
    h = figure();

    plot( t, series(1, t), '-' ,...
          t, series(2, t), '--' ,...
          t, series(3, t), '-.' ,...
          'LineWidth',1);

    legend('noise 0', 'noise 0.5', 'noise 1.0');

    xlabel('Function Evaluations', 'fontsize' , 12);
    ylabel('GBest','fontsize',12,...%20
                  'rotation',0,...
                  'fontangle','italic',...
                  'fontname','Helvetica');%'Bitstream Charter');

    set(gca,'fontsize',14);
  end
  
end


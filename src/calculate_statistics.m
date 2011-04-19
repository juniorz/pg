function [mean_result, std_dev, success_rate, mean_progress] = calculate_statistics( metaheuristic, method, benchmark, noise, jump_name )
%CALCULATE_STATISTICS Summary of this function goes here
%   Detailed explanation goes here

    global Verbosity;
    global NumTrials NumDimensions MaxIterations MaxNoise ETA_p;

    success_rate = 0;
    results = zeros(NumTrials, 1);
    progress_history = NaN(NumTrials, MaxIterations);


    for i=1:NumTrials
      [out, rate, progress] = feval(metaheuristic, benchmark, method, NumDimensions, MaxIterations, noise, jump_name, ETA_p, Verbosity);
      results(i) = out;
      success_rate = success_rate + rate;
      
      progress_history(i,:) = progress(1,:);
      
      if Verbosity > 0
        fprintf('run (%d) = %g [%s %s jump=%s noise=%g]\n', i, out, benchmark, method, jump_name, noise);
      end
    end

    mean_result = mean(results);
    std_dev = std(results);
    mean_progress = mean(progress_history);
    
    success_rate = success_rate / NumTrials;
end


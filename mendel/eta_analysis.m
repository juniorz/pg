function eta_analysis(bench, dim, run, range)

   noise = 1; % nivel maximo

   %  'run' default
%   if nargin < 4 run = 20; end

%   if nargin < 5 
%      step = 1;
%      start = 1;
%      stop = max_percentage;
%      range = start:step:stop;
%   end
   
   
   points = [];
   %list_methods = {'pso' 'lbest' 'fips' 'bbpso'};
   list_methods = {'pso'};


   for met=1:length(list_methods)
      method = char(list_methods(met));
      for eta=range
         media = 0;
         for i=1:run
            [out, rate] = feval('chaoticpso',bench,dim,noise,method,'yes',eta,0);
            media = media + out;
            
            fprintf('run = %d - ans = %g [%s %s eta=%g%%]\n',i,out,bench,method,eta);
         end
         media = media / run;
         points(end+1,:) = [eta media];
      end

      filename = strcat(bench,'_');
      filename = strcat(filename,method);
      filename = strcat(filename,datestr(fix(clock),'_yyyy-mm-dd_HH:MM:SS'));
      filename = strcat(filename,'_eta.pts');

      save(filename, 'points', '-ascii');
   end % list of methods
end % end function


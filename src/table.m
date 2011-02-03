function table(fid, bench, D, max_noise, trials, list_methods)

   eta_p = -1; % default eta percentage (see chaoticpso.m)
   noise = max_noise;

   if (nargin < 6) | (~length(list_methods))
      list_methods = {'pso' 'lbest' 'fips' 'bbpso'};
   end

   fprintf(fid,'%s',bench);

   for met=1:length(list_methods)
      method = char(list_methods(met));

      name = strcat(';',method);
      fprintf(fid,name);

      for noise=0:0.2:max_noise
         [media,desvio,success] = stats(bench,D,noise,method,'no',eta_p,trials);

         fprintf(fid,';%g (%g)',media,desvio);
      end
      
      fprintf(fid,'\n');
   end

   for met=1:length(list_methods)
      method = char(list_methods(met));

      name = strcat(';',method);
      name = strcat(name,'+cj');
      fprintf(fid,name);

      for noise=0:0.2:max_noise
         [media,desvio,success] = stats(bench,D,noise,method,'yes',eta_p,trials);

         fprintf(fid,';%g (%g)',media,desvio);
      end  

      fprintf(fid,'\n');
   end

end

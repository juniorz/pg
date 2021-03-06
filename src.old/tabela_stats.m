% Gera uma tabela
function tabela_stats(fid, bench, D, max_noise, trials, max_iteractions, verbose)

  % Coloquei fixo
  list_methods = { 'HS' };

   eta_p = -1; % default eta percentage (see chaoticpso.m)
   noise = max_noise;

   fprintf(fid,'%s',bench);

  %tabela SEM o jump
   for met=1:length(list_methods)
      method = char(list_methods(met));

      name = strcat(';',method);
      fprintf(fid,name);

      for noise=0:0.2:max_noise
         [media,desvio,success] = calcula_stats(bench, method, D, max_iteractions, noise, 'no', eta_p, trials, verbose);

         fprintf(fid,';%g;%g',media,desvio);
      end
      
      fprintf(fid,'\n');
   end

  % tabela COM o jump
  for met=1:length(list_methods)
     method = char(list_methods(met));

     name = strcat(';',method);
     name = strcat(name,'+cj');
     fprintf(fid,name);

     for noise=0:0.2:max_noise
        [media,desvio,success] = calcula_stats(bench, method, D, max_iteractions, noise, 'yes', eta_p, trials, verbose);

        fprintf(fid,';%g;%g',media,desvio);
     end  

     fprintf(fid,'\n');
  end

end

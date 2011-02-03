function [media,desvio,success] = stats(bench,dim,noise,method,flag_jump,eta_p,run)

   success = 0;
   lista = zeros(run,1);

   for i=1:run
      [out, rate] = feval('chaoticpso',bench,dim,noise,method,flag_jump,eta_p,0);
      lista(i) = out;
      success = success + rate;

      fprintf('run = %d - ans = %g [%s %s %s noise=%g]\n',i,out,bench,method,flag_jump,noise);
   end

   media = mean(lista);
   desvio = std(lista);
   
   success = success / run;
end


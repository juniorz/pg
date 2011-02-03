function gen(filename, arg, list_methods)
   
   D = 30;
   trials = 30;
   max_noise = 0.0;
   
   if nargin < 3
      list_methods = {};
   end
   
   
   fid = fopen(filename,'wb');
   
   fprintf(fid,'Function;Method;Noise\n');
   fprintf(fid,';');

   for noise=0:0.2:max_noise
      fprintf(fid,';%g',noise);
   end
   fprintf(fid,'\n');

   if strcmp(arg,'all')
      table(fid,'schwefel',D,max_noise,trials,list_methods);
      table(fid,'rastrigin',D,max_noise,trials,list_methods);
      table(fid,'ackley',D,max_noise,trials,list_methods);
      table(fid,'griewank',D,max_noise,trials,list_methods);
      table(fid,'schaffer',D,max_noise,trials,list_methods);
      table(fid,'rosenbrock',D,max_noise,trials,list_methods);
      table(fid,'sphere',D,max_noise,trials,list_methods);
   else
      table(fid,arg,D,max_noise,trials,list_methods);
   end


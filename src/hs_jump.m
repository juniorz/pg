function [BestFitness, rate, full_history] = hs_jump(functname, method, Dim, max_iterations, noise, jump_name, eta_percentage, verbose)
% IHM with Jump
% Baseado em: M. Mahdavi, M. Fesanghary, E. Damangir, An improved harmony search
% algorithm for solving optimization problems, Applied Mathematics and Computation 188 (2007) (2007) 1567–1579.
%   - functname (string): name of the benchmark

%   - Dim (int): dimension number

%   - noise (real): noise level to corrupt the function landscape

%   - method (string): method name (IHS, EGHS)

%   - jump_name (string):
%         'no': dont use jump

%   - eta_percentage (real):
%         -1: to consider default eta
%         otherwise: contains the size of the jump as percentage of the search space.

%   - verbose (int):
%         0: silent
%         1: informative
%         2: extremely verbose mode

  %EGHS have no jump
  if strcmp(method, 'EGHS') && false == strcmp(jump_name, 'no')
    BestFitness = NaN;
    rate = 0;
    return
  end

  % rand('twister', sum(100*clock));
  s = RandStream('mt19937ar', 'Seed', sum(100*clock));
  RandStream.setDefaultStream(s);
   
  % - eta_p: default eta value given as percentage of the search space.
  % - inferior and superior: limits of the search space (the same for all dimension).
  if strcmp(functname, 'rosenbrock')
    inferior = -50;
    superior =  50;
    eta_p = 1;
  elseif strcmp(functname, 'rastrigin')
    inferior = -5.12;
    superior =  5.12;
    eta_p = 10;
  elseif strcmp(functname, 'ackley')
    inferior = -32;
    superior =  32;
    eta_p = 1.6;%1;
  elseif strcmp(functname, 'griewank')
    inferior = -600;
    superior =  600;
    eta_p = 0.1;
  elseif strcmp(functname, 'schaffer')
    Dim = 2;
    inferior = -100;
    superior =  100;
    eta_p = 1;
  elseif strcmp(functname, 'sphere')
    inferior = -100;
    superior =  100;
    eta_p = 0.5;
  elseif strcmp(functname, 'schwefel')
    inferior = -500;
    superior =  500;
    eta_p = 3.5;%1.8;
  end

  %asimetrically initialization
  if strcmp(functname, 'schwefel')
    VRmin = ones(Dim, 1) * inferior;
    VRmax = ones(Dim, 1) * inferior/2;
  else
    VRmin = ones(Dim, 1) * superior/2;
    VRmax = ones(Dim, 1) * superior;
  end

  % use default eta percentage (-1: yes)
  if eta_percentage == -1
    eta_percentage = eta_p;
  end
  
  % not default
  % eta = percentage of the search space
  if eta_percentage > 0
    eta = (eta_percentage/100.0) * (superior - inferior);
  end

  % inicializar sequencia caotica com valor diferente de: 0,0.25,0.5,0.75,1.0
  cj=rand;
  while (cj==0) || (cj==0.25) || (cj==0.5) || (cj==0.75) || (cj==1.0)
    cj=rand;
  end

  message = 'HS: %g/%g iterations, GBest = %g\n';
  local_search_count = 5; % numero maximo de falhas antes de chavear para salto

  % Variaveis para o HS
  %global NVAR
%  global NG NH MaxItr HMS HMCR PARmin PARmax bw_min bw_max;
%  global HM NCHV fitness_m BW;
%  global BestIndex WorstIndex BestFit WorstFit BestGen currentIteration;

  WorstFit = NaN;
  WorstIndex = NaN;
  BestIndex = NaN;
  BestFit = NaN;
  BestGen = NaN;
  %currentIteration = NaN;

  % HS general parameters
  MaxItr = max_iterations; % maximum number of iterations
  HMS = 5;           % harmony memory size
	% NG=6;            %number of ineguality constraints
  % NH=0;            %number of eguality constraints

  
  % IHS parameters
  HMCR    = 0.9;       % harmony memory consideration rate  0< HMCR <1
  PARmin  = 0.01;      % minumum pitch adjusting rate
  PARmax  = 0.99;      % maximum pitch adjusting rate
  bw_min   = 0.0001;    % minumum bandwidth
  bw_max   = 1 / (20 * (superior - inferior)); % maxiumum bandwidth

  
  % EGHS parameters
  LUP = 0.9;          %Location Updating Probability
  c_max = 1;
  c_min = 0;

  
  % /**** Initiate Matrix ****/
  HM          = zeros(HMS, Dim);
  NCHV        = zeros(1, Dim);    % New Candidate Harmony Vector??
  BestGen     = zeros(1, Dim);
  fitness_m 	= zeros(1, HMS);

  
  %Iteration-based Coeficients
  PAR         = 0;
  BW          = 0;
  Copt        = 0;                % Coefficient of optimization scale
  %BW          = zeros(1, Dim);
  %gx          = zeros(1, NG);



  % randomly initialize the HM
  HM = zeros(HMS, Dim);
  for d=1:Dim
    HM(1:HMS,d) = unifrnd(VRmin(d)*ones(HMS,1), VRmax(d)*ones(HMS,1));
  end

  for i=1:HMS
    fitness_m(i) = fitness(functname, HM(i,1:Dim), noise);
  end


  %Global Best
  [gbestval, idx] = min(fitness_m);
  gbest = HM(idx,:);
  
  %Global Worst
  [gworstval, idx] = max(fitness_m);
  gworst = HM(idx,:);
  
  % START
  pulou         = 0;                  % flag para controle de salto: sim(1) ou nao(0).   
  fail_count    = 0;                  % contador de falhas
  jump_counter  = zeros(MaxItr, 2);   % contador cumulativo: (1) jumps com sucesso e (2) jumps realizados

  % Para grafico de convergencia
  full_history = NaN(1, MaxItr);

  %Implementation
  currentIteration  = 0;
  while(StopCondition(currentIteration))
    pulou=0;
    
    %Pre Operations
    set_iteration_based_coeficients();
    
    %Main Operations
    improvise_new_harmony_vector();
    new_fitness = fitness(functname, NCHV, noise);

    if verbose > 2
      fprintf('new_fitness = %g\n', new_fitness);
      display(NCHV);
    end

    update_harmony_memory( new_fitness );

    %Post Operations
    update_gbest()
    update_gworst()
    update_jump_counter()

    % mostrar mensangem de acompanhamento do processo
    if verbose > 2
        fprintf(message, currentIteration+1, MaxItr, gbestval);
    end
    
    currentIteration = currentIteration+1;
  end
  
  BestFitness = min(fitness_m);

  % taxa de sucesso com jump.
  rate = jump_counter(end,1)*100/jump_counter(end,2);
  
% /*********************************************/

  function set_iteration_based_coeficients()
    if strcmp(method, 'IHM')
      PAR = PARmin + (PARmax - PARmin) / (MaxItr) * currentIteration;

      coef = log(bw_min / bw_max) / MaxItr;
      BW = bw_max * exp(coef * currentIteration);
    end
    
    if strcmp(method, 'EGHS')
      Copt = c_max - (c_max - c_min) / MaxItr * currentIteration;
    end
  end


  function improvise_new_harmony_vector()

    for i = 1:Dim
      
      %Improved Harmony Search
      if strcmp(method, 'IHM')
        if strcmp(jump_name, 'no') || (fail_count <= local_search_count)
          improvise_new_harmony()
        else
          perform_jump()
        end
      end
    
      %Enhanced Global Harmony Search
      if strcmp(method, 'EGHS')
        if rand(1) <= LUP
          NCHV(i) = gbest(i) - Copt * rand(1) * (gbest(i) - gworst(i));
        else
          NCHV(i) = inferior + rand(1) * (superior - inferior);
        end
      end
      
    end
    
  end



  function improvise_new_harmony()

    if( rand(1) < HMCR ) % memory consideration

      index = randint(1, HMS); %get random memory vector
      NCHV(i) = HM(index, i);

      if( rand(1) < PAR) % pitch adjusting

        result = NCHV(i);

        if( rand(1) < 0.5)
            result = result +  rand(1) * BW(i);

            if( result < superior)
                NCHV(i) = result;
            end
        else
            result = result - rand(1) * BW(i);

            if( result > inferior)
                NCHV(i) = result;
            end
        end
      end
    else
      NCHV(i) = randval( inferior, superior ); % random selection
    end
    
  end


  function perform_jump()
    pulou = 1;

    % IHS + Chaotic-Jump
    % Avaliar! Ele faz o salto em relação as particulas. Eu estou
    % fazendo em relação ao melhor global. Ele dá um salto em uma
    % partícula. Eu terei que dar um salto em uma componente. Ele
    % usa a falha por partícula. EU terei que usar falhas globais.
    % OU eu posso substituir o que ele faz por partícula, fazer por
    % coordenada, mas não deve funcionar muito bem.
    if strcmp(jump_flag, 'cj')
      cj = 4*cj*(1-cj);

      % O valor naquela dimensão é obtido pelo valor naquela
      % dimensão para a melhor solução até o momento MAIS um
      % salto
      NCHV(i) = gbest(i) * (1 + eta*(2*cj-1));

      % Verifica se o valor gerado está dentro dos limites da
      % função
      if NCHV(i) > superior
        NCHV(i) = gbest(i);
      end

      if NCHV(i) < inferior
        NCHV(i) = gbest(i);
      end
    end
  end


  function update_harmony_memory( NewFit )

    if(currentIteration == 0)
      BestFit=fitness_m(1);
      for i = 1:HMS
        if( fitness_m(i) < BestFit )
          BestFit = fitness_m(i);
          BestIndex =i;
        end
      end

      WorstFit = fitness_m(1);
      WorstIndex = 1;
      for i = 1:HMS
        if( fitness_m(i) > WorstFit )
          WorstFit = fitness_m(i);
          WorstIndex =i;
        end
      end
    end

    if (NewFit< WorstFit)

      if( NewFit < BestFit )
        %Deu um erro aqui e eu nao sei porque
        %Subscripted assignment dimension mismatch.
        HM(WorstIndex,:) = NCHV;
        BestGen=NCHV;
        fitness_m(WorstIndex)=NewFit;
        BestIndex=WorstIndex;
      else
        HM(WorstIndex,:) = NCHV;
        fitness_m(WorstIndex)=NewFit;
      end

      WorstFit = fitness_m(1);
      WorstIndex = 1;
      for i = 1:HMS
        if( fitness_m(i) > WorstFit )
          WorstFit = fitness_m(i);
          WorstIndex =i;
        end
      end

    end
  end

  function update_gbest()
    
    %Para fazer o grafico de convergencia
    full_history(1, currentIteration+1) = gbestval;    

    if gbestval > new_fitness
      gbestval = new_fitness;
      gbest = NCHV;

      if pulou == 1
        jump_counter(currentIteration+1, 1) = jump_counter(currentIteration+1, 1) + 1;
      end
    else
      fail_count = fail_count + 1;
    end
  end

  function update_gworst()
    if gworstval < new_fitness
      gworstval = new_fitness;
      gworst = NCHV;
    end
  end

  function update_jump_counter()
    % contador: total jumps
    if pulou == 1
      fail_count = 0;
      jump_counter(currentIteration+1, 2) = jump_counter(currentIteration+1, 2) + 1;
    end
    
    % cumulative sucess sum
    jump_counter(currentIteration+2,:) = jump_counter(currentIteration+1,:);

  end


end %function

% /*****************************************/
function val1=randval(Maxv,Minv)
    val1=rand(1)*(Maxv-Minv)+Minv;
end

function val2=randint(Maxv,Minv)
    val2=round(rand(1)*(Maxv-Minv)+Minv);
end
% /*******************************************/

function val=StopCondition(Itr)
    global MaxItr;
    val=1;
    if(Itr>=MaxItr)
        val=0;
    end
end

% /*******************************************/

function f = fitness(functname, pos, noise)
   f = feval(functname, pos);
   f = f + noise*randn;
end

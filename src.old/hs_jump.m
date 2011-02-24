function [out, rate] = hs_jump(functname, method, Dim, max_iterations, noise, flag_jump, eta_percentage, verbose)
% IHM with Jump
% Baseado em: M. Mahdavi, M. Fesanghary, E. Damangir, An improved harmony search
% algorithm for solving optimization problems, Applied Mathematics and Computation 188 (2007) (2007) 1567–1579.
%   - functname (string): name of the benchmark

%   - Dim (int): dimension number

%   - noise (real): noise level to corrupt the function landscape

%   - method (string): nome do metodo a ser utilizado. Atualmente apenas HS

%   - flag_jump (string):
%         'yes': to use the jump strategy
%         'no': otherwise

%   - eta_percentage (real):
%         -1: to consider default eta
%         otherwise: contains the size of the jump as percentage of the search space.

%   - verbose (int):
%         >0: extremely verbose mode

  % rand('twister', sum(100*clock));
  s = RandStream('mt19937ar', 'Seed', sum(100*clock));
  RandStream.setDefaultStream(s);
   
  % - eta_p: default eta value given as percentage of the search space.
  % - inferior and superior: limits of the search space (the same for all dimension).
  if strcmp(functname,'rosenbrock')
    inferior = -50;
    superior =  50;
    eta_p = 1;
  elseif strcmp(functname,'rastrigin')
    inferior = -5.12;
    superior =  5.12;
    eta_p = 10;
  elseif strcmp(functname,'ackley')
    inferior = -32;
    superior =  32;
    eta_p = 1.6;%1;
  elseif strcmp(functname,'griewank')
    inferior = -600;
    superior =  600;
    eta_p = 0.1;
  elseif strcmp(functname,'schaffer')
    Dim = 2;
    inferior = -100;
    superior =  100;
    eta_p = 1;
  elseif strcmp(functname,'sphere')
    inferior = -100;
    superior =  100;
    eta_p = 0.5;
  elseif strcmp(functname,'schwefel')
    inferior = -500;
    superior =  500;
    eta_p = 3.5;%1.8;
  end

  %asimetrically initialization
  if strcmp(functname,'schwefel')
    VRmin=ones(Dim,1)*inferior;
    VRmax=ones(Dim,1)*inferior/2;
  else
    VRmin=ones(Dim,1)*superior/2;
    VRmax=ones(Dim,1)*superior;
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

  message = sprintf('HS: %%g/%%g iterations, GBest = %%g\n');
  local_search_count = 5; % numero maximo de falhas antes de chavear para salto

  % Variaveis para o HS
  %global NVAR
  global NG NH MaxItr HMS HMCR PARmin PARmax bwmin bwmax;
  global HM NCHV fitness_m BW;
  global BestIndex WorstIndex BestFit WorstFit BestGen currentIteration;

  %MaxItr = 1500;     % maximum number of iterations
  MaxItr = max_iterations;
  HMS = 5;           % harmony memory size

	% NG=6;            %number of ineguality constraints
  % NH=0;            %number of eguality constraints

  % Parameters
  HMCR    = 0.9;       % harmony memory consideration rate  0< HMCR <1
  PARmin  = 0.01;      % minumum pitch adjusting rate
  PARmax  = 0.99;      % maximum pitch adjusting rate
  bwmin   = 0.0001;    % minumum bandwidth
  bwmax   = 1 / (20 * (superior - inferior)); % maxiumum bandwidth

  % /**** Initiate Matrix ****/
  HM          = zeros(HMS, Dim);
  NCHV        = zeros(1, Dim);         % New Candidate Harmony Vector??
  BestGen     = zeros(1, Dim);
  fitness_m 	= zeros(1, HMS);
  BW          = zeros(1, Dim);
  % gx          = zeros(1, NG);

  pulou = 0;              % flag para controle de salto: sim(1) ou nao(0).   
  fail_count = 0;         % contador de falhas
  success_counter = zeros(MaxItr, 2); % contador: jumps com sucesso por iteracao
    
  % randomly initialize the HM
  HM = zeros(HMS, Dim);
  for d=1:Dim
    HM(1:HMS,d) = unifrnd(VRmin(d)*ones(HMS,1), VRmax(d)*ones(HMS,1));
  end

  for i=1:HMS
    %fitness(i) = Fitness(HM(i,:));
    fitness_m(i) = fitness(functname, HM(i,1:Dim), noise);
  end
  
  
  [gbestval, idx] = min(fitness_m);
  gbest = HM(idx,:);  % this is gbest position

   % START
  currentIteration  = 0;
  while(StopCondition(currentIteration))

    PAR=(PARmax-PARmin)/(MaxItr)*currentIteration+PARmin;
    coef=log(bwmin/bwmax)/MaxItr;
    for pp =1:Dim
      BW(pp)=bwmax*exp(coef*currentIteration);
    end

    pulou=0;

    % improvise a new harmony vector
    for i = 1:Dim

      % Default IHS
      if strcmp(flag_jump, 'no') || (fail_count <= local_search_count)

        ran = rand(1);
        if( ran < HMCR ) % memory consideration

          index = randint(1, HMS); %get random memory vector
          NCHV(i) = HM(index, i);

          pvbRan = rand(1);
          if( pvbRan < PAR) % pitch adjusting

            pvbRan1 = rand(1);
            result = NCHV(i);

            if( pvbRan1 < 0.5)
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

      % IHS + Chaotic-Jump
      % Avaliar! Ele faz o salto em relação as particulas. Eu estou
      % fazendo em relação ao melhor global. Ele dá um salto em uma
      % partícula. Eu terei que dar um salto em uma componente. Ele
      % usa a falha por partícula. EU terei que usar falhas globais.
      % OU eu posso substituir o que ele faz por partícula, fazer por
      % coordenada, mas não deve funcionar muito bem.
      else
        pulou = 1;

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

    %newFitness = Fitness(NCHV);
    newFitness = fitness(functname, NCHV, noise);
    UpdateHM( newFitness );


    % Avaliar!
    %----------------------------
    % contador: total jumps
    if pulou == 1
      success_counter(currentIteration+1, 2) = success_counter(currentIteration+1, 2) + 1;
    end

    % contador: reseta contador de falha apos um salto.
    if pulou == 1
      fail_count = 0;
    end

    % atualizar pbest se necessario
    if gbestval > newFitness
      gbestval = newFitness;
      gbest = NCHV;

      % contador: successfull jump
      if pulou==1
        success_counter(currentIteration+1, 1) = success_counter(currentIteration+1, 1) + 1;
      end
    else
      % contador: usado para decidir quando saltar
      fail_count = fail_count + 1;
    end


    % mostrar mensangem de acompanhamento do processo
    if verbose > 2
        fprintf(message, currentIteration, MaxItr, gbestval);
    end

    % cumulative sucess sum
    success_counter(currentIteration+2,:) = success_counter(currentIteration+1,:);

    %---------------------------

    currentIteration = currentIteration+1;
  end
  
  BestFitness = min(fitness_m);
  %END - MainHarmony
   
  % function output
  out = BestFitness;

  % taxa de sucesso com jump.
  rate = success_counter(end,1)*100/success_counter(end,2);

% /*********************************************/

    function UpdateHM( NewFit )
        
        if(currentIteration==0)
            BestFit=fitness_m(1);
            for i = 1:HMS
                if( fitness_m(i) < BestFit )
                    BestFit = fitness_m(i);
                    BestIndex =i;
                end
            end
            
            WorstFit=fitness_m(1);
            for i = 1:HMS
                if( fitness_m(i) > WorstFit )
                    WorstFit = fitness_m(i);
                    WorstIndex =i;
                end
            end
        end
        
        if (NewFit< WorstFit)
            
            if( NewFit < BestFit )
                HM(WorstIndex,:)=NCHV;
                BestGen=NCHV;
                fitness_m(WorstIndex)=NewFit;
                BestIndex=WorstIndex;
            else
                HM(WorstIndex,:)=NCHV;
                fitness_m(WorstIndex)=NewFit;
            end

            WorstFit=fitness_m(1);
            WorstIndex =1;
            for i = 1:HMS
                if( fitness_m(i) > WorstFit )
                    WorstFit = fitness_m(i);
                    WorstIndex =i;
                end
            end
            
        end
    end % main if
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
    if(Itr>MaxItr)
        val=0;
    end
end

% /*******************************************/

function f = fitness(functname, pos, noise)
   f = feval(functname, pos);
   f = f + noise*randn;
end

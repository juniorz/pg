function [out, rate] = hs_jump(functname, Dim, noise, flag_jump, eta_percentage, display_process)
   
   rand('twister',sum(100*clock));

%   - functname (string): name of the benchmark

%   - Dim (int): dimension number

%   - noise (real): noise level to corrupt the function landscape

%   - flag_jump (string):
%         'yes': to use the jump strategy
%         'no': otherwise

%   - eta_percentage (real):
%         -1: to consider default eta
%         otherwise: contains the size of the jump as percentage of the search space.

%   - display_process (int):
%         >0: extremely verbose mode


   
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

%   message = sprintf('HS: %%g/%%g iterations, GBest = %%g\n');

    % inicializar sequencia caotica com valor diferente de: 0,0.25,0.5,0.75,1.0
    cj=rand;
    while (cj==0) || (cj==0.25) || (cj==0.5) || (cj==0.75) || (cj==1.0)
        cj=rand;
    end

    % Variaveis para o HS
    %global NVAR
    global NG NH MaxItr HMS HMCR PARmin PARmax bwmin bwmax;
    global HM NCHV fitness_m PVB BW gx;
    global BestIndex WorstIndex BestFit WorstFit BestGen currentIteration;

    MaxItr=5000;     % maximum number of iterations
%    NVAR=D;         %number of variables
%    NG=6;           %number of ineguality constraints
%    NH=0;           %number of eguality constraints
    HMS=6;           % harmony memory size
    HMCR=0.9;        % harmony consideration rate  0< HMCR <1
    PARmin=0.4;      % minumum pitch adjusting rate
    PARmax=0.9;      % maximum pitch adjusting rate
    bwmin=0.0001;    % minumum bandwidth
    bwmax=1.0;       % maxiumum bandwidth

    % O que seria isso? Uma sequencia caótica?
    % Teria que 
    PVB=[1.0 4;0.6 2;40 80;20 60];   % range of variables

    % /**** Initiate Matrix ****/
    HM=zeros(HMS,Dim);
    NCHV=zeros(1,Dim);         % New Candidate Harmony Vector??
    BestGen=zeros(1,Dim);
    fitness_m=zeros(1,HMS);
    BW=zeros(1,Dim);
    gx=zeros(1,NG);
    % warning off MATLAB:m_warning_end_without_block
   

    %MainHarmony
    %initialize
    % randomly initialize the HM
    for i=1:HMS
        for j=1:Dim
            HM(i,j)=randval(PVB(j,1),PVB(j,2));
        end
        %fitness(i) = Fitness(HM(i,:));
        fitness_m(i) = fitness(functname, HM(i,1:Dim), noise);
    end
    %END - initialize
    currentIteration  = 0;

    while(StopCondition(currentIteration))

        PAR=(PARmax-PARmin)/(MaxItr)*currentIteration+PARmin;
        coef=log(bwmin/bwmax)/MaxItr;
        for pp =1:Dim
            BW(pp)=bwmax*exp(coef*currentIteration);
        end
        % improvise a new harmony vector
        for i =1:Dim
            ran = rand(1);
            if( ran < HMCR ) % memory consideration
                index = randint(1,HMS);
                NCHV(i) = HM(index,i);
                pvbRan = rand(1);
                if( pvbRan < PAR) % pitch adjusting
                    pvbRan1 = rand(1);
                    result = NCHV(i);
                    if( pvbRan1 < 0.5)
                        result =result+  rand(1) * BW(i);
                        if( result < PVB(i,2))
                            NCHV(i) = result;
                        end
                    else
                        result =result- rand(1) * BW(i);
                        if( result > PVB(i,1))
                            NCHV(i) = result;
                        end
                    end
                end
            else
                NCHV(i) = randval( PVB(i,1), PVB(i,2) ); % random selection
            end
        end
        %newFitness = Fitness(NCHV);
        newFitness = fitness(functname, NCHV, noise);
        UpdateHM( newFitness );

        currentIteration=currentIteration+1;
    end
    BestFitness = min(fitness_m);
   %END - MainHarmony

   % Teria que adicionar os JUMPs no perform_HS

   % function output
   out = BestFitness;
   
   % salvar o contador de jump apenas quando estiver em modo verbose.
   if display_process > 0
%      save_sucess(success_counter, functname, flag_method);
   end

   % taxa de sucesso com jump.
   %rate = success_counter(end,1)*100/success_counter(end,2);
   rate = 1;

% /*********************************************/

    function UpdateHM( NewFit )
        % global Dim MaxItr HMS ;
        % global HM NCHV BestGen fitness ;
        % global BestIndex WorstIndex BestFit WorstFit currentIteration;
        
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


function save_sucess(success_counter, bench, method)

   sucess_rate = success_counter(:,1)*100 ./ success_counter(:,2);
   
   sucess_rate(isnan(sucess_rate)) = 100;

   filename = strcat(bench,'_');
   filename = strcat(filename,method);
   filename = strcat(filename,datestr(fix(clock),'_yyyy-mm-dd_HH:MM:SS'));
   filename = strcat(filename,'_success.pts');

   save(filename,'sucess_rate','-ascii');
   
end
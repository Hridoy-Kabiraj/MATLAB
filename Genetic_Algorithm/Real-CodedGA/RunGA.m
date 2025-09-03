function out = RunGA(problem, params)

    % Problem
    CostFunction = problem.CostFunction;
    nVar = problem.nVar;
    VarMin = problem.VarMin;
    VarMax = problem.VarMax;
    VarSize = [1, nVar];

    % Params
    MaxIt = params.MaxIt;
    nPop = params.nPop;
    pC = params.pC;
    nC = round(pC*nPop/2)*2;   % We need even number for offsprings
    gamma = params.gamma;
    mu = params.mu;
    beta = params.beta;
    sigma = params.sigma;

    % Template for empty individuals
    empty_individuals.Position = [];
    empty_individuals.Cost = [];

    % Best solution ever found
    bestsol.Cost = inf;

    %Initialization
    pop = repmat(empty_individuals, nPop, 1);
    for i = 1:nPop

        % Generate Random Solution
        pop(i).Position = unifrnd(VarMin, VarMax, VarSize);

        % Evaluate Solution
        pop(i).Cost = CostFunction(pop(i).Position);
    
        % Compare solution to best solution ever found
        if pop(i).Cost < bestsol.Cost
            bestsol = pop(i);
        end
    end

    % Best Cost of Iteration
    bestcost = nan(MaxIt, 1);

    % Main Loop
    for it = 1:MaxIt

        % Initialize offsprings population
        popc = repmat(empty_individuals, nC/2, 2);

        % Crossover
        for k = 1:nC/2

            % Selection probabilities
            c = [pop.Cost];
            avgc = mean(c);
            if avgc ~= 0
                c = c/avgc;
            end
            probs = exp(-beta*c);

            % Select parents
            p1 = pop(RouletteWheelSelection(probs));
            p2 = pop(RouletteWheelSelection(probs));

            % Perfrom corssover
            [popc(k, 1).Position, popc(k, 2).Position] = UniformCrossover(p1.Position, p2.Position, gamma);

        end

        % Convert popc to Single column matrix
        popc = popc(:);

        % Mutation
        for l = 1:nC

            % Perform mutation
            popc(l).Position = Mutate(popc(l).Position, mu, sigma);

            %Check for variable bounds
            popc(l).Position = max(popc(l).Position, VarMin);
            popc(l).Position = min(popc(l).Position, VarMax);

            % Evaluate
            popc(l).Cost = CostFunction(popc(l).Position);

            % Compare solution to best solution ever found
            if popc(l).Cost < bestsol.Cost
                bestsol = popc(l);
            end

        end

       
        % Merge population
        pop = [pop; popc];

        % Sort population
        pop = SortPopulation(pop);

        % remove extra individuals
        pop = pop(1:nPop);

        % update best Cost of Iteration
        bestcost(it) = bestsol.Cost;

        % Display Iteration Information
        disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(bestcost(it))]);

    end

    % Results
    out.pop = pop;
    out.bestsol = bestsol;
    out.bestcost = bestcost;

end
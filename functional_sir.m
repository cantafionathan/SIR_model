% ADAPTED FROM MODEL IN THIS PAPER
% Turkyilmazoglu M. An extended epidemic model with vaccination: Weak-immune SIRVI. 
% Physica A. 2022;598:127429. doi:10.1016/j.physa.2022.127429

rng(404);


%-----MODEL-----%

% 3^4 factorial. 
isolation = [0, 0.5, 1]; % describe what "Isolation" means in the model - low, medium, high
vac_rate = [0, 0.01, 0.02]; % based on historical data of percentage of population vaccinated per day, 
% the absolute most was china who managed 1.57% per day, with current
% knowledge we think countries could certainly get this up to 2% per day if needed
quar_dur = [0, 7, 14]; % based on the fact that 14 day quarantine was recommended 
num_daily = [15, 30, 45]; % this seems like a reasonable average number for low, middle, and high density populations

[A, B, C, D] = ndgrid(isolation, vac_rate, quar_dur, num_daily);
Av = A(:);
Bv = B(:);
Cv = C(:);
Dv = D(:);

for i = 1:length(Av)
    % randomly generated levels
    prob_spread = 0.1*rand(1)+0.05; % random number between (0.05, 0.15)
    recovery_rate = 1/(14*rand(1) + 7); % random number between (1/7, 1/21)
    vac_eff = 0.004*rand(1) + 0.001; % random number between (0.001, 0.005) ~ (99.5, 99.9)
    prob_sympt = 0.6*rand(1); % random number between (0, 0.6)

    load(i) = pandemic(Dv(i), Av(i), prob_spread, recovery_rate, Bv(i), vac_eff, prob_sympt, Cv(i));
end

M = [load', Av, Bv, Cv, Dv];
Labels = ["load", "soc.iso", "rate.vac", "quar.dur", "num.daily"];
writetable(array2table(M, 'VariableNames', Labels), 'pandemic_data.csv')


function load = pandemic(num_daily, isolation, prob_spread, recovery_rate, vac_rate, vac_eff, prob_sympt, quar_dur)
    
    % SIVQR model
    N = 1000000; % total population
    I0 = 0.01*N; % number of people sick once the health authority notices an issue
    S = N - I0; % total number of susceptible individuals
    I = I0; % total number of infected individuals
    V = 0; % total number of vaccinated individuals
    Q = 0; % total number of quarantined individuals
    R = 0; % total number of recovered individuals
    
    %-----Parameters-----%
    
    % number of daily interactions if there was no policies
    % proxy for population density
    base_num_inter = num_daily; 
    
    % social isolation policy ranging from zero isolation to total isolation
    % reduces number of daily interactions
    isolation_factor = 0.9*isolation; % note the factor of 0.9
    
    % average number of daily interactions
    daily_interactions = base_num_inter.*(1-isolation_factor); 
    
    % nu = vaccination rate of susceptible individuals
    nu = vac_rate; % levels 0, 0.1, 0.2, 0.3
    
    % vaccine effectivness
    vac_eff = vac_eff;
    
    % probability of getting infected from interacting w/ infected invdividual
    probability_of_spread_unvaccinated = prob_spread; 
    probability_of_spread_vaccinated = vac_eff.*probability_of_spread_unvaccinated; 
    
    % beta = transmission rate
    beta_unvac = daily_interactions.*probability_of_spread_unvaccinated;
    beta_vac = daily_interactions.*probability_of_spread_vaccinated;
    
    % gamma = recovery rate
    gamma = recovery_rate; % covid recovery in 2 weeks
    
    %probability that an infectious individual is sympomatic
    prob_symptomatic = prob_sympt; % covid
    
    % quarantining sympomatic individuals
    % should always be set to 1, if no quarantine set quarantine duration to 0
    quarantine = 1; 
    
    qd = quar_dur;
    quarantine_duration = min(qd, 1./gamma); % qd days or recovery
    % rate at which we take people out of quarantine
    if quarantine_duration == 0
        quarantine_rate = 0;
        quarantine = 0;
    else
        quarantine_rate = 1/quarantine_duration;
    end
    
    %rate of recovery for quarantined individuals:
    if quarantine_rate == 0
        gamma_q = gamma;
    elseif quarantine_duration >= 1./gamma
        gamma_q = 1;
    else
        gamma_q = 1./(1./gamma - 1./quarantine_rate);
    end
    
    %-----Solver-----%
    
    % ds/dt = -beta_unvac*s*i - nu*s
    % di/dt beta_unvac*s*i + beta_vac*v*i - quarantine*prob_sympomatic*i + (1-gamma_q)*quarantine_rate*q - gamma*i
    % dv/dt = nu*s- beta_vac*v*i
    % dq/t = quarantine*prob_symptomatic*i - quarantine_rate*q
    % dr/dt = gamma*i + gamma_q*quarantine_rate*q \\ can be solved for from s + i + v + r = 1
    
    % diff equations with forward euler:
    k = 0.001; % lower time step for convergence, set to k = 1 for interpretability
    Tf = 90; % days
    numsteps = ceil(Tf/k);
    k = Tf/numsteps;
    
    s = zeros(1, numsteps);
    i = zeros(1, numsteps);
    v = zeros(1, numsteps);
    q = zeros(1, numsteps);
    r = zeros(1, numsteps);
    temp = zeros(1, numsteps);
    
    s(1) = S/N; % prop of population susceptible
    i(1) = I/N; % prop of population infected
    v(1) = V/N; % prop of vaccinated individuals
    q(1) = Q/N; % prop of quarantined individuals
    r(1) = R/N; % prop of population recovered
    
    for n = 1:(numsteps-1)
        s(n+1) = s(n) + k.*(-beta_unvac.*s(n).*i(n) - nu.*s(n));
        i(n+1) = i(n) + k.*(beta_unvac.*s(n).*i(n) + beta_vac.*v(n).*i(n) - quarantine.*prob_symptomatic.*i(n) + (1-gamma_q).*quarantine_rate.*q(n) - gamma.*i(n));
        v(n+1) = v(n) + k.*(nu.*s(n)- beta_vac.*v(n).*i(n));
        q(n+1) = q(n) + k.*(quarantine.*prob_symptomatic.*i(n) - quarantine_rate.*q(n));
        r(n+1) = r(n) + k.*(gamma.*i(n) + gamma_q.*quarantine_rate.*q(n));
    end
    
    % plot(1:numsteps, s, 'b', ...
    %     1:numsteps, i, 'r', ...
    %     1:numsteps, v, 'g', ...
    %     1:numsteps, q, 'm', ...
    %     1:numsteps, r, 'k')
    % xticks(0:10/k:numsteps)
    % xticklabels(0:10:Tf)
    % legend('Susceptible', 'Infected', 'Vaccinated', 'Quarantined', 'Recovered')
    
    % fprintf('PARAMETERS/TREATMENTS: \n')
    % fprintf('average number of daily interactions (population density) = %d \n', base_num_inter);
    % fprintf('social isolation ranging from 0 (no isolation), to 1 (max isolation) = %d \n',isolation_factor);
    % fprintf('probability of spreading to unvaccinated = %d \n',probability_of_spread_unvaccinated);
    % fprintf('rate of recovery = %d days \n', 1/gamma);
    % fprintf('rate at which susceptible individuals are vaccinated = %d \n',nu);
    % fprintf('vaccine effectivness (how much does it reduce prob(spread)?) = %d \n',vac_eff);
    % fprintf('prob_symptomatic = %d \n',prob_symptomatic);
    % fprintf('quarantine_duration = %d days \n',quarantine_duration);
    % fprintf('RESPONSE: \n')
    % fprintf('peak of infected population = %d%% \n', round(max(i)*100));
    
    load = max(i + q); % people in quarantine are not getting anyone sick but they are burden on healthcare resources

end


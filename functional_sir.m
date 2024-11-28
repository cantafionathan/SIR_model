% ADAPTED FROM MODEL IN THIS PAPER
% Turkyilmazoglu M. An extended epidemic model with vaccination: Weak-immune SIRVI. 
% Physica A. 2022;598:127429. doi:10.1016/j.physa.2022.127429

%-----MODEL-----%

% 4 x 4 x 4 x 4 factorial. Can try to use generator table on the 4 level
% factors to reduce number of runs

% Consider as a 2^8 factorial with factors 
% A = isolation factor; pseudo-factors with 2 levels A1, A2
% B = vaccination rate; pseudo-factors with 2 levels B1, B2
% C = quarantine duration; pseudo-factors with 2 levels C1, C2
% D = # daily interactions; pseudo-factors with 2 levels D1, D2
%
% Reduce number of runs from 256 to 64 by using generators 
% D1 = A1A2B1B2; D2 = A1A2C1C2 this is a resolution 5 design
% so 2-way interactions will aliased with at worst 3 way interactions and
% we can still estimate those effects (which we need to since A1A2, etc.
% are like main effects)

% use generators to determine levels of num_daily

L = {[1, -1], [1, -1], [1, -1], [1, -1], [1, -1], [1, -1]};
n = length(L);
[L{:}] = ndgrid(L{end:-1:1});
L = cat(n+1,L{:});
L = fliplr(reshape(L,[],n));

A1 = L(1:end, 1); % 1
A2 = L(1:end, 2); % 3
B1 = L(1:end, 3); % 2
B2 = L(1:end, 4); % 5
C1 = L(1:end, 5); % 4
C2 = L(1:end, 6); % 6
D1 = A1.*B1.*A2.*C1; % 7 = 1234
D2 = A1.*B1.*B2.*C2; % 8 = 1256

% compute levels of all the factors that aren't random

isolation = [0, 0.27, 0.4, 0.81];
vac_rate = [0, 0.003, 0.006, 0.009];
quar_dur = [0, 8, 16, 24];
num_daily = [10, 20, 30, 40];

[A, B, C] = ndgrid(isolation, vac_rate, quar_dur);
Av = A(:);
Bv = B(:);
Cv = C(:);

Dv = zeros(length(L), 1);
for i = 1:length(L)
    if (D1(i) == 1) && (D2(i) == 1)
        Dv(i) = num_daily(1);
    elseif (D1(i) == 1) && (D2(i) == -1)
        Dv(i) = num_daily(2);
    elseif (D1(i) == -1) && (D2(i) == 1)
        Dv(i) = num_daily(3);
    else
        Dv(i) = num_daily(4);
    end
end

for i = 1:length(L)
    % randomly generated levels
    prob_spread = 0.1*rand(1)+0.05; % random number between (0.05, 0.15)
    recovery_rate = 1/(14*rand(1) + 7); % random number between (1/7, 1/21)
    vac_eff = 0.2*rand(1) + 0.4; % random number between (0.4, 0.6)
    prob_sympt = 0.6*rand(1); % random number between (0, 0.6)

    load(i) = pandemic(Dv(i), Av(i), prob_spread, recovery_rate, Bv(i), vac_eff, prob_sympt, Cv(i));
end

M = [load, Av, Bv, Cv, Dv];
Labels = ["peak_infected", "social_isolation", "rate_vaccinated", "quarantine_duration", "daily_interactions"];
writetable(array2table(M, 'VariableNames', Labels), 'fractional_pandemic.csv')

% num_daily = [10, 20, 30, 40]; % 4 level factor
% isolation = [0, 0.3, 0.6, 0.9]; % 4 level factor
% prob_spread = [0.05, 0.1, 0.15]; % consider randomizing and not recording
% recovery_rate = [1/7, 1/14, 1/21]; % consider randomizing and not recording
% vac_rate = [0, 0.003, 0.006, 0.009]; % 4 level factor
% vac_eff = [0, 0.005, 0.01, 0.015, 0.02]; % consider randomizing and not recording
% prob_sympt = [0, 0.2, 0.4, 0.6]; % consider randomizing and not recording
% quar_dur = [0, 8, 16, 24]; % 4 level factor

% % [A, B, C, D, E, F, G, H] = ndgrid(num_daily, isolation, prob_spread, recovery_rate, vac_rate, vac_eff, prob_sympt, quar_dur);
% % Av = A(:);
% % Bv = B(:);
% % Cv = C(:);
% % Dv = D(:);
% % Ev = E(:);
% % Fv = F(:);
% % Gv = G(:);
% % Hv = H(:);
% %
% % load = zeros(length(Av), 1);
% % for i = 1:length(Av)
% %     load(i) = pandemic(Av(i), Bv(i), Cv(i), Dv(i), Ev(i), Fv(i), Gv(i), Hv(i));
% % end
% %
% % M = [load, Av, Bv, Cv, 1./Dv, Ev, Fv, Gv, Hv];
% % Labels = ["peak_infected", "pop_density", "social_isolation", "prob_spread", "recovery_duration", "rate_vaccinated", "vaccine_effectivness", "prob_symptomatic", "quarantine_duration"];
% % writetable(array2table(M, 'VariableNames', Labels), 'pandemic.csv')

% [A, B, C] = ndgrid(isolation, vac_rate, quar_dur);
% Av = A(:);
% Bv = B(:);
% Cv = C(:);
% 
% load = zeros(length(Av), 1);
% for i = 1:length(Av)
%     num_daily = round(25*rand(1)+5); % random number between (5, 30) % CONSIDER NOT RANDOMIZING
%     prob_spread = 0.1*rand(1)+0.05; % random number between (0.05, 0.15)
%     recovery_rate = 1/(14*rand(1) + 7); % random number between (1/7, 1/21)
%     vac_eff = 0.2*rand(1) + 0.4; % random number between (0.4, 0.6)
%     prob_sympt = 0.6*rand(1); % random number between (0, 0.6)
%     load(i) = pandemic(num_daily, Av(i), prob_spread, recovery_rate, Bv(i), vac_eff, prob_sympt, Cv(i));
% end
% 
% M = [load, Av, Bv, Cv];
% Labels = ["peak_infected", "social_isolation", "rate_vaccinated", "quarantine_duration"];
% writetable(array2table(M, 'VariableNames', Labels), 'reduced_pandemic.csv')


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
    isolation_factor = isolation; 
    
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
    
    load = max(i);

end


clc;
clear;
close all;

%% Parameters Initialization
n = 100; 
xm = 100; 
ym = 100; 
sink = [50, 50]; 
p = 0.1; 
Eo = 0.5; 
ETX = 50e-9; 
ERX = 50e-9; 
Efs = 10e-12; 
Emp = 0.0013e-12; 
EDA = 5e-9; 
k = 4000; 
do = sqrt(Efs / Emp); 
rounds = 1000; 

nodes = initialize_nodes(n, xm, ym, Eo);

distance_to_sink = calculate_distances(nodes, sink);
death_round = 0; 
death_count = 0;

dead_nodes_history = zeros(1, rounds);
total_energy_history = zeros(1, rounds);
cluster_heads_history = zeros(1, rounds);

for r = 1:rounds
    if death_count == n
        break;
    end

    
    nodes = reset_nodes(nodes);

    % Preallocate cluster_heads with maximum possible size
    cluster_heads = zeros(1, n); 
    ch_index = 0; 

    % Select Cluster Heads
    for i = 1:n
        if nodes.E(i) > 0 && nodes.G(i) == 0
            prob = p / (1 - p * mod(r, round(1 / p)));
            if rand < prob
                nodes.type{i} = 'CH';
                nodes.G(i) = round(1 / p);
                ch_index = ch_index + 1;
                cluster_heads(ch_index) = i;
            end
        end
    end
    cluster_heads = cluster_heads(1:ch_index); 

    nodes = form_clusters(nodes, cluster_heads);

    % Energy Dissipation
    nodes = dissipate_energy(nodes, cluster_heads, distance_to_sink, do, Efs, Emp, ETX, EDA, k);

    prev_death_count = death_count;
    epsilon = 1e-9;
    death_count = sum(double(nodes.E <= epsilon));  
    if death_count > prev_death_count && death_round == 0
        death_round = r;
    end

    % Record data for plots
    dead_nodes_history(r) = round(death_count); 
    total_energy_history(r) = sum(nodes.E);
    cluster_heads_history(r) = length(cluster_heads);

    
    if r <= 10 || r == 220 || death_count == n
        fprintf('Round %d: Dead nodes: %d, Total energy: %.4f, Cluster heads: %d\n', ...
                r, death_count, total_energy_history(r), cluster_heads_history(r));
    end
end


fprintf('First node died at round: %d\n', death_round);
fprintf('Total dead nodes after %d rounds: %d\n', rounds, death_count);

%% Visualization

plot_network(nodes, sink);

% Plot Dead Nodes vs. Rounds
figure;
plot(1:rounds, dead_nodes_history, 'r-', 'LineWidth', 2);
xlabel('Rounds');
ylabel('Number of Dead Nodes');
title('Dead Nodes Over Time');
grid on;

% Plot Total Energy vs. Rounds
figure;
plot(1:rounds, total_energy_history, 'b-', 'LineWidth', 2);
xlabel('Rounds');
ylabel('Total Energy (J)');
title('Total Energy Depletion Over Time');
grid on;

% Plot Cluster Heads vs. Rounds
figure;
plot(1:rounds, cluster_heads_history, 'g-', 'LineWidth', 1);
xlabel('Rounds');
ylabel('Number of Cluster Heads');
title('Cluster Heads Over Time');
grid on;

%% Functions

% Function to initialize sensor nodes
function nodes = initialize_nodes(n, xm, ym, Eo)
    nodes.x = rand(1, n) * xm;  
    nodes.y = rand(1, n) * ym;  
    nodes.E = Eo * ones(1, n);  
    nodes.type = repmat({'N'}, 1, n);  
    nodes.G = zeros(1, n);  
    nodes.cluster = zeros(1, n); 
end

% Function to calculate distances of nodes to the sink
function distances = calculate_distances(nodes, sink)
    distances = sqrt((nodes.x - sink(1)).^2 + (nodes.y - sink(2)).^2);
end

% Function to reset cluster head status and decrement timers
function nodes = reset_nodes(nodes)
    nodes.G = max(nodes.G - 1, 0);
    nodes.type = repmat({'N'}, 1, length(nodes.x));
end

% Function to assign nodes to the nearest cluster head
function nodes = form_clusters(nodes, cluster_heads)
    n = length(nodes.x);
    for i = 1:n
        if strcmp(nodes.type{i}, 'N') && nodes.E(i) > 0
            min_dist = inf;
            for j = cluster_heads
                d = sqrt((nodes.x(i) - nodes.x(j))^2 + (nodes.y(i) - nodes.y(j))^2);
                if d < min_dist
                    min_dist = d;
                    nodes.cluster(i) = j;
                end
            end
        end
    end
end

% Function to dissipate energy during communication
function nodes = dissipate_energy(nodes, cluster_heads, distance_to_sink, do, Efs, Emp, ETX, EDA, k)
    n = length(nodes.x);
    for i = cluster_heads
        num_members = sum(nodes.cluster == i);
        nodes.E(i) = max(nodes.E(i) - (ETX + EDA) * k * num_members + Emp * k * (distance_to_sink(i)^4), 0);
        for j = 1:n
            if nodes.cluster(j) == i && nodes.E(j) > 0
                d = sqrt((nodes.x(i) - nodes.x(j))^2 + (nodes.y(i) - nodes.y(j))^2);
                if d < do
                    nodes.E(j) = max(nodes.E(j) - (ETX + EDA) * k + Efs * k * (d^2), 0);
                else
                    nodes.E(j) = max(nodes.E(j) - (ETX + EDA) * k + Emp * k * (d^4), 0);
                end
            end
        end
    end
end

% Function to plot the network
function plot_network(nodes, sink)
    figure;
    hold on;
    for i = 1:length(nodes.x)
        if nodes.E(i) > 0
            plot(nodes.x(i), nodes.y(i), 'go', 'MarkerSize', 8, 'LineWidth', 1.5);
        else
            plot(nodes.x(i), nodes.y(i), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
        if strcmp(nodes.type{i}, 'CH')
            plot(nodes.x(i), nodes.y(i), 'bs', 'MarkerSize', 10, 'LineWidth', 2);
        end
    end
    plot(sink(1), sink(2), 'kp', 'MarkerSize', 12, 'LineWidth', 2);
    legend('Alive Nodes', 'Dead Nodes', 'Cluster Heads', 'Sink');
    title('WSN Node Deployment and Status');
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    grid on;
    hold off;
end
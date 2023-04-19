% Plot values in optimization process

%% Figure1: objective function
figure(1)
clf
width =  500;
height = 380;
set(gcf, "Position", [1000 500 width height])
fontsize(gca,15,"pixels")
box;
hold on
xlim([0,30])
ylim([0,10])
xlabel('Iteration')
ylabel('Objective Function Value')

% objective function at each iteration
load("results.mat")
plot(0:num_iter-1, J_list(:,1), 'b-x', 'DisplayName','J' ) 
plot(0:num_iter-1, J_list(:,3), 'm-x','DisplayName','J_C' ) 

% with fixed rho
load("results_fixed_rho.mat")
plot(0:num_iter-1, J_list(:,1), 'b--', 'DisplayName','J (with fixed \rho)' ) 
plot(0:num_iter-1, J_list(:,3), 'm--', 'DisplayName','J_C (with fixed \rho)' ) 

legend('Location','northeast')
saveas(gcf,'plot/iteration_objective','epsc')

%% Figure2: parameter rho
figure(2);
clf
width =  500;
height = 380;
set(gcf, "Position", [1000 60 width height])
fontsize(gca,15,"pixels")
box;
hold on
xlim([0,30])
ylim([1,2.5])
xlabel('Iteration')
ylabel('Design Parameter \rho')

load("results.mat")
plot(0:num_iter-1, rho_list, 'r-x' ) 
saveas(gcf,'plot/iteration_parameter','epsc')
clear; close all; clc;
%% Genetic Programming Test
%the following script contains the main script and functions to set up and
%run a genetic programming process to optimize a control law for an
%inverted pendulum system

%define system constants in a structure
P.m = 2;        % Pendulum bob mass [kg]
P.M = 10;       % Cart mass [kg]
P.g = 9.81;     % gravitational acceleration [m/s/s]
P.L = .66;      % Pendulum length [m]

inputs = {'x(1)','x(2)','x(3)','x(4)', 'x(5)', 'x(6)','1','2','-1','0.5',...
          '7','1/5','-1.85','1/10', '3.14','-3.14','exp(1)','exp(1.5)',...
          'tan(x(3))', 'cot(x(3))', 'asin(x(3))', 'acos(x(3))',...
          'cos(x(3))','cos(x(4))','sin(x(3))','sin(x(4))'};
%           'sign(x(1))','sign(x(2))','sign(x(3))','sign(x(4))', 'sign(x(5))', 'sign(x(6))',...

ops = {'plus','minus','times'};

% inputs = {'x(1)','x(2)','x(3)','x(4)','1','2','-1','0.5', 'x(5)', 'x(6)'};
% ops = {'plus','minus','times', 'cos(plus','sin(plus','cos(minus',...
%     'sin(minus', 'sign(plus','sign(minus','cos(times','sin(times',...
%     'sign(times'};

%Create an initial population
%these create random trees with varying numbers node levels
i1 = treePop(10,3, inputs, ops);
i2 = treePop(10,4, inputs, ops);
i3 = treePop(10,5, inputs, ops);
i4 = treePop(10,6, inputs, ops);

%Create gain matricies, for statefeedback controllers
K1 = [0.5  2.5  200  30]; 
K2 = [10   35  500   75];
K3 = [0.001 0.03  120  15];
%create trees that represent the control law u = -K*x for the gain matrices
t1 = stateFBtree(K1);
t2 = stateFBtree(K2);
t3 = stateFBtree(K3);

%construct the combined initial population
initPop = [i1 , i2 , i3 , i4,{t1 t2 t3}];
%% Parameters to Change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define reproduction statistics
crossover = 0.2; mutation = 0.6; 
eliteVals= 3; numGens = 20;
%Define simulation timespan when simulating and evaluating control cost
tspan = [0 10];
%Define system initial conditions
x0 = [0 0 pi/10 0 0 0];
%define cost matricies
statePenalty = [0 30 100 20 0 0];
Q = diag(statePenalty);
R = 0.001;
%Set timeout value for integration. Makes code not get stuck solving bad
%simulations
timeout = 2; %seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Run the genetic algorithm function for this problem
[bestTree,bestFcn,bestCost] = treeGA(initPop, crossover , mutation ,...
    eliteVals , numGens , tspan , x0 , P, Q ,R,timeout, inputs, ops);
 
%an odeset option must be specified because I am too lazy to make a special
%function, so I have chosen Refine
options = odeset('Refine',4);
%simulate the top result
[tout , xout , uout] = fcnSim(bestFcn,tspan,x0,P,options);
figure;
subplot(3,1,1)
plot(tout,xout(:,1))
ylabel('Cart Pos , x [m]')
subplot(3,1,2)
plot(tout,xout(:,3))
ylabel('Pendulum Ang , \theta [rad]')
subplot(3,1,3)
plot(tout,uout)
ylabel('Input Force , f [N]')
xlabel('Time , t [seconds]')
title('GA Optimum Control Response')

%% now plot the controller for longer timespan and different IC's
tspan = [0 50];
x0 = [0 0 pi/10 0 0 0];
[tout , xout , uout] = fcnSim(bestFcn,tspan,x0,P,options);
figure;
subplot(3,1,1)
plot(tout,xout(:,1))
ylabel('Cart Pos , x [m]')
subplot(3,1,2)
plot(tout,xout(:,3))
ylabel('Pendulum Ang , \theta [rad]')
subplot(3,1,3)
plot(tout,uout)
ylabel('Input Force , f [N]')
xlabel('Time , t [seconds]')
title('GA Optimum Control Response (extended)')

%display the winning tree
disp('The best tree is ...')
disp(bestTree.tostring)
disp('fuk ya')

%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ipop ] = treePop(PopSz,nodeLvls, inputs, ops)
%creates a population of function trees
%ipop is initial population of trees
%PopSz is number of trees in the population
%nodeLvls is the number of binary node levels. 
% inputs = {'x(1)','x(2)','x(3)','x(4)','1','2','-1','0.5', 'x(5)', 'x(6)'};
% ops = {'plus','minus','times', 'cos(plus','sin(plus','cos(minus',...
%     'sin(minus', 'sign(plus','sign(minus','cos(times','sin(times',...
%     'sign(times'};
ipop = cell(1,PopSz);

%determine number of non leaf nodes
nodes = 2^(nodeLvls) - 1;

for n = 1:PopSz
    
    t = tree(ops{randi(length(ops))});
    for element = 1:nodes
        if element < 2^(nodeLvls - 1)
            t = t.addnode(element,ops{randi(length(ops))});
            t = t.addnode(element,ops{randi(length(ops))});
        else
            t = t.addnode(element,inputs{randi(length(inputs))});
            t = t.addnode(element,inputs{randi(length(inputs))});
        end    
    end
ipop{n} = t;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ str ] = tree2str(tree)
%get all nodes of tree
iterator = tree.nodeorderiterator;
%if still multiple nodes in tree, 
    if length(iterator) ~= 1
        for idx = 1
            
            str = char([tree.Node{iterator(idx)},'(']);
            child = tree.getchildren(iterator(idx));
            term1 = tree2str(tree.subtree(child(1)));
            term2 = tree2str(tree.subtree(child(2)));
            if(length(str)>7)
            str = char([str , term1 ,',' term2 , '))']); 
            else
            str = char([str , term1 ,',' term2 , ')']);   
          
            end
        end
    else
    str = char(tree.Node{iterator});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func] = tree2func(tree)
%computes a function from a tree 
str = tree2str(tree);
str = char(['@(x)',str]);
func = str2func(str);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ntree1, ntree2] = treeCrossover(tree1,tree2)
    % first, define the number of nodes that can be chopped
    maxnode1 = length(tree1.Parent);
    maxnode2 = length(tree2.Parent);
    % next randomly decide which node to chop/graft beginning at node 2
    node1 = randi([2 maxnode1]);
    node2 = randi([2 maxnode2]);
    % determine the parent which this node comes from. This is useful when
    % grafting the new section to the tree
    parent1 = tree1.Parent(node1);
    parent2 = tree2.Parent(node2);
    % Next, find the portion that will be grafted/chopped from the tree
    sec1 = tree1.subtree(node1);
    sec2 = tree2.subtree(node2);
    %  chop the where the new subtree will be grafted, then graft
    ntree1 = tree1.chop(node1);
    ntree1 = ntree1.graft(parent1,sec2);
    ntree2 = tree2.chop(node2);
    ntree2 = ntree2.graft(parent2,sec1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ntree = treeMutation(tree, inputs, ops)
%     inputs = {'x(1)','x(2)','x(3)','x(4)','1','2','-1','0.5', 'x(5)', 'x(6)'};
%     ops = {'plus','minus','times', 'cos(plus','sin(plus','cos(minus',...
%         'sin(minus', 'sign(plus','sign(minus','cos(times','sin(times',...
%         'sign(times'};
    % find the max node that can be mutated
    maxnode = tree.getchildren(tree.Parent(end));
    maxnode = maxnode(end);
    % next, randomly decide which node to mutate
    mutnode = randi([1 maxnode]);
    
    % if the mutation node is a leaf, mutate it with a random input
    if tree.isleaf(mutnode)
        ntree = tree.set(mutnode,inputs{randi(length(inputs))});
    else
        % change the node with an operation if it is not a leaf
        ntree = tree.set(mutnode,ops{randi([1 3])});
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function repMethod = sexStyle(mutPercent, crossPercent)
% The following function decides the method to reproduce when the
% percentages are given in decimal form.
% The output values are the following:
% 1: Mutation
% 2: Crossover
% 3: Replication
    randNum = rand();
    if randNum < mutPercent
        repMethod = 1;
    elseif randNum < (mutPercent + crossPercent)
        repMethod = 2;
    else
        repMethod = 3;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = costCalc(t,x,u,Q,R)
% The following function calculates the cost for the simulation based on
% the Q and R cost matrices
% Inputs:
% x = State Vector from ode15s
% u = control effort from simulation
% Q = State cost vector (make it the same as the LQR)
% R = Control effort cost vector (also make same as LQR)
    [rows, cols] = size(x);
    cost_vec = zeros(rows,1);
    for i = 1:rows
        cost_vec(i) = x(i,:)*Q*x(i,:)' + u(i)'*R*u(i);
    end
    cost = trapz(t,cost_vec); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [percentSelect, sortedtreeList] = selectionPercent(treeList, costList)
% The following function claculates the fitness and creates the selection
% percentage array based on the costList after all simulations have been
% run
    % Sort the costList from least to greatest
    [costList, idx] = sort(costList,'ascend');
    % Next, sort the treelist to correspod to the costList
    sortedtreeList = treeList(idx);
    % Calculate the fitness
    fitnessList = sum(costList)./costList;
    % calulate the selection percentage matrix
    percentSelect = fitnessList./sum(fitnessList);
    
    % This final loop organizes the percentSelect array such that the
    % function rand() can be used to determine which individual to use
    percentSelect = [0 percentSelect];
    for i = 2:length(percentSelect)
        percentSelect(i) = percentSelect(i) + percentSelect(i-1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indiv = indivSelect(sortedtreeList,percentSelect)
% percentSelect is the selection percentage matrix
% sortedtreeList is the treeList that corresponds to the selection matrix
    perc = rand();
    idx = find(percentSelect<perc);
    indiv = sortedtreeList(idx(end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tout, xout , uout] = fcnSim(fcn_hand , tspan , x0 ,P ,options)
%simulates the system with given control law function
%outputs state vector and input for entire time span
%define constants from parameter stucture
m = P.m;
M = P.M;
L = P.L;
g = P.g;
%create function handle for entire system
sys_hand = @(t,x)cartpend(t,x,m,M,L,g,fcn_hand);
[tout , xout] = ode15s(sys_hand,tspan,x0,options);
uout = zeros(length(tout),1);
for n = 1:length(tout)
  uout(n,1) = fcn_hand(xout(n,:));  
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bestTree , bestFcn, bestCost] = treeGA(initPop , crossover ,...
    mutation , eliteVals, genNum, tspan , x0, P , Q, R , timeout,...
    inputs, ops)

%genetic algorthims for the control laws represented by trees

pop = initPop;
for gen = 1:genNum
    disp(['Generation', num2str(gen)])
    %determine sim data for each member of the population 
    for k = 1:length(pop)
       %get function from current tree
       u_fcn = tree2func(pop{k});
       %reset warnings 
       %warnings are used to assign infinite cost to controllers that are
       %bad, ie singularities in system or just take a long time to
       %simulate
       lastwarn('','');
       %start timer
       tstart = tic;
       %set options for current timer
       options = odeset('Events',@(t,y)myEvent(t,y,tstart,timeout));
       [tout{k}, xout{k} , uout{k}] = fcnSim(u_fcn, tspan , x0 , P,options);
       [warnMsg , warnID] = lastwarn();
       if isempty(warnMsg)
           %evaluate cost of control law
           cost(k) = costCalc(tout{k} , xout{k} , uout{k}, Q , R);
       else
           cost(k) = inf;
           warning('Controller given Inf cost')
       end
   end
   %create percentage selection array and sorted tree lest
   [percentSelect, sortedtreeList] = selectionPercent(pop, cost);
   
   for k = 1:length(pop)
       %replicate elite values
       if k <= eliteVals
           pop1{k} = sortedtreeList{k};    
       else
           %Choose a reproduction type
           repMethod = sexStyle(mutation , crossover);
           switch repMethod
               case 1 %Mutation
                   tree = indivSelect(sortedtreeList,percentSelect);
                   tree = tree{1};
                   pop1{k} = treeMutation(tree, inputs, ops);
               case 2 %Crossover
                   tree1 = indivSelect(sortedtreeList,percentSelect);
                   tree1 = tree1{1};
                   tree2 = indivSelect(sortedtreeList,percentSelect);
                   tree2 = tree2{1};           
                   [ntree{1} , ntree{2}] = treeCrossover(tree1,tree2);
                   %randomly choose one of the children to be surviving
                   %offspring
                   pop1{k} = ntree{randi(2)};              

               case 3 %Replication
                   tree = indivSelect(sortedtreeList,percentSelect);
                   pop1{k} = tree{1};

           end
        end
   end
   pop = pop1;  
end
if gen == genNum
    %determine best tree from last generation
    for k = 1:length(pop)
           %get function from current tree
           u_fcn = tree2func(pop{k});
           %reset warnings
           lastwarn('','');
           %start timer
           tstart = tic;
           %set options for current timer
           options = odeset('Events',@(t,y)myEvent(t,y,tstart,timeout));
           [tout{k}, xout{k} , uout{k}] = fcnSim(u_fcn, tspan , x0 ,P,options);
           [warnMsg , warnID] = lastwarn();
           if isempty(warnID)
               %evaluate cost of control law
               cost(k) = costCalc(tout{k} , xout{k} , uout{k}, Q , R);
           else
               cost(k) = Inf;
               warning('Controller given Inf cost')
           end
    end
   %create percentage selection array and sorted tree lest
   [percentSelect, sortedtreeList] = selectionPercent(pop, cost);
   bestTree = sortedtreeList{1};
   bestFcn = tree2func(bestTree);
   
   [tout, xout , uout] = fcnSim(bestFcn, tspan , x0, P,options);
   
   bestCost = costCalc(tout,xout,uout,Q,R);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,isterminal,direction] = myEvent(t,y,tstart,timeout)
%event function that limits the amount of time odesolver can run
value = toc(tstart) < timeout;
isterminal = 1; %halt integration
direction = []; 
if ~value
    w = sprintf('Elapsed time exceeds %0.2f seconds, integration stopped',timeout);
    warning(w);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [treeK] = stateFBtree(K)
%creates a control law tree for a given state feedback matrix K
treeK = tree('plus');
treeK = treeK.addnode(1,'plus'); treeK = treeK.addnode(1,'plus');
treeK = treeK.addnode(2,'times'); treeK = treeK.addnode(2,'times');
treeK = treeK.addnode(3,'times'); treeK = treeK.addnode(3,'times');
treeK = treeK.addnode(4,'x(1)'); treeK = treeK.addnode(4,num2str(K(1))); 
treeK = treeK.addnode(5,'x(2)'); treeK = treeK.addnode(5,num2str(K(2))); 
treeK = treeK.addnode(6,'x(3)'); treeK = treeK.addnode(6,num2str(K(3))); 
treeK = treeK.addnode(7,'x(4)'); treeK = treeK.addnode(7,num2str(K(4)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

m = 2;        % Pendulum bob mass [kg]
M = 10;       % Cart mass [kg]
g = 9.81;     % gravitational acceleration [m/s/s]
L = .66;      % Pendulum length [m]


%test ga function
i1 = treePop(10,2);
i2 = treePop(10,3);
i3 = treePop(10,4);
initPop = [i1 , i2 , i3];


crossover = 0.8;
mutation = 0.2;
NumElites = 3;
genNum = 20;
tsim = 5;

[bestTree,bestFcn] = treeGA(initPop, crossover ,mutation , NumElites , genNum, tsim);
sys = @(t,x)cartpend(x,m,M,L,g,bestFcn);
[tout , xout] = ode45(sys,[0 tsim],[0 0 pi/10 0]);
figure;
hold on
plot(tout,xout(:,1))
plot(tout,xout(:,3))
legend('x','theta')
%%
function [ ipop ] = treePop(PopSz,nodeLvls)
%creates a population of function trees
%ipop is initial population of trees
%PopSz is number of trees in the population
%nodeLvls is the number of binary node levels. 
inputs = {'x(1)','x(2)','x(3)','x(4)','1','2','-1','0.5'};
ops = {'plus','minus','times'};
ipop = cell(1,PopSz);

%determine number of non leaf nodes
nodes = 2^(nodeLvls) - 1;

for n = 1:PopSz
    t = tree(ops{randi([1 3])});
    for element = 1:nodes
        if element < 2^(nodeLvls - 1)
            t = t.addnode(element,ops{randi([1 3])});
            t = t.addnode(element,ops{randi([1 3])});
        else
            t = t.addnode(element,inputs{randi(length(inputs))});
            t = t.addnode(element,inputs{randi(length(inputs))});
        end    
    end
ipop{n} = t;
end
end

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
        str = char([str , term1 ,',' term2 , ')']); 
    end
else
    str = char(tree.Node{iterator});
end
end

function [func] = tree2func(tree)
%computes a function from a tree 
str = tree2str(tree);
str = char(['@(x)',str]);
func = str2func(str);
end

function [ntree1, ntree2] = treeCrossover(tree1,tree2)
    % first, define the number of nodes that can be chopped
    maxnode1 = max(tree1.Parent);
    maxnode2 = max(tree2.Parent);
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

function ntree = treeMutation(tree)
    inputs = {'x(1)','x(2)','x(3)','x(4)','1','2','-1','0.5'};
    ops = {'plus','minus','times'};
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

function cost = costCalc(t,x,u,Q,R)
% The following function calculates the cost for the simulation based on
% the Q and R cost matrices
% Inputs:
% x = State Vector from ode45
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

function indiv = indivSelect(sortedtreeList,percentSelect)
% percentSelect is the selection percentage matrix
% sortedtreeList is the treeList that corresponds to the selection matrix
    perc = rand();
    idx = find(percentSelect<perc);
    indiv = sortedtreeList(idx(end));
end

function [tout, xout , uout] = fcnSim(fcn_hand , tspan , x0 ,options, M , m , L , g)
%simulates the system with given control law function
%outputs state vector and input for entire time span



sys_hand = @(t,x)cartpend(x,m,M,L,g,fcn_hand);
[tout , xout] = ode45(sys_hand,tspan,x0,options);
uout = zeros(length(tout),1);
for n = 1:length(tout)
  uout(n,1) = fcn_hand(xout(n,:));  
end


end

function [bestTree , bestFcn] = treeGA(initPop , crossover , mutation , eliteVals, genNum, tsim)
%genetic algorthims for the control laws represented by trees
%constants I dont want to put as inputs right now
m = 2;        % Pendulum bob mass [kg]
M = 10;       % Cart mass [kg]
g = 9.81;     % gravitational acceleration [m/s/s]
L = .66;      % Pendulum length [m]
pen = [0, 500, 100, 500];  Q = diag(pen); 
R = .005; 

pop = initPop;
for gen = 1:genNum
   %determine sim data for each member of the population 
   for k = 1:length(pop)
       %get function from current tree
       u_fcn = tree2func(pop{k});
       %reset warnings
       lastwarn('','');
       %start timer
       tstart = tic;
       %set options for current timer
       options = odeset('Events',@(t,y)myEvent(t,y,tstart));
       [tout{k}, xout{k} , uout{k}] = fcnSim(u_fcn, [0 tsim] , [0 0 pi/10 0],options , M , m , L , g );
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
               pop1{k} = treeMutation(tree);
           case 2 %Crossover
               tree1 = indivSelect(sortedtreeList,percentSelect);
               tree1 = tree1{1};
               %ensure tree has enough nodes for crossover
               while max(tree1.Parent) < 2
                   tree1 = indivSelect(sortedtreeList,percentSelect);
                   tree1 = tree1{1};
               end
               tree2 = indivSelect(sortedtreeList,percentSelect);
               tree2 = tree2{1};
               %ensure tree has enough nodes for crossover
               while max(tree2.Parent) < 2
                   tree2 = indivSelect(sortedtreeList,percentSelect);
                   tree2 = tree2{1};
               end
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
%determine best tree from last generation
for k = 1:length(pop)
       %get function from current tree
       u_fcn = tree2func(pop{k});
       %reset warnings
       lastwarn('','');
       %start timer
       tstart = tic;
       %set options for current timer
       options = odeset('Events',@(t,y)myEvent(t,y,tstart));
       [tout{k}, xout{k} , uout{k}] = fcnSim(u_fcn, [0 tsim] , [0 0 pi/10 0],options , M , m , L , g );
       [warnMsg , warnID] = lastwarn();
       if isempty(warnID)
       %evaluate cost of control law
       cost(k) = costCalc(tout{k} , xout{k} , uout{k}, Q , R);
       else
       cost(k) = inf;
       warning('Controller given Inf cost')
       end
end
   %create percentage selection array and sorted tree lest
   [percentSelect, sortedtreeList] = selectionPercent(pop, cost);
   bestTree = sortedtreeList{1};
   bestFcn = tree2func(bestTree);
end


function [value,isterminal,direction] = myEvent(t,y,tstart)
%event function that limits the amount of time odesolver can run
value = toc(tstart) < 5;
isterminal = 1; %halt integration
direction = []; 
if value == 0
    warning('Elapsed time exceeds 5 seconds, integration stopped')
end
end

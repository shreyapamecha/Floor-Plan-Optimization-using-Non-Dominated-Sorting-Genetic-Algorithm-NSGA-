% Assignment 3

% Non-Domination Sorting Genetic Algorithm - II (NSGA-II)
% Optimization Problem: Maximizing the area of laboratory (f1) | Minimizing the
% total cost of designing teh workspace (f2)

%clear all;
%clear;
%clc;

particles = 10; %number of particles
CP = 0.8; %crossover probability
equip_room = 450; %unit pricing of equipment room
lab_room = 600; %unit pricing of laboratory
seating_room = 300; %unit pricing of meeting place
pop_for_graph=[]; %matrix to store functional values and rank, so that it can be used to draw pareto fronts later

% 10 <= x <= 15 
% 7 <= y <= 15

e_length = 0; %initial length of the equipment room
e_width = 0; %initial width of the equipment room
l_length = 0; %initial length of the laboratory
l_width = 0; %initial width of the laboratory
s_length = 25; % Fixed length of the seating/meeting place
s_width = 0; %initial width of the seating place

%for binary values of x and y (8 columns)
pop = hardlim(rand(particles,8)-0.5); %1 for values greater than zero and 0 for values less than and equal to zero.

%Calculating Functional Values: area of the laboratory (maximize) and the cost for
%designing the workplace (minimize)
pop = functional_values(pop,1,particles,equip_room,lab_room,seating_room);

%deciphering Ranks of all the particles and storing them in the 13th column
%of 'pop' matrix
[pop,total_rank,rank] = RankList(pop);

%Storing functional values and ranks of all the iterations in this matrix
%so that it can be used to draw pareto fronts later
pop_for_graph=[pop_for_graph;pop(:,11:13)];

for iteration=1:30
    % Crossover
    crossover=CP*particles; 

    %Determine if sorting of particles within a rank is required of not
    %rank_i => the rank which requires sorting using Crowding Distance
    %method
    %if crowding_distance_required==1, then it's required otherwise it's
    %not
    [crowding_distance_required,rank_i] = Crowding_reqd(total_rank,rank,crossover);

    if crowding_distance_required==1
        %Calculates Crowding distance & stores the sorted particles in the
        %'pop' matrix
        pop=CalculateCrowdingDistance(pop,rank_i);
    end

    % Adding 8 Children to the existing 'pop' matrix
    pop = AddChildren(pop,CP,particles);

    % Calculating functional values of the children particles
    pop = functional_values(pop,particles+1,particles+CP*particles,equip_room,lab_room,seating_room);
    % Determining the new ranks for all teh particles
    [pop,total_rank,rank] = RankList(pop);

    %From here we want to take out top 10 particles for further process
    crossover=particles;
    [crowding_distance_required,rank_i] = Crowding_reqd(total_rank,rank,crossover);

    if crowding_distance_required==1
        %Calculates Crowding distance and stores the sorted particles in
        %'pop' matrix
        pop=CalculateCrowdingDistance(pop,rank_i);
    end

    %Removing the rest 8 and keeping the strong 10 particles
    pop=pop(1:particles,:);
    
    %Storing functional values and ranks of all the iterations in this matrix
    %so that it can be used to draw pareto fronts later
    pop_for_graph=[pop_for_graph;pop(:,11:13)];
end

%Sort Rows with repect to the 1st functional values
pop=sortrows(pop,11);


% User-Defined Functions

% This function randomly selects 2 out of 8 top particles (4 times) and generates
% children using crossover
function pop=AddChildren(pop,CP,particles)
    size_of_pop=size(pop);
    rows=size_of_pop(1);
    crossover=CP*particles;
    cross=randperm(crossover); 
    next_row=rows;
    for i=1:crossover/2
        parent1=pop(cross(2*i-1),1:8);
        parent2=pop(cross(2*i),1:8);
        cross_pos=randi(7);
        child1=[parent1(1:cross_pos) parent2(cross_pos+1:8)];
        child2=[parent2(1:cross_pos) parent1(cross_pos+1:8)];
        next_row=next_row+1;
        pop(next_row,1:8)=child1(:,:);
        next_row=next_row+1;
        pop(next_row,1:8)=child2(:,:);
    end

end

%Calculates Crowding distance and stores the sorted particles in
%'pop' matrix
function pop = CalculateCrowdingDistance(pop,rank_i)
    size_of_pop=size(pop);
    rows=size_of_pop(1);
    front=[];
    crowding_distance=[];
    k=0; %the number of rows with rank = rank_i
    for j=1:rows
        if pop(j,13)==rank_i
           k=k+1;
           front=[front;pop(j,:)];
           front(k,11)= -front(k,11); % so now we have to minimize both the functions;
           crowding_distance=[crowding_distance,0];
        end
    end
   
    size_of_front = size(front);
    rows_of_front = size_of_front(1);
    
    min_function1 = min(front(:,11));
    min_function2 = min(front(:,12));
    pareto_front=[];
    index_front = [];
    
    for j=1:rows_of_front
       if front(j,11)==min_function1 
          pareto_front=[pareto_front;front(j,:)];
          front(j,:)=[];
          break
       end
    end
    
    size_of_front = size(front);
    rows_of_front = size_of_front(1);
    
    for j=1:rows_of_front
       if front(j,12)==min_function2 
          pareto_front=[pareto_front;front(j,:)];
          front(j,:)=[];
          break
       end
    end
    
    size_of_front = size(front);
    rows_of_front = size_of_front(1);
    
    for j=1:rows_of_front
        pareto_front=[pareto_front;front(j,:)];
    end
    
    size_of_pareto_front=size(pareto_front);
    rows_of_pareto_front=size_of_pareto_front(1);
    
    for j=1:rows
        if pop(j,13)==rank_i
            start_index_of_rank_i = j;
            disp(start_index_of_rank_i)
            break
        end
    end
    
    if rows_of_pareto_front>3
        %pareto_front_swap_second_last = pareto_front;
        %rows_to_swap = [2,rows_of_pareto_front];
        
        %pareto_front_swap_second_last(rows_to_swap,:)=pareto_front_swap_second_last(rows_to_swap([2,1]),:);
        
        crowding_distance(1)=Inf;
        crowding_distance(k)=Inf;
    
        sorted1=sortrows(pareto_front,11);
        sorted2=sortrows(pareto_front,12);
        
        % Crowding Distance Formula
        for j=2:k-1 
           crowding_distance(j)=crowding_distance(j)+((sorted1(j+1,11)-sorted1(j-1,11))/(max(sorted1(:,11))-min(sorted1(:,11))));
        end
        
        for j=2:k-1
            crowding_distance(j)=crowding_distance(j)+((sorted2(j+1,12)-sorted2(j-1,12))/(max(sorted2(:,12))-min(sorted2(:,12))));
        end
        %
        
        for j=1:rows_of_pareto_front
            pareto_front(j,11)=-pareto_front(j,11);
            sorted1(j,11)=-sorted1(j,11);
        end    
        
        for j=1:2
           pop(start_index_of_rank_i,:)=pareto_front(j,:);
           start_index_of_rank_i=start_index_of_rank_i+1;
        end
        
        crowding_distance(1)=0;
        crowding_distance(rows_of_pareto_front)=0;
               
        for j=3:rows_of_pareto_front
            index_max = find(crowding_distance==max(crowding_distance(2:rows_of_pareto_front-1)));%finding index which has the maximum element
            size_of_index_max=size(index_max);
            cols_index_max=size_of_index_max(2);
            if cols_index_max>1
                pop(start_index_of_rank_i,:)=sorted1(index_max(1),:);
                crowding_distance(index_max(1))=0;
            else
                pop(start_index_of_rank_i,:)=sorted1(index_max,:);
                crowding_distance(index_max)=0;
            end
            start_index_of_rank_i=start_index_of_rank_i+1;
        end    
    end

    
    if rows_of_pareto_front==3 || rows_of_pareto_front==2
        for j=1:rows_of_pareto_front
            pareto_front(j,11)=-pareto_front(j,11);
        end 
    end
    if rows_of_pareto_front==3  
        pareto_front=sortrows(pareto_front,11);
        pop(start_index_of_rank_i,:)=pareto_front(1,:);
        pop(start_index_of_rank_i+1,:)=pareto_front(3,:);
        pop(start_index_of_rank_i+2,:)=pareto_front(2,:);
    end

end

% deciding if we need to calculate crowding distance or not
function [crowding_distance_required,rank_i] = Crowding_reqd(total_rank,rank,crossover)
    for i=1:rank
        crossover=crossover-total_rank(i);
        if crossover==0
            crowding_distance_required=0;
            break
        elseif crossover<0
            % To whatever value (i) it has stopped, there we have to apply crowding distance sorting within a pareto front or rank   
            crowding_distance_required=1;
            break 
        end
    end
    rank_i=i; % The rank of which crowding distance needs to be calculated
end

% Giving ranks to the particles and finding Pareto Fronts
function [pop,total_rank,rank] = RankList(pop)
    rank=1;

    pop(:,13)=0;
    cut_column = pop(:,13);
    size_matrix=size(pop);
    rows = size_matrix(1);
    total_rank=[]; %a list showing total number of ranks
  
    % This loop runs until every particle is assigned a rank
    while(any(cut_column==0))
        no_rank_yet = [];
        total=0;
        
        for i=1:rows
            if(pop(i,13)==0)
                no_rank_yet=[no_rank_yet,i];
            end
        end
        
        size_no_rank = size(no_rank_yet);
        cols1 = size_no_rank(2);
        
        for i=1:cols1
       
            ith_row = no_rank_yet(i);
            both_not_satisfied = 0;
            for j=1:cols1
                if no_rank_yet(i)==no_rank_yet(j)
                    continue
                end
                jth_row = no_rank_yet(j);
               
                if (pop(ith_row,11)<pop(jth_row,11)) && (pop(ith_row,12)>pop(jth_row,12))
                    both_not_satisfied = 1; 
                    break
                end
               
            end
            
            if both_not_satisfied==0
                pop(ith_row,13)=rank;
                total=total+1;
            end
        end
        
        cut_column=pop(:,13);
        total_rank=[total_rank,total]; %Stores the total particles which has a particular rank
        rank=rank+1;
    end
    
    pop=sortrows(pop,13); % Sort Rows of pop
    
end

% Calculating the functional values and storing them in the main matrix 'pop'
function pop = functional_values(pop,start,stop,equip_room,lab_room,seating_room)
    for i=start:stop
        pop(i,9)=bin2dec(num2str(pop(i,1:4)));
        pop(i,10)=bin2dec(num2str(pop(i,5:8)));
    
        while (pop(i,9)<10)
            pop(i,1:4)=hardlim(rand(1,4)-0.5);
            pop(i,9)=bin2dec(num2str(pop(i,1:4)));
        end
    
        while (pop(i,10)<7)
            pop(i,5:8)=hardlim(rand(1,4)-0.5);
            pop(i,10)=bin2dec(num2str(pop(i,5:8)));
        end
    
        l_length = pop(i,9);
        l_width = 25 - pop(i,10);
        e_length = 25 - pop(i,9);
        e_width = pop(i,10);
        s_length = 25;
        s_width = pop(i,9);
    
        pop(i,11)= l_length*l_width;
        pop(i,12) = (equip_room*e_length*e_width)+(lab_room*l_length*l_width)+(seating_room*s_length*s_width);
    
        pop(i,13) = 0; % Initializing the RANK
    end

end


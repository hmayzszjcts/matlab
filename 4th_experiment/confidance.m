%% load
clc;  clear all;  close all;
load DATA_RFFP_MAY15_LTE_WLAN_7
data = DATA_RFFP_MAY15_LTE_WLAN_7;

%% wrong xy Position delete
%non_zero_xy= find(data(:,17)~= 0& data(:,18)~= 0);
non_zero_xy= find(data(:,17)~= 0 & data(:,18)~= 0 & data(:,7)~= 0  & data(:,8)~= 0 & data(:,9)~= 0 & data(:,17)<= 6.488201490499051e+05 & data(:,18)<= 6.821629319075700e+06);
updated_data = data(non_zero_xy,:);
%min(updated_data(:,3)); max(updated_data(:,3));

%% Grid dividation

% x_max = max(updated_data(:,17));
% x_min = min(updated_data(:,17));
% 
% y_max=max(updated_data(:,18));
% y_min=min(updated_data(:,18));

x_max =  6.488201490499051e+05; % max(updated_data(:,17));
x_min = min(updated_data(:,17));


y_max= 6.821629319075700e+06; %max(updated_data(:,18));
y_min=min(updated_data(:,18));



%x_max-x_min; y_max-y_min;

n=25; % amount of grid along xy coordinate
x_div=(x_max-x_min)/n; y_div = (y_max-y_min)/n;

grid_x_y_center = [ x_min,y_min,x_min+0.5*x_div,y_min+0.5*y_div];

ii=0;
for ii=1:n
    
    grid_x_y_center(ii+1,1) =  x_min+ ii*x_div;
    grid_x_y_center(ii+1,2) =  y_min+ ii*y_div;
    grid_x_y_center(ii+1,3) =  x_min+ (ii*x_div) + 0.5*x_div  ;
    grid_x_y_center(ii+1,4) =  y_min+ (ii*y_div) + 0.5*y_div  ;
end


%% UE and non UE dividation
whitout_UE= updated_data;

m=100; % UE amount
step = floor(size(updated_data,1)/m) - 20;

l=1;
u=20;

for jj = 1:m
    
    if jj==1
        UE{jj}=  updated_data(l:u,:);
        whitout_UE(l:u,:)=[];
    else
        
        l=l+step;
        u=u+step;
        UE{jj}=  updated_data(l:u,:);
        whitout_UE(l:u,:)=[];
        
    end
    
end

%% Putting into grid

ii=0;
jj=0;
counter=0;
Training_data=[];
for ii=1:n
    
    for jj=1:n
        
        rows=[];
        sample=[];
        counter = counter+1;
        
        if(jj==n)
            
            rows=  find(whitout_UE(:,17)>= grid_x_y_center(ii,1)  & whitout_UE(:,17)<=grid_x_y_center(ii+1,1)  & whitout_UE(:,18)>= grid_x_y_center(jj,2)  & whitout_UE(:,18)<=grid_x_y_center(jj+1,2));
            
        else
            
            rows=  find(whitout_UE(:,17)>= grid_x_y_center(ii,1)  & whitout_UE(:,17)<grid_x_y_center(ii+1,1)  & whitout_UE(:,18)>= grid_x_y_center(jj,2)  & whitout_UE(:,18)<grid_x_y_center(jj+1,2));
            
        end
        
        sample=whitout_UE(rows,:);
        [r c] = size(sample);
        
        if(r~=0)
            sample(1,c+1)=[grid_x_y_center(ii,3)];
            sample(1,c+2)=[grid_x_y_center(jj,4)];
        end
        
        Training_data{counter} = sample;
        
    end
end

%% Training Signature creation


for jj=1:size(Training_data,2)
    
    temp=[];c=[];u=[];I=[];J=[];v=[];K=[];L=[];
    
    if size(Training_data{jj},1)~=0
        
        block=[];
        
        temp= cell2mat(Training_data(jj));
        c = sort(temp(:,3:9),2);
        [u,I,J] = unique(c, 'rows', 'first');
        
        for pp=1:size(u,1)
            
            signature=[];
            same_row_index = ismember(c,u(pp,:),'rows');
            
            checker=sum(same_row_index);
            
            if checker>=2
                same_row_value = temp(same_row_index,:);
                
                rss_avg=[] ;
                for ii=10:18
                    avg = mean( same_row_value(:,ii));
                    rss_avg=horzcat(rss_avg,avg);
                end
                signature = horzcat(u(pp,:),rss_avg);
                block=vertcat(block,signature);
            end
            
            
            
        end
        
        
        
        grid_signature= block;
        grid_signature(1,19)=temp(1,19);
        grid_signature(1,20)=temp(1,20);
        
        
    else
        grid_signature=[];
    end
    
    Training_signature{jj}=grid_signature;
    
    
    
end




% num = length( Training_signature);
% rows = 0;
% for i=1:num
%     [r c] = size( Training_signature{i});
%     rows = rows+r
% end


%% Testing signature creation

for jj=1:size(UE,2)
    
    temp=[];c=[];u=[];I=[];J=[];v=[];ixDupRows=[];dupRowValues=[];K=[];L=[];
    
    if size(UE{jj},1)~=0
        
        temp= cell2mat(UE(jj));
        c = sort(temp(:,3:9),2);
        [u,I,J] = unique(c, 'rows', 'first');
        
        for pp=1:size(u,1)
            
            signature=[];
            same_row_index = ismember(c,u(pp,:),'rows');
            
            checker = sum(same_row_index);
            
            if checker >=2
                same_row_value = temp(same_row_index,:);
                rss_avg=[] ;
                for ii=10:18
                    avg = mean( same_row_value(:,ii));
                    rss_avg=horzcat(rss_avg,avg);
                end
                
                
                signature = horzcat(u(pp,:),rss_avg);
                block=vertcat(block,signature);
                
                
            end
            
            
        end
        
        grid_signature=block;
        
        block=[];
        
    else
        grid_signature=[];
    end
    
    Testing_signature{jj}=grid_signature;
    
end

% Testing_signature{10}(1:2,:)=[]; %Error correction

% num = length(Testing_signature);
% rows = 0;
% for i=1:num
%     [r c] = size(Testing_signature{i});
%     rows = rows+r
% end


%% Test signature matching
counter = 0;
zero(:,1:16)=0;
for ii = 1:length(Testing_signature)
    
    p = Testing_signature{ii};
    
    for jj = 1:size(p,1)
        
        counter = counter +1;
        single_test_signature_id = p(jj,1:7);
        single_test_signature_value = p(jj,1:16);
        from_all_block=[];
        
        for mm=1:length(Training_signature)
            
            if size(Training_signature{mm},1)~=0
                
                block=Training_signature{mm};
                
                find_similar=ismember(block(:,1:7),single_test_signature_id);
                
                checking_index = sum(find_similar,2);
                
                alldata=[];
                
                for kk=1:size(checking_index,1)
                    
                    if checking_index(kk,1) >=7 % variable?
                        member = block(kk,1:16) ;
                        alldata = cat(1, alldata, member);
                    else
                        alldata=[];
                    end
                    
                end
                
                if size(alldata,1)~=0
                    
                    from_all_block=cat(1, from_all_block, alldata);
                    
                end
                
            end
        end
        
        
        MATCHED{counter}= vertcat(single_test_signature_value,zero,from_all_block); %testing
        
        
    end
end

%% error correction
for ii=1:length(MATCHED)
    MATCHED{ii}( all(~MATCHED{ii},2), : ) = [];
end


%% distance calculation
counter = 0;
for jj=1:length(MATCHED)
    
    counter=counter+1;
    
    model = cell2mat(MATCHED(jj));
    
    for ii = 2:size(model,1) %3 because zero has gone
        
        multiple_rss_distance=0;
        
        for kk=1:7
            
            for ll = 1:7
                if model(1,kk) == model(ii,ll)
                    
                    single_rss_distance = (model(1,kk+7)-model(ii,ll+7))^2;
                    
                    multiple_rss_distance= multiple_rss_distance+single_rss_distance;
                    
                else
                    
                end
            end
        end
        
        model(ii,19) = sqrt(multiple_rss_distance);
        
        model(ii,20) = sqrt((model(1,15)-model(ii,15))^2 + (model(1,16)-model(ii,16))^2 );
        
        model_distance{counter} = model;
        
    end
end

%% confidance level

confidance_model=[];
counter =0;
con=0;
for ii = 1:length(model_distance)
    counter= counter+1;
    if size(model_distance{ii},1)~=0
        
        temporary = cell2mat(model_distance(ii));
        
        
        for jj = 2:size(temporary,1)
            
            find_member=ismember(temporary(1,1:7),temporary(jj,1:7));
            con = sum(find_member) / 7; % 7 is the max possible
            temporary(jj,21)=con;
            
        end
        
        confidance_model{counter}=temporary;
    end
    
    
    
end

%% descend according confidance level

counter = 0;
for ii = 1:length(confidance_model)
    counter= counter+1;
    if size(confidance_model{ii},1)~=0  %%I change it
        
        % sort only the first column, return indices of the sort
        [~,sorted_inds] = sort( confidance_model{ii}(:,21),'descend' );
        
        % reorder the rows based on the sorted indices
        sorted = confidance_model{ii}(sorted_inds,:);
        
        confidance_sorted_modle{counter}= sorted;
        
    end
    
end

%% TOP confidance level lower rss distance

counter = 0;
for ii = 1:length(confidance_sorted_modle)
    counter= counter+1;
    if size(confidance_sorted_modle{ii},1)~=0  %%I change it
                
        maxi = max(confidance_sorted_modle{ii}(:,21));
        
        indx = find(confidance_sorted_modle{ii}(:,21)==maxi);
        
        dist_rows= confidance_sorted_modle{ii}(indx,:);
        
        lowest_dist=[];
        line=[];
        
        if size(dist_rows,1)==1
            
            lowest_dist=dist_rows;
            
        elseif size(dist_rows,1)> 1
            
            minrss= min(dist_rows(:,19));
            
            local_index=find(dist_rows(:,19)==minrss);
            
            line = dist_rows(local_index,:);
            
            if size(line,1)==1
                
                lowest_dist = line;
                
            elseif size(line,1)>1
                
                value_mean= mean(line(:,20));
                
                line(1,20)=value_mean;
                
                lowest_dist=line(1,:);
            end
            
        end
        
        low=vertcat(confidance_sorted_modle{ii}(end,:),lowest_dist);
        
        
        low_dist_matched{counter}=low;
    end
    
end




%%
collection_shortest_distance=[];
for ii= 1:length(low_dist_matched)
    
    if size(low_dist_matched{ii},1)>1
        
        shortest_distance = low_dist_matched{ii}(2,20);
        
        collection_shortest_distance = vertcat(collection_shortest_distance,shortest_distance);
        
    end
    
end
Parcentile_68_min_rss = prctile(collection_shortest_distance(:,1),68)
Parcentile_95_min_rss = prctile(collection_shortest_distance(:,1),95)
%% for top minimum 5
%
% logic1=[];
% logic2=[];
% collection_min5_distance=[];
%
% n=2;
%
% for ii= 1:length(sorted_modle)
%
%     if size(sorted_modle{ii},1) >= n+1
%
%         shortest_5 = sorted_modle{ii}(2:n+1,20); %%3
%         result = mean(shortest_5);
%         logic1= vertcat(logic1,result);
%
%     elseif size(sorted_modle{ii},1)>=2 && size(sorted_modle{ii},1)<n+1 %%
%
%         shortest_5 = sorted_modle{ii}(2:end,20);%%
%         result = mean(shortest_5);
%         logic2= vertcat(logic2,result);
%
%     end
% end
%
%
% collection_min5_distance = vertcat(logic1,logic2);
%
% Parcentile_68_avg5_rss = prctile(collection_min5_distance(:,1),68)
% Parcentile_95_avg5_rss = prctile(collection_min5_distance(:,1),95)

%% plotting
%
% cdfplot(collection_shortest_distance);
% hold on
% cdfplot(collection_min5_distance);
% hold off
% axis([0 200 0 1]);
% ax = gca;
% ax.YTick = [0.68 0.98];
% legend('Min','Avg5');


%% load
clc;  clear all;  close all;
load DATA_RFFP_MAY15_LTE_WLAN_7
data = DATA_RFFP_MAY15_LTE_WLAN_7;

%% wrong xy Position delete
% non_zero_xy= find(data(:,17)~= 0 & data(:,18)~= 0 & data(:,7)~= 0  & data(:,8)~= 0 & data(:,9)~= 0);
% updated_data = data(non_zero_xy,:);

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

n= 25; % amount of grid along xy coordinate
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
        %         grid_signature(1,19)=temp(1,19);
        %         grid_signature(1,20)=temp(1,20);
        
        
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


counter=0;
for ii=1:length(Training_signature)
    
    counter = counter+1;
    gcu_ids=[];
    
    if size(Training_signature{ii},1)~=0 && size(Training_signature{ii},1)~=1
        
        gcu= unique(Training_signature{ii}(:,1:7),'sorted');
        gcu_ids=gcu';
        
    elseif size(Training_signature{ii},1)==1

        gcu= Training_signature{ii}(:,1:7);
        gcu_ids=gcu;
         
            
     end
    
    Traing_gcu_ids_only{counter} =  gcu_ids;
end

%% for putting id anc corresponding value

IDS= Traing_gcu_ids_only;

for ii=1:length(Traing_gcu_ids_only)
    
    if size(Traing_gcu_ids_only{ii},1)~=0
        
        n=size(Traing_gcu_ids_only{ii},2);
        line_value=[];
        
        for jj=1:size(Traing_gcu_ids_only{ii},2)
            
            value_block=[];
            
            for kk=1:size(Training_signature{ii},1)
                
                for ll=1:size(Training_signature{ii},2)
                    
                    if Traing_gcu_ids_only{ii}(1,jj)==Training_signature{ii}(kk,ll)
                        value= Training_signature{ii}(kk,ll+7);
                        value_block=vertcat(value_block,value);
                    end
                end
            end
            
            id=Traing_gcu_ids_only{ii}(1,jj);
            mean_val=mean(value_block);
            tt=find(IDS{ii}(1,:)==id);
            IDS{ii}(1,(tt+200))= mean_val; % to put far known place    
            IDS{ii}(1,400)= 0;
            
        end
    end
end


%% for puting xy mean
id_val_gcu=IDS;

for ii=1:length(id_val_gcu)
    
    if size(id_val_gcu{ii},1)~=0
        
        id_val_gcu{ii}(1,end+1)= mean(Training_signature{ii}(:,15));
        id_val_gcu{ii}(1,end+1)= mean(Training_signature{ii}(:,16));
    end
    
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

num = length(Testing_signature);
rows = 0;
for i=1:num
    [r c] = size(Testing_signature{i});
    rows = rows+r
end

%% Test signature matching
counter = 0;
for ii = 1:length(Testing_signature)
    
    p = Testing_signature{ii};
    
    for jj = 1:size(p,1)
        
        counter = counter +1;
        single_test_signature_id = p(jj,1:7);
        single_test_signature_value = p(jj,1:16);
        single_test_signature_value(jj,402)=0;
        
        from_all_block=[];
        
        
        
        for mm=1:length(Traing_gcu_ids_only)
            
            if size(id_val_gcu{mm},1)~=0
                
                block=id_val_gcu{mm};
                
                find_similar=ismember(block(1,:),single_test_signature_id);
                
                checking_index = sum(find_similar,2);
                
                sig=[];
                
                if checking_index >=7  % variable?
                    
                    sig=vertcat(sig,block(1,:));

                    
                else
                    sig=[];
                end
                
                if size(sig,1)~=0
                    
                    from_all_block=vertcat(from_all_block, sig);
                    
                 %   from_all_block=cat(1, from_all_block, sig);
                end
                
            end
        end
        
        
        MATCHED{counter}=vertcat(single_test_signature_value,from_all_block);
        
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
    
    for ii = 2:size(model,1)
        
        multiple_rss_distance=0;
        
        for kk=1:7
            
            for ll = 1:200
                if model(1,kk) == model(ii,ll)
                    
                    single_rss_distance = (model(1,kk+7)-model(ii,ll+200))^2;
                    
                    multiple_rss_distance= multiple_rss_distance+single_rss_distance;
                    
                else
                    
                end
            end
        end
        
        model(ii,403) = sqrt(multiple_rss_distance);

        model(ii,404) = sqrt((model(1,15)-model(ii,401))^2 + (model(1,16)-model(ii,402))^2 );
        
        model_distance{counter} = model;
        
    end
end


%% Ascending

counter = 0;
for ii = 1:length(model_distance)
    counter= counter+1;
    if size(model_distance{ii},1)~=0  %%I change it
        
        % sort only the first column, return indices of the sort
        [~,sorted_inds] = sort( model_distance{ii}(:,403) );
        
        % reorder the rows based on the sorted indices
        sorted = model_distance{ii}(sorted_inds,:);
        
        sorted_modle{counter}= sorted;
        
    end
    
end


%% result collection for shortest one

collection_shortest_distance=[];
for ii= 1:length(sorted_modle)
    
    if size(sorted_modle{ii},1)>=2
        
        shortest_distance = sorted_modle{ii}(2,404);
        
        collection_shortest_distance = vertcat(collection_shortest_distance,shortest_distance);
        
    end
    
end
Parcentile_68_min_rss = prctile(collection_shortest_distance(:,1),68)
Parcentile_95_min_rss = prctile(collection_shortest_distance(:,1),95)
% save PE_min_rss_dist collection_shortest_distance;
% mean(collection_shortest_distance)


%% for top minimum 5

logic1=[];
logic2=[];
collection_min5_distance=[];
for ii= 1:length(sorted_modle)
    
    if size(sorted_modle{ii},1) >= 6
        
        shortest_5 = sorted_modle{ii}(2:6,404);
        result = mean(shortest_5);
        logic1= vertcat(logic1,result);
        
    elseif size(sorted_modle{ii},1)>=2 && size(sorted_modle{ii},1)<6
        
        shortest_5 = sorted_modle{ii}(2:end,404);
        result = mean(shortest_5);
        logic2= vertcat(logic2,result);
        
    end
end

collection_min5_distance = vertcat(logic1,logic2);
Parcentile_68_avg5_rss = prctile(collection_min5_distance(:,1),68)
Parcentile_95_avg5_rss = prctile(collection_min5_distance(:,1),95)

%% test
% cdfplot(collection_shortest_distance);
% hold on
% cdfplot(collection_min5_distance);
% hold off
% axis([0 200 0 1]);
% ax = gca;
% ax.YTick = [0.68 0.98];
% legend('Min','Avg5');


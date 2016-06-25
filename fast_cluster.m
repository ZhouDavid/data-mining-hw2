%Algorithm:fast_cluster   
%Author:Zhou Jianyu  2013011326 Tsinghua University

%step1:load data
clear;
raw_data=load('example_distances.dat');
record_num = size(raw_data,1);%数据的行数
max_id1 = max(raw_data(:,1));
max_id2 = max(raw_data(:,2));
max_id = max([max_id1,max_id2]);
dist = zeros(max_id);
dc = 0.033;   %设置dc
for i = 1:record_num   %遍历数据，给矩阵赋值，dist[i][j]表示id为i和j两点的距离
    dist(raw_data(i,1),raw_data(i,2)) = raw_data(i,3);
    dist(raw_data(i,2),raw_data(i,1)) = raw_data(i,3);
end

%计算每个点的p和tau
p = zeros(max_id,1);
tau = zeros(max_id,1);
for i=1:max_id-1
    for j =i+1:max_id
       if(dist(i,j)<=dc&&dist(i,j)~=0)
           p(i)=p(i)+1;
           p(j)=p(j)+1;
       end
    end
end
min = 1000;
for i = 1:max_id
    for j=1:max_id
        if(p(j)>p(i))
            if(min>dist(i,j))
                if(dist(i,j)==0)
                    fprintf('%d,%d',i,j);
                    a=input('hehe:');
                end
                min = dist(i,j);
            end
        end
    end
    tau(i)=min;
end
subplot(2,1,1)
plot(p(:),tau(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')


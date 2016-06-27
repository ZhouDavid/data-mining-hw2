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

d = sort(raw_data(:,3));
percent = 0.03;
dc = d(record_num*percent);   %设置dc
for i = 1:record_num   %遍历数据，给矩阵赋值，dist[i][j]表示id为i和j两点的距离
    dist(raw_data(i,1),raw_data(i,2)) = raw_data(i,3);
    dist(raw_data(i,2),raw_data(i,1)) = raw_data(i,3);
end

%计算每个点的p和tau
p = zeros(max_id,1);
tau = zeros(max_id,1);
for i=1:max_id-1
    for j =i+1:max_id
       p(i)=p(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
       p(j)=p(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
    end
end
% min = 1000;
% for i = 1:max_id-1
%     for j=i+1:max_id
%         if(p(j)>p(i))
%             if(min>dist(i,j))
%                 if(dist(i,j)==0)
%                     fprintf('%d,%d',i,j);
%                     a=input('hehe:');
%                 end
%                 min = dist(i,j);
%             end
%         end
%     end
%     tau(i)=min;
% end



[psort,pos]=sort(p,'descend');
tau(pos(1))=-1.;
nneigh(pos(1))=0;

for i=2:max_id
   tau(pos(i))=1000;  
   for j=1:i-1
     if(dist(pos(i),pos(j))<tau(pos(i)))
        tau(pos(i))=dist(pos(i),pos(j));
        nneigh(pos(i))=pos(j);
     end
   end
end
tau(pos(1)) = max(tau);
subplot(2,1,1)
plot(p(:),tau(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')

y = zeros(max_id,1);
y = p.*tau;
%设置y的阈值
ysort = sort(y,'descend');
ythreshold = ysort(6);
cluster_num = 0;

tausort=sort(tau,'descend');
tauthreshold = tausort(percent*max_id);
%找到聚类中心
for i=1:max_id
  id2cluster(i)=-1;
end
for i=1:max_id
    if(y(i)>ythreshold)
        cluster_num=cluster_num+1;
        cluster_center(cluster_num) = i;  %记录第i个聚类中心的id
        id2cluster(i) = cluster_num;
    end
end


for i =1:max_id
    if(id2cluster(pos(i))==-1)
        id2cluster(pos(i)) = id2cluster(nneigh(pos(i)));%一个点所处的聚类中心编号与比该点density高且距离最近的点的编号相同
    end
end
cmap=colormap;

for i=1:cluster_num
    ii=int8((i*64.)/(cluster_num*1.));
    subplot(2,1,1)
    hold on
    plot(p(cluster_center(i)),tau(cluster_center(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ii,:),'MarkerEdgeColor',cmap(ii,:));
end

%识别噪声点
halo = id2cluster;
pb = zeros(cluster_num,1);
for i=1:max_id-1
    for j = i+1:max_id
        if((id2cluster(i)~=id2cluster(j))&&(dist(i,j)<dc))
            border_density = (p(i)+p(j))/2;
            if(p(i)>pb(id2cluster(i)))
                pb(id2cluster(i)) = p(i);
            end
            if(p(j)>pb(id2cluster(j)))
                pb(id2cluster(j)) = p(j);
            end
        end
    end
end

for i = 1:max_id
    if(p(i)<pb(id2cluster(i)))
        halo(i)=0;
    end
end

subplot(2,1,2);
Y = mdscale(dist, 2, 'criterion','metricstress');
plot(Y(:,1),Y(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
for i=1:max_id
    cluster_id = halo(i);
    color_id = int8((cluster_id*64.)/(cluster_num*1.));
    if(color_id==0) 
        continue;
    end
    hold on
    plot(Y(i,1),Y(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(color_id,:),'MarkerEdgeColor',cmap(color_id,:));
end
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')

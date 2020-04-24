clc,clear,close all
warning off
format longG
%% SOA 参数
maxiter = 20;  % 迭代次数
sizepop = 10;  % 种群数量
Umax=0.9500;   % 最大隶属度值
Umin=0.0111;   % 最小隶属度值
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
fai1 = 0.5;
fai2 = 0.5;
w1 = 0.5;
% w2 = 0.5;
w2max = 0.7;
w2min = 0.2;
%% 初始化种群
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
    fitness(i) = fun([x1,x2]);
end
%% 记录一组最优值
[bestfitness,bestindex]=min(fitness); % 索引最小值
zbest=pop(bestindex,:);   % 全局最佳
gbest=pop;                % 个体最佳
fitnessgbest=fitness;     % 个体最佳适应度值
fitnesszbest=bestfitness; % 全局最佳适应度值
Frontpop = pop;
%% 迭代寻优
for i=1:maxiter
    for j=1:sizepop
       dego = gbest(j,:)-pop(j,:);
       dalt = zbest-pop(j,:);
       dpro = Frontpop(j,:) - pop(j,:);
            Frontpop = pop;  % 记录上一时刻的最优种群
       dij = sign( w1*dpro+fai1*dego+fai2*dalt );   % 确定方向
       % 确定步长
       [maxfitness,index] = sort(fitness,'descend');
       u = Umax - index(j)/sizepop *(Umax-Umin);
       u = u.*rand(1,2);          % 2个变量
       T = sqrt(-log(u) );
       xmin = pop(index(1),:);    % 最差的个体
       xmax = pop(index(end),:);  % 最优的个体
      
       w2 = w2max - j/sizepop*(w2max-w2min);
       delta = w2*abs( xmin-xmax );
      
       alpha = delta.*T;
      
       % 步长进行限制
       for k=1:length(alpha)
           if alpha(k)>0.6
               alpha(k)=0.6;
           elseif alpha(k)<-0.6
               alpha(k) = -0.6;
           end
       end
      
       % 位置更新
       pop(j,:) = pop(j,:) + dij.*alpha;        

        % x1  越界限制
        if pop(j,1)>popmax1
            pop(j,1)=popmax1;
        end
        if pop(j,1)<popmin1
            pop(j,1)=popmin1;
        end
        % x2  越界限制
        if pop(j,2)>popmax2
            pop(j,2)=popmax2;
        end
        if pop(j,2)<popmin2
            pop(j,2)=popmin2;
        end
        
        % 适应度更新
        fitness(j) = fun(pop(j,:));
        
        % 比较  个体间比较
        if fitness(j)<fitnessgbest(j)
            fitnessgbest(j) = fitness(j);
            gbest(j,:) = pop(j,:);
        end
        if fitness(j)<bestfitness
            bestfitness = fitness(j);
            zbest =  pop(j,:);
        end
        
    end
    fitness_iter(i) = bestfitness;
   
end
disp('最优解')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight
grid on
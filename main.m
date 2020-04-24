clc,clear,close all
warning off
format longG
%% SOA ����
maxiter = 20;  % ��������
sizepop = 10;  % ��Ⱥ����
Umax=0.9500;   % ���������ֵ
Umin=0.0111;   % ��С������ֵ
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
fai1 = 0.5;
fai2 = 0.5;
w1 = 0.5;
% w2 = 0.5;
w2max = 0.7;
w2min = 0.2;
%% ��ʼ����Ⱥ
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
    fitness(i) = fun([x1,x2]);
end
%% ��¼һ������ֵ
[bestfitness,bestindex]=min(fitness); % ������Сֵ
zbest=pop(bestindex,:);   % ȫ�����
gbest=pop;                % �������
fitnessgbest=fitness;     % ���������Ӧ��ֵ
fitnesszbest=bestfitness; % ȫ�������Ӧ��ֵ
Frontpop = pop;
%% ����Ѱ��
for i=1:maxiter
    for j=1:sizepop
       dego = gbest(j,:)-pop(j,:);
       dalt = zbest-pop(j,:);
       dpro = Frontpop(j,:) - pop(j,:);
            Frontpop = pop;  % ��¼��һʱ�̵�������Ⱥ
       dij = sign( w1*dpro+fai1*dego+fai2*dalt );   % ȷ������
       % ȷ������
       [maxfitness,index] = sort(fitness,'descend');
       u = Umax - index(j)/sizepop *(Umax-Umin);
       u = u.*rand(1,2);          % 2������
       T = sqrt(-log(u) );
       xmin = pop(index(1),:);    % ���ĸ���
       xmax = pop(index(end),:);  % ���ŵĸ���
      
       w2 = w2max - j/sizepop*(w2max-w2min);
       delta = w2*abs( xmin-xmax );
      
       alpha = delta.*T;
      
       % ������������
       for k=1:length(alpha)
           if alpha(k)>0.6
               alpha(k)=0.6;
           elseif alpha(k)<-0.6
               alpha(k) = -0.6;
           end
       end
      
       % λ�ø���
       pop(j,:) = pop(j,:) + dij.*alpha;        

        % x1  Խ������
        if pop(j,1)>popmax1
            pop(j,1)=popmax1;
        end
        if pop(j,1)<popmin1
            pop(j,1)=popmin1;
        end
        % x2  Խ������
        if pop(j,2)>popmax2
            pop(j,2)=popmax2;
        end
        if pop(j,2)<popmin2
            pop(j,2)=popmin2;
        end
        
        % ��Ӧ�ȸ���
        fitness(j) = fun(pop(j,:));
        
        % �Ƚ�  �����Ƚ�
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
disp('���Ž�')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight
grid on
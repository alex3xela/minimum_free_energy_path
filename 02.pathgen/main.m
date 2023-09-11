
clc;
clear all;
close all;

%%%%%%%%%%%% parameters %%%%%%%%%%%%%%
% number of clusters, heuristically num of samples/200
NC=20;
% resample the string 
newNC = 30;
% smoothing parameter. The higher the smoother
s=50;
% selects here the indices of the start and the target frame
boundaries_id=[1,998];
subsample = 1;
% dist matrix file name
dist = 'mat.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state',0);
randn('state',0);

dst_mtx = load(dist);
dst_mtx = dst_mtx.^2;
N = size(dst_mtx,1);
for i=1:N
    dst_mtx(i,i)=0;
end

if (subsample>1)
    rr = unique(randi(N,[1,round(N/subsample)]));
    dst_mtx = dst_mtx(rr,rr);
    N = size(dst_mtx,1);
    boundaries_id(2)=N;
end

if (boundaries_id(1)>N || boundaries_id(2)>N)
    error('Start and stop indices must be < N');
end

%% kernel 

%sigma2 = 2*mean(mean(dst_mtx));
%sigma2 = 100;
sigma2 = 10*max(max(dst_mtx))/(sqrt(2*N));
k_mtx = exp(-dst_mtx/sigma2);

%% init path and minimize functional

%%%%
%run the regularized path in kernel space algorithm
%%%%
init_medoids_id = rpksInitMedoids('kpp',k_mtx,boundaries_id,NC);
init_medoids_id = rpksInitPath('greedy',k_mtx,init_medoids_id);
[medoids_id,labels,cost,card]=rpks(k_mtx, init_medoids_id, s);


%%%%
%plot cost vs iterations
%%%%
h_cost=figure;
plot(cost);

%% projection in R3
display(['Projecting in R3..']);
% project data into R3, assume R3 is good for representing data
D = zeros(N,N);
for i=1:N
    for j=i:N
        D(i,j) = k_mtx(i,i)+k_mtx(j,j)-2*k_mtx(i,j);
        D(j,i) = D(i,j);
    end
end
M = zeros(N,N);
for i=1:N
    for j=i:N
        M(i,j) = D(1,j)+D(i,1)-D(i,j);
        M(j,i) = M(i,j);
    end
end
M = 0.5*M;
[U S V] = svd(M);
X = U*sqrt(S);
D = 3;

%% plot and save section %%

figure;
medoids = X(medoids_id,:);
plot3(medoids(:,1),medoids(:,2),medoids(:,3),'r-x','MarkerSize',20);
hold on;
plot3(X(boundaries_id(1),1),X(boundaries_id(1),2),X(boundaries_id(1),3),'m-*','MarkerSize',20);
plot3(X(boundaries_id(2),1),X(boundaries_id(2),2),X(boundaries_id(2),3),'k-*','MarkerSize',20);
axis equal;

fp = fopen('medoids_indices','w');
for i=1:length(medoids_id)
    fprintf(fp,'%d\n',medoids_id(i));
end
NC = length(card)+2;
fclose(fp);
[Xnew,combo] = overSample(k_mtx,NC,newNC,medoids_id,X);
scatter3(Xnew(:,1),Xnew(:,2),Xnew(:,3),50*ones(size(Xnew,1),1),50*ones(size(Xnew,1),1),'go'), hold on;
ss = scatter3(X(:,1),X(:,2),X(:,3),20*ones(size(X,1),1),0.5*ones(size(X,1),1) , '.');
alpha(ss,0.01);
grid on;
axis equal

i1 = boundaries_id(1);
i2 = boundaries_id(2);
combo = [[1 0 i1 1]; combo; [1 0 i2 1]];
fileID = fopen('medoids_oversampled.txt','w');
for i=1:size(combo,1)
    fprintf(fileID,'%f %f %.0f %.0f\n',combo(i,1),combo(i,2),combo(i,3)-1,combo(i,4)-1);
end
fclose(fileID);

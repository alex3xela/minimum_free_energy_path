clc;
clear all;
close all;

%%%%%%%%%%%% parameters %%%%%%%%%%%%%%
% number of clusters, heuristically num of samples/200
NC= 25;
% resample the string to get equispaced medoids
newNC = 50;
% smoothing parameter. The higher the smoother
s=400;
% selects here the indices of the start and the target frame
boundaries_id=[1,379];
% input matrix (concatenated coordinates of the ligand)
inputMatrix = './ligandMatrix.txt';
showText = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state',0);
randn('state',0);

X = load(inputMatrix);
N = size(X,1);
NA = size(X,2)/3;

% build MSD dist matrix for later projection
distMatrix = zeros(N,N);
display('Computing MSD matrix');
for i=1:N
    for j=i:N
        dist = 1.0/NA*sum((X(i,:)-X(j,:)).^2);
        distMatrix(i,j)=dist;
        distMatrix(j,i)=dist;
    end
end

if (boundaries_id(1)>N || boundaries_id(2)>N)
    error('Start and stop indices must be < N');
end

k_mtx = X*X';

init_medoids_id = rpksInitMedoids('kpp',k_mtx,boundaries_id,NC);
init_medoids_id = rpksInitPath('projection',k_mtx,init_medoids_id);
init_centroids = X(init_medoids_id,:);

[medoids,labels,cost]=rpls(X, init_centroids, s);

% search the nearest samples 
medoids_id = [];
for j=1:size(medoids,1)
    c = medoids(j,:);
    minDist = 1e20;
    winner = -1;
    for i=1:size(X,1)
        x = X(i,:);
        dist = sum((x-c).^2);
        if (dist<minDist)
            minDist = dist;
            winner = i;
        end
    end
    medoids_id = [medoids_id winner];
end
%% projection in R3 
display(['Projecting in R3..']);

D = distMatrix;
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
X = X(:,1:3);
D = 3;

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
NC = length(medoids_id);
fclose(fp);
[Xnew,combo] = overSample(k_mtx,NC,newNC,medoids_id,X);
hold on;
scatter3(Xnew(:,1),Xnew(:,2),Xnew(:,3),50*ones(size(Xnew,1),1),50*ones(size(Xnew,1),1),'go'), hold on;
prova1=X(:,1);
pr=prova1(:);
pr(1:length(pr))=0;
pr(1:145)=2;
ss = scatter3(X(:,1),X(:,2),X(:,3),20*ones(size(X,1),1),pr, '.');



ss = scatter3(X(:,1),X(:,2),X(:,3),20*ones(size(X,1),1),0.5*ones(size(X,1),1), '.');
if (showText==1)
    text(X(:,1),X(:,2),X(:,3),cd);
end
%ss = scatter3(X(:,1),X(:,2),X(:,3),20*ones(size(X,1),1),0.5*ones(size(X,1),1) , '.');
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

%% search the nearest physical points to the equidistant points

winners = [];
for i=1:size(combo,1)
    winner = -1;
    minDist = 1e20;
    convex = combo(i,1)*X(combo(i,3),:)+combo(i,2)*X(combo(i,4),:);
    for j=1:size(X,1)
        dist = sum((convex-X(j,:)).^2);
        if (dist<minDist)
            minDist = dist;
            winner = j;
        end       
    end
    winners = [winners winner];
end

fileID = fopen('phyisical_medoids.txt','w');
for i=1:size(combo,1)
    fprintf(fileID,'%d\n',winners(i)-1);
end
fclose(fileID);

hold on;
%ss = scatter3(X(winners,1),X(winners,2),X(winners,3),20*ones(size(X(winners,:),1),1),0.5*ones(size(X(winners,:),1),1), 'k*');
plot3(X(winners,1),X(winners,2),X(winners,3),'k-*','MarkerSize',10);

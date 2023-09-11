

function [Xnew,combo,distStations] = overSample(K,Nc,numProt,medoids,X)

% compute the new string given the desired number of
% proptotypes. returns the points in the embedded space (Xnew)
% the matrix of pairwise optimal weighting to get the oversample clusters
% (combo); the first and second columns are the weight of respectively the third
% and fourth column where there are the meodoids indices they refer lastly it gives the vector of distances in the kernel spaces among the
% consecutive new oversampled clusters. 
% To get the new clusters 

% update D/cluDist auxiliary vectors
D = zeros(Nc,1);
cluDist = zeros(Nc,1);
ss = zeros(numProt,1);
Xnew = [];
D(1)=0;
cluDist(1)=0;
ss(1)=0;
N = size(K,1);
combo = [];
combo2 = [];

for hh=2:Nc
    D(hh) = 0;
    for q=2:hh
       % dist q-1, q
       iqm1 = medoids(q-1);
       iq = medoids(q);
       aa = K(iqm1,iqm1);      
       bb = K(iq,iq);
       cc = K(iq,iqm1);
       dist2 = aa+bb-2*cc;
       dist = sqrt(dist2);
       cluDist(q)=dist;
       D(hh)=D(hh)+dist;
    end
end

% delta dist along the string
delta = D(Nc)/(numProt-1);
for hh=2:numProt
    ss(hh)=(hh-1)*delta;
end

for j=2:numProt-1        
   % j = p in the paper
   % post = q in the paper    
   pre = 1;
   post = 2;
   test=0;             

   % select the vector on which to place the new cluster
   % that is we are moving along the path until we can reach
   % the right distance on the path             
   % can be optimized doing a map once and for all
   % after this loop we have localized the two points
   % in between the new cluster must be placed
   while(test==0)
      if (ss(j)<=D(post) && ss(j)>D(pre))
          test=1;
      else
          post = post+1;
          pre = pre+1;
      end                     
   end

   alpha = (ss(j)-D(pre));
   gamma = alpha/cluDist(post);
  
   linearCombo = zeros(N,1);
   linearCombo(medoids(pre)) = (1-gamma);   
   linearCombo(medoids(post)) = (gamma); 
   combo = [combo; 1-gamma gamma medoids(pre) medoids(post)];
   combo2 = [combo2; linearCombo(:)'];
   % also we want to have a look in the prjected space how
   % the convex combination looks like
   Xnew = [Xnew; linearCombo(:)'*X];   
end


% test the quality of linear combinations by checking inter stations distances
% in the original kernel space
val = [];
for i=1:size(Xnew,1)-1
   mixed = 0;
   mixed = combo2(i,:)*K*combo2(i+1,:)';
   mixed = -2*mixed;  
   val(i) = combo2(i,:)*K*combo2(i,:)'+combo2(i+1,:)*K*combo2(i+1,:)'+mixed;
end

distStations = val'
avgStat = mean(distStations)
stdDev = std(distStations)
function Random_Forrest(Tree_num, Depth, X)
sizeX=size(X);
N=sizeX(1);
P=sizeX(2)-1;
DATA=X(:,2:end);
CLASS=X(:,1);


for i=1:Tree_num
    Train_DT(P,Depth,CLASS,DATA);
end


end
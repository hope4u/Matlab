%%%%%%%%%%%%%%%%%%%%%
%   removed lines   %
%%%%%%%%%%%%%%%%%%%%%


%% generate Data
numOfSymb = N;
bitsPerSymb = 2;
numOfBits = 2*N;

data = randi([0,1],numOfBits,1)*2-1;
x = reshape(data,[],2) * [1;1i] * 1/sqrt(2);


%% Titles
        fprintf(1,'\n%s:\n',Type{i});
        fprintf(1,'%s\n',repmat('-',1,length(Type{i})+1));
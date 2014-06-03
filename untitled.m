[a,b,c] = size(Error)
%     for cc=1:c
%         if sum(sum(Error))
%             cc
%         end
%     end
error = [];
for cc=1:c    
    if R_sum(1,1,6,cc)<8
        error(end+1)=cc
    elseif R_sum(2,1,6,cc)<20
        error(end+1)=cc
    elseif R_sum(3,1,6,cc)<40
        error(end+1)=cc
    elseif R_sum(4,1,6,cc)<60
        error(end+1)=cc
    end
end
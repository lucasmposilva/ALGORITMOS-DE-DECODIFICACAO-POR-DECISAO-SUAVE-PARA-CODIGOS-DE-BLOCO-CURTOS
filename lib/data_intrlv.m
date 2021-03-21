function data_intrlv(ifile,ofile,N)

load(ifile,"data");

nrow = size(data,1);
ncol = size(data,2);

out                 = zeros(nrow,N.*ncol);
out(1:nrow,1:ncol)  = data;

for k = 1:N
    out([1:nrow],[(k.*ncol)+1:((k+1).*ncol)]) = intrlv(data,randperm(nrow));
end
clear data; data = out; clear out;
save(ofile,"data");
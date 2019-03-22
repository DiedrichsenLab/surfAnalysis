function names=surf_getGiftiColumnNames(G)
% function names=getGiftiColumnNames(this)
% returns the column names from a functional gifti file
this = builtin('struct',G);
N=length(this.data); 
for n=1:N 
    for i=1:length(this.data{n}.metadata)
        if (strcmp(this.data{n}.metadata(i).name,'Name'))
            names{n}=this.data{n}.metadata(i).value;
            break;
        end;
    end; 
end; 
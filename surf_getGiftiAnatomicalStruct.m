function anatomicalStruct=surf_getGiftiAnatomicalStruct(G)
% function names=getGiftiAnatomicalStruct(this)
% returns the Primary anatomical structure for a gifti object
this = builtin('struct',G);
N=length(this.metadata); 
anatomicalStruct=[];
for i=1:N 
    if (strcmp(this.metadata(i).name,'AnatomicalStructurePrimary'))
        anatomicalStruct=this.metadata(i).value;
        return; 
    end;
end; 
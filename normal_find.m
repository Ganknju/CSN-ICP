%normal_find()作用是将坐标点与对应的法向量对应

function [outpoint_normal] = normal_find(outpoint,normal)
outpoint_normal = zeros(size(outpoint,1),6);
for i = 1:size(outpoint,1)                               %遍历行数           
    normal_x =  normal(:,1);                            %取第一列中元素
    xindex = find(normal_x==outpoint(i,1)) ; %找到x轴对应关键点，find（）返回位置
    normal_y = normal(xindex,2);
    yindex = find(normal_y==outpoint(i,2)) ;
    normal_z = normal(xindex(yindex),3);
    zindex = find(normal_z==outpoint(i,3)) ;
    index = xindex(yindex(zindex));               %同时满足xyz的关键点
    outpoint_normal(i,:) = normal(index,:);
end

end
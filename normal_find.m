%normal_find()�����ǽ���������Ӧ�ķ�������Ӧ

function [outpoint_normal] = normal_find(outpoint,normal)
outpoint_normal = zeros(size(outpoint,1),6);
for i = 1:size(outpoint,1)                               %��������           
    normal_x =  normal(:,1);                            %ȡ��һ����Ԫ��
    xindex = find(normal_x==outpoint(i,1)) ; %�ҵ�x���Ӧ�ؼ��㣬find��������λ��
    normal_y = normal(xindex,2);
    yindex = find(normal_y==outpoint(i,2)) ;
    normal_z = normal(xindex(yindex),3);
    zindex = find(normal_z==outpoint(i,3)) ;
    index = xindex(yindex(zindex));               %ͬʱ����xyz�Ĺؼ���
    outpoint_normal(i,:) = normal(index,:);
end

end
%% start
load('Tables.mat');
QR = imread('input/Hello_phone15.jpg');
gray_QR =rgb2gray(QR);
total_gray = sum(sum(gray_QR));
total_gray_unit = size(find(gray_QR),1);
mean_total_gray = total_gray/total_gray_unit;
level = graythresh(QR);
QR_binary_invert = im2bw(255-gray_QR,1-level);
QR_binary = im2bw(gray_QR,level);

Harris_corner_detect_num = round(sqrt(size(gray_QR,1)*size(gray_QR,2)));

temp_QR_binary_invert = QR_binary_invert;
CC = bwconncomp(temp_QR_binary_invert);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
temp_QR_binary_invert(CC.PixelIdxList{idx}) = 0;


L = bwlabel(temp_QR_binary_invert);
R = regionprops(L,'Area','BoundingBox','PixelList');
NR = numel(R);

cell_R = struct2cell(R);
mat_R_Area = cell2mat(cell_R(1,:));
[~,mat_R_Area_index] = sort(mat_R_Area,'descend');

ALL_Point_outer = zeros(3,8);
ALL_Point_outer_format2 = zeros(4,6);
ALL_Point_outer_index = zeros(3,5);

for i=1:NR
BBmax = R(i).BoundingBox;
DIAG1 = sum(R(i).PixelList,2);
DIAG2 = diff(R(i).PixelList,[],2);

[~,dUL] = min(DIAG1);    [~,dDR] = max(DIAG1);
[~,dDL] = min(DIAG2);    [m,dUR] = max(DIAG2);

pts = R(i).PixelList([dUL dDL dDR dUR dUL],:);
ALL_Point_outer(i,:) = [pts(1,:),pts(2,:),pts(3,:),pts(4,:)];
ALL_Point_outer_format2(:,2*i-1:2*i) = pts(1:4,:);
ALL_Point_outer_index(i,:) = [i dUL dDL dDR dUR];
end
detect_point_set = [vec2mat(ALL_Point_outer_format2(1,:),2);vec2mat(ALL_Point_outer_format2(2,:),2);...
    vec2mat(ALL_Point_outer_format2(3,:),2);vec2mat(ALL_Point_outer_format2(4,:),2)];
find_num = 0; 
ALL_Point_outer_format2(5,:)= ALL_Point_outer_format2(1,:);
Point_outer_format2 = zeros(4,6);

outer_Area = zeros(1,3);
for i=1:size(mat_R_Area_index,2)
    edge_shape_index = mat_R_Area_index(i);
    edge_set = [ALL_Point_outer_format2(:,edge_shape_index*2-1),ALL_Point_outer_format2(:,edge_shape_index*2)];
    [in,on] = inpolygon(detect_point_set(:,1),detect_point_set(:,2),...
        edge_set(:,1),edge_set(:,2));
    if sum(in)-sum(on) == 4      
        fon = find(on);
        in(fon) = 0;
        in_points = detect_point_set(in,:);
        edge_aera = polyarea(edge_set(1:4,1)',edge_set(1:4,2)');
        in_aera = polyarea(in_points(1:4,1)',in_points(1:4,2)');
        area_ratio = edge_aera/in_aera;
        if area_ratio >4 && area_ratio<6
            find_num = find_num+1;
            hold on; 
            plot(edge_set(:,1),edge_set(:,2),'y','linewidth',3);
            Point_outer_format2(:,find_num*2-1:find_num*2) = edge_set(1:4,:);
            outer_Area(find_num) = mat_R_Area(edge_shape_index);
            if find_num == 3
                break
            end
        end
    end
end 

Center_outer_block = mean(Point_outer_format2);
Center_outer_block = reshape(Center_outer_block,[2,3])';
triangle_length = pdist(Center_outer_block);
Center_outer_triangle = mean(Center_outer_block);
max_outer_Area = max(outer_Area);
min_outer_Area = min(outer_Area);
for i=1:NR
    if  R(i).Area >max_outer_Area || R(i).Area < min_outer_Area
        if pdist([R(i).BoundingBox(:,1:2);Center_outer_triangle])> mean(triangle_length)
            temp_QR_binary_invert(R(i).PixelList(:,2),R(i).PixelList(:,1))=0;
            %disp(i);
            %hold on;
            %plot(R(i).PixelList(:,1),R(i).PixelList(:,2),'y*');
        end
    end
end


%% order 3 block
% some value needed 
Edge_temp_QR_invert_b = edge(temp_QR_binary_invert);
total_count_Edge_QR_invert_b = sum(sum(Edge_temp_QR_invert_b));
P2X_E_QR_invert_b = sum(Edge_temp_QR_invert_b);
P2Y_E_QR_invert_b = sum(Edge_temp_QR_invert_b,2);
PrL_E_Y = find(P2Y_E_QR_invert_b,1,'last')- find(P2Y_E_QR_invert_b,1,'first');
PrL_E_X = find(P2X_E_QR_invert_b,1,'last')- find(P2X_E_QR_invert_b,1,'first');
Thread2X_find_boundary = total_count_Edge_QR_invert_b/PrL_E_Y;
Thread2Y_find_boundary = total_count_Edge_QR_invert_b/PrL_E_X;
F_bounadry_X = [find(P2X_E_QR_invert_b>Thread2X_find_boundary,1,'first'),find(P2X_E_QR_invert_b>Thread2X_find_boundary,1,'last')];
F_bounadry_Y = [find(P2Y_E_QR_invert_b>Thread2Y_find_boundary,1,'first'),find(P2Y_E_QR_invert_b>Thread2Y_find_boundary,1,'last')];
F_boundry_Center = [mean(F_bounadry_X),mean(F_bounadry_Y)];
hold on;
plot(F_boundry_Center(1),F_boundry_Center(2),'r*');
%plot(F_bounadry_X(1),10,'y*');
%plot(F_bounadry_X(2),10,'y*');
%plot(10,F_bounadry_Y(1),'y*');
%plot(10,F_bounadry_Y(2),'y*');

% change row to colum
%Center_outer_block = reshape(Center_outer_block,[2,3])';
dis = [pdist([mean(Center_outer_block(2:3,:));F_boundry_Center]);...
    pdist([mean([Center_outer_block(1,:);Center_outer_block(3,:)]);F_boundry_Center]);...
    pdist([mean(Center_outer_block(1:2,:));F_boundry_Center])];

CenterXY_index_dis_outer_block = [Center_outer_block,[1;2;3],dis];
[~,index_Temp] = sort(CenterXY_index_dis_outer_block(:,4));
Sorted_CenterXY_index_dis_outer_block = CenterXY_index_dis_outer_block(index_Temp,:);
center_longedge = sum(Sorted_CenterXY_index_dis_outer_block(2:3,1:2))'/2;
largest_diff_center = Sorted_CenterXY_index_dis_outer_block(1,1:2)'-center_longedge;
largest_diff_center = sign(largest_diff_center);
if largest_diff_center(1)<0 && largest_diff_center(2)<0
    QR_Orient = 1;
elseif largest_diff_center(1)>0 && largest_diff_center(2)<0
    QR_Orient = 2;
elseif largest_diff_center(1)<0 && largest_diff_center(2)>0
    QR_Orient = 4;
else 
    QR_Orient = 3;
end
position = sign( (Sorted_CenterXY_index_dis_outer_block(1,1)-center_longedge(1))*(Sorted_CenterXY_index_dis_outer_block(2,2)-center_longedge(2))...
    - (Sorted_CenterXY_index_dis_outer_block(1,2)-center_longedge(2))*(Sorted_CenterXY_index_dis_outer_block(2,1)-center_longedge(1)));
if position == 1
    Sorted_CenterXY_index_dis_outer_block([2 3],:)=Sorted_CenterXY_index_dis_outer_block([3 2],:);
end

ordered_index = 2*Sorted_CenterXY_index_dis_outer_block(:,3);
ordered_index = [ordered_index(1)-1,ordered_index(1),ordered_index(2)-1,ordered_index(2),ordered_index(3)-1,ordered_index(3)];
ordered_Point_outer_format2 = Point_outer_format2(:,ordered_index);


%% improve the corner point with Haris 
template = zeros(size(temp_QR_binary_invert));
Harris_area = [min(ordered_Point_outer_format2(:,1))-5,max(ordered_Point_outer_format2(:,1))+5,min(ordered_Point_outer_format2(:,2))-5,max(ordered_Point_outer_format2(:,2))+5;...
    min(ordered_Point_outer_format2(:,3))-5,max(ordered_Point_outer_format2(:,3))+5,min(ordered_Point_outer_format2(:,4))-5,max(ordered_Point_outer_format2(:,4))+5;...
    min(ordered_Point_outer_format2(:,5))-5,max(ordered_Point_outer_format2(:,5))+5,min(ordered_Point_outer_format2(:,6))-5,max(ordered_Point_outer_format2(:,6))+5];
temp_QR_binary_invert_H1 = template;
temp_QR_binary_invert_H1(Harris_area(1,3):Harris_area(1,4),Harris_area(1,1):Harris_area(1,2)) = 1;
temp_QR_binary_invert_H2 = template;
temp_QR_binary_invert_H2(Harris_area(2,3):Harris_area(2,4),Harris_area(2,1):Harris_area(2,2)) = 1;
temp_QR_binary_invert_H3 = template;
temp_QR_binary_invert_H2(Harris_area(3,3):Harris_area(3,4),Harris_area(3,1):Harris_area(3,2)) = 1;
QR_Corner_1 = corner(temp_QR_binary_invert.*temp_QR_binary_invert_H1,'Harris','QualityLevel',0.3);
QR_Corner_2 = corner(temp_QR_binary_invert.*temp_QR_binary_invert_H2,'Harris','QualityLevel',0.3);
QR_Corner_3 = corner(temp_QR_binary_invert.*temp_QR_binary_invert_H3,'Harris','QualityLevel',0.3);
QR_Corner = [QR_Corner_1;QR_Corner_2;QR_Corner_3];
sum_QR_Corner = sum (QR_Corner,2);
new_QR_Corner = repmat(QR_Corner,1,4);
I_O_points_outer_block = zeros(4,6);

for i=1:3
temp_row = ordered_Point_outer_format2(:,2*i-1:2*i)';
temp_row = temp_row(:)';
new_4_point = repmat(temp_row,size(new_QR_Corner,1),1);
result = abs(new_QR_Corner - new_4_point);
result=[sum(result(:,1:2),2),sum(result(:,3:4),2),sum(result(:,5:6),2),sum(result(:,7:8),2)];
[~,I1]=min(result(:,1));
[~,I2]=min(result(:,2));
[~,I3]=min(result(:,3));
[~,I4]=min(result(:,4));
new_points_box2 = [I1,I2,I3,I4,I1];
I_O_points_outer_block(:,2*i-1:2*i) = QR_Corner(new_points_box2(1:4),:);
temp_plot_data = QR_Corner(new_points_box2,:);
end




%% try to find the far point
%hold on; plot(I_O_points_outer_block(:,5),I_O_points_outer_block(:,6),'g*');
I_O_center = sum(I_O_points_outer_block(1:4,:))/4;
I_O_center = sum([I_O_center(3:4);I_O_center(5:6)])/2;

P_farpoint = [F_bounadry_X(2),F_bounadry_Y(2)];

if QR_Orient ==1
    x1 = [I_O_points_outer_block(4,3),I_O_points_outer_block(3,3)];
    y1 = [I_O_points_outer_block(4,4),I_O_points_outer_block(3,4)];
    x2 = [I_O_points_outer_block(2,5),I_O_points_outer_block(3,5)];
    y2 = [I_O_points_outer_block(2,6),I_O_points_outer_block(3,6)];
    x3 = [I_O_points_outer_block(1,1),I_O_points_outer_block(3,1)];
    y3 = [I_O_points_outer_block(1,2),I_O_points_outer_block(3,2)];
    Lx = min(x3);
    Ly = min(y3);
    Lx_2 = max(x3);
    Ly_2 = max(y3);
    Bx = I_O_points_outer_block(4,3);
    By = I_O_points_outer_block(4,4);
    Rx = I_O_points_outer_block(2,5);
    Ry = I_O_points_outer_block(2,6);
elseif QR_Orient ==2
    x1 = [I_O_points_outer_block(1,3),I_O_points_outer_block(4,3)];
    y1 = [I_O_points_outer_block(1,4),I_O_points_outer_block(4,4)];
    x2 = [I_O_points_outer_block(3,5),I_O_points_outer_block(4,5)];
    y2 = [I_O_points_outer_block(3,6),I_O_points_outer_block(4,6)];
    x3 = [I_O_points_outer_block(2,1),I_O_points_outer_block(4,1)];
    y3 = [I_O_points_outer_block(2,2),I_O_points_outer_block(4,2)];
    Lx = max(x3);
    Ly = min(y3);
    Lx_2 = min(x3);
    Ly_2 = max(y3);
    Bx = I_O_points_outer_block(1,3);
    By = I_O_points_outer_block(1,4);
    Rx = I_O_points_outer_block(3,5);
    Ry = I_O_points_outer_block(3,6);
elseif QR_Orient ==3
    x1 = [I_O_points_outer_block(2,3),I_O_points_outer_block(1,3)];
    y1 = [I_O_points_outer_block(2,4),I_O_points_outer_block(1,4)];
    x2 = [I_O_points_outer_block(4,5),I_O_points_outer_block(1,5)];
    y2 = [I_O_points_outer_block(4,6),I_O_points_outer_block(1,6)];
    x3 = [I_O_points_outer_block(1,1),I_O_points_outer_block(3,1)];
    y3 = [I_O_points_outer_block(1,2),I_O_points_outer_block(3,2)];
    Lx = max(x3);
    Ly = max(y3);
    Lx_2 = min(x3);
    Ly_2 = min(y3);
    Bx = I_O_points_outer_block(2,3);
    By = I_O_points_outer_block(2,4);
    Rx = I_O_points_outer_block(4,5);
    Ry = I_O_points_outer_block(4,6);
elseif QR_Orient ==4
    x1 = [I_O_points_outer_block(3,3),I_O_points_outer_block(2,3)];
    y1 = [I_O_points_outer_block(3,4),I_O_points_outer_block(2,4)];
    x2 = [I_O_points_outer_block(1,5),I_O_points_outer_block(2,5)];
    y2 = [I_O_points_outer_block(1,6),I_O_points_outer_block(2,6)];
    x3 = [I_O_points_outer_block(2,1),I_O_points_outer_block(4,1)];
    y3 = [I_O_points_outer_block(2,2),I_O_points_outer_block(4,2)];
    Lx = min(x3);
    Ly = max(y3);
    Lx_2 = max(x3);
    Ly_2 = min(y3);
    Bx = I_O_points_outer_block(3,3);
    By = I_O_points_outer_block(3,4);
    Rx = I_O_points_outer_block(1,5);
    Ry = I_O_points_outer_block(1,6);
else 
    disp('direction wrong');
end
    

%line1, line2, line3
if x1(1)==x1(2)
    x1(2)=x1(1)+0.00001;
elseif x2(1)==x2(2)
    x2(2)=x2(1)+0.00001;
elseif x3(1)==x3(2)
    x3(2)=x3(1)+0.00001;
end
p1 = polyfit(x1,y1,1);
p2 = polyfit(x2,y2,1);
p3 = polyfit(x3,y3,1);
%calculate intersection
x_intersect_1 = fzero(@(x) polyval(p1-p2,x),3);
y_intersect_1 = polyval(p1,x_intersect_1);

x_intersect_2 = fzero(@(x) polyval(p3-p1,x),3);
y_intersect_2 = polyval(p3,x_intersect_2);

x_intersect_3 = fzero(@(x) polyval(p3-p2,x),3);
y_intersect_3 = polyval(p3,x_intersect_3);

x_intersect_4 = Lx+2*(I_O_center(1)-Lx);
y_intersect_4 = Ly+2*(I_O_center(2)-Ly);

x_intersect_5 = mean([x_intersect_1,x_intersect_2,x_intersect_3,x_intersect_4]);
y_intersect_5 = mean([y_intersect_1,y_intersect_2,y_intersect_3,y_intersect_4]);

x_intersect_F = round([x_intersect_5,y_intersect_5]);
%hold on;
%line(x1,y1);
%line(x2,y2);
%line(x3,y3)
%plot(round(x_intersect_F(1)),round(x_intersect_F(2)),'r*');
%plot(x_intersect_4,y_intersect_4,'b*');
%plot(x_intersect_1,y_intersect_1,'g*');
%plot(x_intersect_3,y_intersect_3,'g*');

%% rotation correction 
QR_Length_mean = mean([pdist([Lx,Ly;Bx,By]);pdist([Bx,By;x_intersect_F(1),x_intersect_F(2)]);...
    pdist([x_intersect_F(1),x_intersect_F(2);Rx,Ry]);pdist([Rx,Ry;Lx,Ly])]);

I_O_points_outer_block(5,:) = I_O_points_outer_block(1,:);
QR_outer_length =zeros(4,3);
for i=1:size(QR_outer_length,1)
QR_outer_length(i,:) = [pdist(I_O_points_outer_block(i:i+1,1:2)),...
    pdist(I_O_points_outer_block(i:i+1,3:4)),...
    pdist(I_O_points_outer_block(i:i+1,5:6))];
end
QR_outer_length_mean = sum(sum(QR_outer_length))/12;
version_size = round(QR_Length_mean/round(QR_outer_length_mean/7));
%check version size 
remind = mod(version_size,4);
if mod(version_size,4)~=1
    inter = fix(version_size/4);
    version_size = inter*4+1;
end
old_homo_set = [[Lx;Ly;1], [Bx;By;1], [x_intersect_F(1);x_intersect_F(2);1],[Rx;Ry;1]];
new_homo_set = [[1;1;1],[1;version_size*10;1],[version_size*10;version_size*10;1],[version_size*10;1;1]];
H = homography2d(old_homo_set,new_homo_set);
T = projective2d(H');
Correct_QR = imwarp(255-QR,T);
imshow(Correct_QR);
%% calculate the thread value
Gray_Correct_QR = rgb2gray(Correct_QR);
total_Correct_QR = sum(sum(Gray_Correct_QR));
total_unit_number = size(find(Gray_Correct_QR),1);
%Thread_Correct_QR = total_Correct_QR/total_unit_number/255;
Thread_Correct_QR = graythresh(Correct_QR);
%% convert to bw again 

Correct_QR_b = im2bw(Gray_Correct_QR,Thread_Correct_QR);

temp_C_QR_binary_invert = Correct_QR_b;
CCC = bwconncomp(temp_C_QR_binary_invert,8);
numPixels = cellfun(@numel,CCC.PixelIdxList);
[~,idx] = max(numPixels);
temp_C_QR_binary_invert(CCC.PixelIdxList{idx}) = 0;
%temp_C_QR_binary_invert_Corner_P = corner(temp_C_QR_binary_invert,'Harris');

%C_patch_center = round([mean(temp_C_QR_binary_invert_Corner_P(:,1)), mean(temp_C_QR_binary_invert_Corner_P(:,2))]);
%C_patch_X = (max(temp_C_QR_binary_invert_Corner_P(:,1))- min(temp_C_QR_binary_invert_Corner_P(:,1)))/10;
%C_patch_Y = (max(temp_C_QR_binary_invert_Corner_P(:,2))- min(temp_C_QR_binary_invert_Corner_P(:,2)))/10;
%temp_C_QR_binary_invert(C_patch_center(1):size(temp_C_QR_binary_invert,2),C_patch_center(2):size(temp_C_QR_binary_invert,1)) = 0;

L = bwlabel(temp_C_QR_binary_invert);
R = regionprops(L,'Area','BoundingBox','PixelList');
NR = numel(R);

cell_R = struct2cell(R);
mat_R_Area = cell2mat(cell_R(1,:));
[~,mat_R_Area_index] = sort(mat_R_Area,'descend');

ALL_Point_outer = zeros(3,8);
ALL_Point_outer_format2 = zeros(4,6);
ALL_Point_outer_index = zeros(3,5);

for i=1:NR
BBmax = R(i).BoundingBox;
DIAG1 = sum(R(i).PixelList,2);
DIAG2 = diff(R(i).PixelList,[],2);

[~,dUL] = min(DIAG1);    [~,dDR] = max(DIAG1);
[~,dDL] = min(DIAG2);    [m,dUR] = max(DIAG2);

pts = R(i).PixelList([dUL dDL dDR dUR dUL],:);
ALL_Point_outer(i,:) = [pts(1,:),pts(2,:),pts(3,:),pts(4,:)];
ALL_Point_outer_format2(:,2*i-1:2*i) = pts(1:4,:);
ALL_Point_outer_index(i,:) = [i dUL dDL dDR dUR];
end
detect_point_set = [vec2mat(ALL_Point_outer_format2(1,:),2);vec2mat(ALL_Point_outer_format2(2,:),2);...
    vec2mat(ALL_Point_outer_format2(3,:),2);vec2mat(ALL_Point_outer_format2(4,:),2)];
find_num = 0; 
ALL_Point_outer_format2(5,:)= ALL_Point_outer_format2(1,:);
C_Point_outer_format2 = zeros(4,6);
for i=1:size(mat_R_Area_index,2)
    edge_shape_index = mat_R_Area_index(i);
    edge_set = [ALL_Point_outer_format2(:,edge_shape_index*2-1),ALL_Point_outer_format2(:,edge_shape_index*2)];
    [in,on] = inpolygon(detect_point_set(:,1),detect_point_set(:,2),...
        edge_set(:,1),edge_set(:,2));
   
    if sum(in)-sum(on) == 4      
        fon = find(on);
        in(fon) = 0;
        in_points = detect_point_set(in,:);
        edge_aera = polyarea(edge_set(1:4,1)',edge_set(1:4,2)');
        in_aera = polyarea(in_points(1:4,1)',in_points(1:4,2)');
        area_ratio = edge_aera/in_aera;
        
        if area_ratio >4 && area_ratio<6
            disp(area_ratio);
            find_num = find_num+1;
             plot(edge_set(:,1),edge_set(:,2),'y','linewidth',3);
            C_Point_outer_format2(:,find_num*2-1:find_num*2) = edge_set(1:4,:);
            if find_num == 3
                break
            end
        end
    end
end 

%% improve the size again 
%C_old_homo_set = [[C_Point_outer_format2(1,1);C_Point_outer_format2(1,2);1],[C_Point_outer_format2(2,1);C_Point_outer_format2(2,2);1],[C_Point_outer_format2(3,1);C_Point_outer_format2(3,2);1],[C_Point_outer_format2(4,1);C_Point_outer_format2(4,2);1]];
%C_new_homo_set = [[1;1;1],[70;1;1],[70;70;1],[70;1;1]];
%C_H = homography2d(C_old_homo_set,C_new_homo_set);
%C_T = projective2d(C_H');
%new_Correct_QR = imwarp(Correct_QR,T);
%imshow(new_Correct_QR);

%% find LU point
mean_C_Point_outer_format2 = mean(C_Point_outer_format2);
mean_C_Point_outer_format2 = vec2mat(mean_C_Point_outer_format2,2);
[~,LU_index] = min(sum(mean_C_Point_outer_format2,2));
[~,LD_index] = max(mean_C_Point_outer_format2(:,2));
[~,RU_index] = max(mean_C_Point_outer_format2(:,1));
C_Point_outer_format2 = C_Point_outer_format2(:,[2*LU_index-1,2*LU_index,2*LD_index-1,2*LD_index,2*RU_index-1,2*RU_index]);

%% detect the size of the gird need improve already know the version size
C_QR_bound_box_pos = [C_Point_outer_format2(1,1:2);...
    C_Point_outer_format2(2,5:6);...
    C_Point_outer_format2(4,3:4)];
C_QR_bound_box = [C_Point_outer_format2(1,1:2),pdist([C_Point_outer_format2(1,1:2);C_Point_outer_format2(2,5:6)]),...
    pdist([C_Point_outer_format2(1,1:2);C_Point_outer_format2(4,3:4)])];
C_QR_bound_box(1,3:4) = round([mean(C_QR_bound_box(1,3:4)),mean(C_QR_bound_box(1,3:4))]);
C_Point_outer_format2(5,:) = C_Point_outer_format2(1,:);
C_outer_length =zeros(4,3);
for i=1:size(C_outer_length,1)
C_outer_length(i,:) = [pdist(C_Point_outer_format2(i:i+1,1:2)),pdist(C_Point_outer_format2(i:i+1,3:4)),pdist(C_Point_outer_format2(i:i+1,5:6))];
end
C_outer_length_mean = round(sum(sum(C_outer_length))/12);
C_grid_size = round(C_outer_length_mean/7);
C_check = mod(C_QR_bound_box(3),C_grid_size);
C_grid_num = round(C_QR_bound_box(3)/C_grid_size);
if C_check~=0
    C_QR_bound_box(3) = round(C_grid_num)*round(C_grid_size);
end

%% test area
%C_QR_Corner = corner(Croped_Correct_QR_b,'Harris',Harris_corner_detect_num);

%imshow(Croped_Correct_QR_b);
%hold on
%plot(C_QR_Corner(:,1), C_QR_Corner(:,2), 'r*');
%C_boundary_edge_set = [C_QR_bound_box(1,1:2);C_QR_bound_box(1),C_QR_bound_box(2)+C_QR_bound_box(4);...
%    C_QR_bound_box(1)+C_QR_bound_box(3),C_QR_bound_box(2)+C_QR_bound_box(4);C_QR_bound_box(1)+C_QR_bound_box(3),C_QR_bound_box(2);...
%    C_QR_bound_box(1,1:2)];

%[in_C_QR_Corner,on_C_QR_Corner] = inpolygon(C_QR_Corner(:,1),C_QR_Corner(:,2),...
%        C_boundary_edge_set(:,1),C_boundary_edge_set(:,2));
%  in_C_QR_Corner_list =  C_QR_Corner(in_C_QR_Corner,:) ;
%  hold on;
%plot(C_boundary_edge_set(:,1),C_boundary_edge_set(:,2),'y','linewidth',3);
%plot(in_C_QR_Corner_list(:,1), in_C_QR_Corner_list(:,2), 'g*');
% 1
%sum_in_C_QR_Corner_list = sum(in_C_QR_Corner_list,2);
%[~,max_index] = max(sum_in_C_QR_Corner_list);
% 2
%sum_C_QR_Corner_list = sum(C_QR_Corner,2);
%[~,max_index] = max(sum_C_QR_Corner_list);
%whatIwant = C_QR_Corner(max_index,:);
%plot(whatIwant(:,1), whatIwant(:,2), 'c*');
%vi = 1:21;
%vi_points = vi*10;
%vi_points_x = repmat(vi_points,1,21)';
%vi_points_y = repmat(vi_points',1,21)';
%vi_points_y = vi_points_y(:);
%plot(vi_points_x,vi_points_y, 'y*');
%imshow(CB);
%[imagePoints,boardSize] = detectCheckerboardPoints(CB);
%hold on ;
%plot(imagePoints(:,1), imagePoints(:,2), 'y*');

%% generate grid
C_grid_size_list = 1:fix(C_grid_size/3);
C_grid_1 = zeros(C_grid_size);
C_grid_2 = zeros(C_grid_size);
C_grid_3 = zeros(C_grid_size);
C_grid_4 = zeros(C_grid_size);
C_grid_5 = zeros(C_grid_size);
C_grid_1(C_grid_size_list(1)*3,C_grid_size_list(1)*3) = 1;
C_grid_2(C_grid_size_list(2)*3,C_grid_size_list(2)*3) = 1;
C_grid_3(C_grid_size_list(3)*3,C_grid_size_list(3)*3) = 1;
C_grid_4(C_grid_size_list(1)*3,C_grid_size_list(3)*3) = 1;
C_grid_5(C_grid_size_list(3)*3,C_grid_size_list(1)*3) = 1;
C_grid_sampler_1 = repmat(C_grid_1,C_grid_num);
C_grid_sampler_2 = repmat(C_grid_2,C_grid_num);
C_grid_sampler_3 = repmat(C_grid_3,C_grid_num);
C_grid_sampler_4 = repmat(C_grid_4,C_grid_num);
C_grid_sampler_5 = repmat(C_grid_5,C_grid_num);

%% extract bit information 
offset = C_QR_bound_box(3)-1;
Croped_Correct_QR_b = Correct_QR_b(C_QR_bound_box(2):C_QR_bound_box(2)+offset,C_QR_bound_box(1):C_QR_bound_box(1)+offset);
S_mat_1 = (Croped_Correct_QR_b+1).*C_grid_sampler_1;
S_mat_2 = (Croped_Correct_QR_b+1).*C_grid_sampler_2;
S_mat_3 = (Croped_Correct_QR_b+1).*C_grid_sampler_3;
S_mat_4 = (Croped_Correct_QR_b+1).*C_grid_sampler_4;
S_mat_5 = (Croped_Correct_QR_b+1).*C_grid_sampler_5;
new_S_mat_1 = S_mat_1(any(S_mat_1,2),:);
final_S_mat_1 = new_S_mat_1(:,any(S_mat_1))-1;
new_S_mat_2 = S_mat_2(any(S_mat_2,2),:);
final_S_mat_2 = new_S_mat_2(:,any(S_mat_2))-1;
new_S_mat_3 = S_mat_3(any(S_mat_3,2),:);
final_S_mat_3 = new_S_mat_3(:,any(S_mat_3))-1;
new_S_mat_4 = S_mat_4(any(S_mat_4,2),:);
final_S_mat_4 = new_S_mat_4(:,any(S_mat_4))-1;
new_S_mat_5 = S_mat_5(any(S_mat_5,2),:);
final_S_mat_5 = new_S_mat_5(:,any(S_mat_5))-1;
final_S_mat = final_S_mat_2;
%% find error correction and mask pattern 
r_mask_S_mat = final_S_mat;
r_m_size = size(r_mask_S_mat,1);
err_mask_1 = [final_S_mat(9,1:6),final_S_mat(9,8:9),final_S_mat(9,r_m_size-6:r_m_size)];
err_mask_2 = [final_S_mat(1:6,9);final_S_mat(8:9,9);final_S_mat(r_m_size-6:r_m_size,9)];
err_mask_2 = fliplr(err_mask_2');
err_mask = fix(mean([err_mask_1;err_mask_2]));
type_mask = bi2de(xor(err_mask(:,3:5),[1,0,1]),'left-msb');
type_err = bi2de(xor(err_mask(:,1:2),[1,0]),'left-msb');
if type_mask == 0
    for i=1:size(r_mask_S_mat,1)
        for j = 1:size(r_mask_S_mat,1)
            if mod(i-1+j-1,2) == 0
                r_mask_S_mat(i,j) = 1-r_mask_S_mat(i,j);
            end
        end
    end
elseif type_mask == 1
    for i = 1:size(r_mask_S_mat,1)
        if mod(i-1,2) ==0
            r_mask_S_mat(1,:) = 1-r_mask_S_mat(1,:) ;
        end
    end
elseif type_mask == 2
    for j = 1:size(r_mask_S_mat,1)
        if mod(j-1,3) ==0
            r_mask_S_mat(:,j) = 1-r_mask_S_mat(:,j) ;
        end
    end
elseif type_mask == 3
    for i=1:size(r_mask_S_mat,1)
        for j = 1:size(r_mask_S_mat,1)
            if mod(i-1+j-1,3) == 0
                r_mask_S_mat(i,j) = 1-r_mask_S_mat(i,j);
            end
        end
    end
elseif type_mask == 4
    for i=1:size(r_mask_S_mat,1)
        for j = 1:size(r_mask_S_mat,1)
            if mod(floor(i-1/2)+floor(j-1/3),2) == 0
                r_mask_S_mat(i,j) = 1-r_mask_S_mat(i,j);
            end
        end
    end
elseif type_mask == 5
    for i=1:size(r_mask_S_mat,1)
        for j = 1:size(r_mask_S_mat,1)
            if mod((i-1)*(j-1),2)+mod((i-1)*(j-1),3) == 0
                r_mask_S_mat(i,j) = 1-r_mask_S_mat(i,j);
            end
        end
    end
elseif type_mask == 6
    for i=1:size(r_mask_S_mat,1)
        for j = 1:size(r_mask_S_mat,1)
            if mod((mod((i-1)*(j-1),2)+mod((i-1)*(j-1),3)),2) == 0
                r_mask_S_mat(i,j) = 1-r_mask_S_mat(i,j);
            end
        end
    end
elseif type_mask == 7
    for i=1:size(r_mask_S_mat,1)
        for j = 1:size(r_mask_S_mat,1)
            if mod((mod(i+j-2,2)+mod((i-1)*(j-1),3)),2) == 0
                r_mask_S_mat(i,j) = 1-r_mask_S_mat(i,j);
            end
        end
    end
end

%% encoding type and err correction 

 block_Table = [1,0,16,10;1,0,19,7;1,0,9,17;1,0,13,13;...
        1,0,28,16;1,0,34,10;1,0,16,28;1,0,22,22;...
        1,0,44,26;1,0,55,15;2,0,13,44;2,0,17,36;...
        2,0,32,36;1,0,80,20;4,0,9,64;2,0,24,52;...
        2,0,43,48;1,0,108,26;2,1,11,88;2,1,15,72;...
        4,0,27,64;2,0,68,36;4,0,15,112;4,0,19,96;];
 
GP_Table_index = {'7';'10';'13';'15';'16';'17';'18';'20';'22';'24';'26';'28'};
GP_Table_content = {[1,127,122,154,164,11,68,117];...
[1,216,194,159,111,199,94,95,113,157,193];...
[1,137,73,227,17,177,17,52,13,46,218,83,132,128];...
[1,29,196,111,163,112,74,10,105,237,132,151,32,134,26];...
[1,59,13,104,189,68,209,30,8,163,65,41,229,98,50,36,59];...
[1,119,66,83,120,119,22,197,83,249,41,143,134,85,53,125,99,79];...
[1,239,251,183,113,149,175,199,215,240,220,73,82,173,75,32,67,217,146];...
[1,152,185,240,5,111,99,6,220,112,150,69,36,187,22,228,198,121,121,165,174];...
[1,89,179,247,176,182,244,19,189,69,40,28,137,29,123,67,253,86,218,230,26,145,245];...
[1,112,118,169,70,178,237,216,102,115,150,229,73,130,72,61,43,206,1,237,247,127,217,144,117];...
[1,246,51,183,4,136,98,199,152,77,56,206,24,145,40,209,117,233,42,135,68,70,144,146,77,43,94];...
[1,252,9,28,13,18,251,208,150,103,174,100,41,167,12,247,56,117,119,233,127,181,100,121,147,176,74,58,197]};
GP_Table = table(GP_Table_content,'RowNames',GP_Table_index);


version = round((size(r_mask_S_mat,1)-21)/4+1);

block_type = (version -1)*4+type_err+1;
block_info = block_Table(block_type,:);

type_UP = 1;
type_DOWN = 0;
r_m_t_S_mat = r_mask_S_mat;
r_m_t_S_mat(7,:)=[];
r_m_t_S_mat(:,7)=[];

r_m_t_size = size(r_m_t_S_mat,1);
already_added = 1;
if version == 1
    Alignment = 0;
    data_corr_code = zeros(1,208);
    
    for i=1:10
        direction = mod(i,2);
        column_temp = r_m_t_S_mat(:,21-2*i:22-2*i);
        if i<5
            column_temp = column_temp(9:20,:);
        elseif i>6
            column_temp = column_temp(9:12,:);
        end
        if direction == type_UP
            data_code = column_temp';
            data_code = fliplr(data_code(:)');
        elseif direction == type_DOWN
            data_code = [column_temp(:,2),column_temp(:,1)]';
            data_code = data_code(:)';
        end 
        endpos = already_added+size(data_code,2);
        data_corr_code(already_added:endpos-1) = data_code;
        already_added = endpos;
    end
elseif version <7 && version >1
    Alignment = 1;
    
    Align_r_t_m = [r_m_t_size-9,r_m_t_size-3];
    code_length = r_m_size*r_m_size - 9*9 -8*9*2 - 2*version*4;
    data_corr_code = zeros(1,code_length);
    already_added = 1;
    iter_num = (r_m_size-1)/2;
    for i=1:iter_num
        direction = mod(i,2);
        column_temp = r_m_t_S_mat(:,r_m_size-2*i:r_m_size+1-2*i);
        if i~=5
            if i < 3
                column_temp = column_temp(9:r_m_t_size,:);
            elseif i < 5
                column_temp = [column_temp(9:Align_r_t_m(1),:);column_temp(Align_r_t_m(2):r_m_t_size,:)];
            elseif i > iter_num - 4
                column_temp = column_temp(9:r_m_t_size-8,:);
            end
            if direction == type_UP
                data_code = column_temp';
                data_code = fliplr(data_code(:)');
            elseif direction == type_DOWN
                data_code = [column_temp(:,2),column_temp(:,1)]';
                data_code = data_code(:)';
            end 
        elseif i == 5
                column_temp_1 = column_temp(Align_r_t_m(2):r_m_t_size,:);
                column_temp_2 = column_temp(Align_r_t_m(1)+1:Align_r_t_m(2)-1,1);
                column_temp_3 = column_temp(1:Align_r_t_m(1),:);
                data_code_1 = column_temp_1';
                data_code_1 = fliplr(data_code_1(:)');
                data_code_2 = column_temp_2';
                data_code_2 = fliplr(data_code_2);
                data_code_3 = column_temp_3';
                data_code_3 = fliplr(data_code_3(:)');
                data_code = [data_code_1,data_code_2,data_code_3];
        end
        endpos = already_added+size(data_code,2);
        data_corr_code(already_added:endpos-1) = data_code;
        already_added = endpos;
    end
else
    disp('version large than 6 is not support yet');
end

if block_info(1) == 1
    data_need_corr = data_corr_code(:,1:block_info(3)*8);
    data_ece = data_corr_code(:,block_info(3)*8+1:block_info(3)*8+block_info(4)*8);
    data_need_corr_mat = vec2mat(data_need_corr,8);
    data_ece_mat = vec2mat(data_ece,8);
    data_need_corr_de = bi2de(data_need_corr_mat,'left-msb');
    data_ece_de = bi2de( data_ece_mat,'left-msb');
    
    
    GP_info = GP_Table({num2str(block_info(4))},:);
    GP_info_mat = cell2mat(table2cell(GP_info));
    % err correction 
    m = 8; % Number of bits in each symbol
    n = block_info(3)+block_info(4); k = block_info(3); % Codeword length and message length
    % Simplest syntax for encoding
    data_and_ece_code_de = bi2de(vec2mat(data_corr_code(:,1:n*8),8),'left-msb');
    hDec = comm.RSDecoder(n,k);
    release(hDec)
    hDec.PrimitivePolynomialSource = 'Property';
    hDec.GeneratorPolynomialSource = 'Property';
    hDec.GeneratorPolynomial = GP_info_mat;
    hDec.PrimitivePolynomial       = [1 0 0 0 1 1 1 0 1];
    data_form_ece = step(hDec, data_and_ece_code_de);
    if isequal(data_and_ece_code_de(1:block_info(3),:),data_form_ece(1:block_info(3),:)) ~= 1
        disp('correct data');
        data_from_ece_bi = de2bi(data_form_ece(1:block_info(3),:),'left-msb');
        data_from_ece_bi_T = data_from_ece_bi';
        data_from_ece_bi_vec = data_from_ece_bi_T(:);
        data_corr_code(:,1:block_info(3)*8) = data_from_ece_bi_vec';
    end
    
elseif block_info(1) > 1
%   block_type = (version -1)*4+type_err+1;
%    block_info = block_Table(block_type,:);    
        if block_info(2)>0
            reformat_data_size = block_info(1)*block_info(3)+block_info(1)*(block_info(3)+1);
            old_index_list = 1:reformat_data_size;
            old_index_list_1 = old_index_list(1:reformat_data_size-2);
            old_index_list_2 = old_index_list(reformat_data_size-1:reformat_data_size);
            new_index_mat = vec2mat(old_index_list_1,block_info(1)*2);
            new_index_mat_1 = new_index_mat(:,1:block_info(1));
            new_index_mat_2 = new_index_mat(:,2*block_info(1)-1,2*block_info(1));
            new_index_mat_2 = [new_index_mat_2;old_index_list_2];
            new_index_1 = new_index_mat_1(:);
            new_index_2 = new_index_mat_2(:);
            new_index_list = [new_index_1;new_index_2];
            new_data_mat = vec2mat(data_corr_code(:,1:reformat_data_size*8),8);
            new_data =new_data_mat( new_index_list,:); % useful mat 
            data_need_corr_de = bi2de(new_data,'left-msb');
            data_need_corr = new_data; 
            
            %vector
            
            reformat_ece_size = block_info(4);
            new_ece_vec = data_corr_code(:,reformat_data_size*8+1:reformat_data_size*8+reformat_ece_size*8);
            new_data_mat_ece = vec2mat(new_ece_vec,8);
            old_index_list_ece = 1:reformat_ece_size;
            new_index_mat_ece = vec2mat(old_index_list_ece,block_info(1)*2);
            new_index_list_ece = new_index_mat_ece(:);
            new_data_ece = new_data_mat_ece(new_index_list_ece,:);
            data_ece_mat = new_data_ece; %mat
            data_ece_de = bi2de(new_data_ece,'left-msb');
            
            F_T_data_index = [1,block_info(3);...
                block_info(3)+1,block_info(3)+block_info(3);...
                block_info(3)+block_info(3)+1, block_info(3)+block_info(3)+block_info(3)+1;...
                block_info(3)+block_info(3)+block_info(3)+2,block_info(3)+block_info(3)+block_info(3)+block_info(3)+2];
            ece_codeword_length = block_info(4)/(block_info(1)*2);
            F_T_ece_index = [1,ece_codeword_length;...
                ece_codeword_length+1,ece_codeword_length*2;...
                ece_codeword_length*2+1,ece_codeword_length*3;...
                ece_codeword_length*3+1,ece_codeword_length*4];
            
            % error correction 
            GP_info = GP_Table({num2str(ece_codeword_length)},:);
            GP_info_mat = cell2mat(table2cell(GP_info));
            
            for i = block_info(1)*2
                msg = [data_need_corr_de(F_T_data_index(i,1):F_T_data_index(i,2),:);data_ece_de(F_T_ece_index(i,1):F_T_ece_index(i,2),:)];
                m = 8; % Number of bits in each symbol
                n = size(msg,1); k = size(data_need_corr_de(F_T_data_index(i,1):F_T_data_index(i,2),:),1); % Codeword length and message length
                % Simplest syntax for encoding

                hDec = comm.RSDecoder(n,k);
                release(hDec)
                hDec.PrimitivePolynomialSource = 'Property';
                hDec.GeneratorPolynomialSource = 'Property';
                hDec.GeneratorPolynomial = GP_info_mat;
                hDec.PrimitivePolynomial       = [1 0 0 0 1 1 1 0 1];
                data_form_ece = step(hDec, msg);
                if isequal(data_form_ece(1,k,:),...
                        data_need_corr_de(F_T_data_index(i,1):F_T_data_index(i,2),:)) ~= 1
                    data_from_ece_bi = de2bi(data_form_ece,'left-msb');
                    data_need_corr(F_T_data_index(i,1):F_T_data_index(i,2),:) = data_from_ece_bi(1,k,:);
                end
            end
            
            data_need_corr_T = data_need_corr'; 
            data_corr_code = data_need_corr_T(:)';
            
            
        elseif block_info(2) == 0
            reformat_data_size = block_info(1)*block_info(3);
            reformat_ece_size = block_info(4);
            old_index_list = 1:reformat_data_size;
            old_index_list_ece = 1:reformat_ece_size;
            new_index_mat = vec2mat(old_index_list,block_info(1));
            new_index_mat_ece = vec2mat(old_index_list_ece,block_info(1));
            new_index_list = new_index_mat(:);
            new_index_list_ece = new_index_mat_ece(:);
            new_data_mat = vec2mat(data_corr_code(:,1:reformat_data_size*8),8);
            new_data_mat_ece = vec2mat(data_corr_code(:,reformat_data_size*8+1:reformat_data_size*8+reformat_ece_size*8),8); 
            new_data =new_data_mat( new_index_list,:);
            new_data_ece = new_data_mat_ece(new_index_list_ece,:);
            data_need_corr_mat = new_data;
            data_ece_mat = new_data_ece;
            
            
            % err correction

            data_need_corr_de = bi2de(data_need_corr_mat,'left-msb')';
            data_ece_de = bi2de( data_ece_mat,'left-msb')';
            data_need_corr_de_mat = vec2mat(data_need_corr_de,block_info(3));
            data_ece_de_mat = vec2mat(data_ece_de,block_info(4)/block_info(1));
            data_and_ece_code_de = [data_need_corr_de_mat,data_ece_de_mat]';
            %data_form_ece = zeros(size(data_and_ece_code_de));
            
            GP_info = GP_Table({num2str(block_info(4)/block_info(1))},:);
            GP_info_mat = cell2mat(table2cell(GP_info));
            % err correction 
            m = 8; % Number of bits in each symbol
            n = block_info(3)+block_info(4)/block_info(1); k = block_info(3); % Codeword length and message length
            % Simplest syntax for encoding
            
            hDec = comm.RSDecoder(n,k);
            release(hDec)
            hDec.PrimitivePolynomialSource = 'Property';
            hDec.GeneratorPolynomialSource = 'Property';
            hDec.GeneratorPolynomial = GP_info_mat;
            hDec.PrimitivePolynomial       = [1 0 0 0 1 1 1 0 1];
            for i = 1:size(data_and_ece_code_de,2)
                data_form_ece = step(hDec, data_and_ece_code_de(:,i));                
                if isequal(data_need_corr_de_mat(i,:)',data_form_ece) ~= 1
                    disp(i);
                    data_from_ece_bi = de2bi(data_form_ece,'left-msb');
                    data_need_corr_mat(i*block_info(3)-block_info(3)+1:i*block_info(3),:) = data_from_ece_bi;
                end
            end
            
            data_need_corr_mat_T = data_need_corr_mat';
            data_ece_mat_T = data_ece_mat';
            ece_data = data_ece_mat_T(:)';
            data_corr_code = data_need_corr_mat_T(:)';
        end

end 
%% the data code has been changed !!!!

%% err correction test
%data_cd_len = size(data_corr_code,2)/8-block_info(4);
%ec_cd_len = block_info(4);
%data_cd = data_corr_code(1:data_cd_len*8);
%ec_cd = data_corr_code((data_cd_len*8+1):size(data_corr_code,2));
%data_cd_de = bi2de(vec2mat(data_cd,8),'left-msb');
%ec_cd_de = bi2de(vec2mat(ec_cd,8),'left-msb');
%hDec = comm.RSDecoder(255,size(data_cd_de,1));
%[dc nerrs] = step(hDec, ec_cd_de);
%% after get the bit array 
encoding_type = data_corr_code(:,1:4);
%encoding_type = fliplr(encoding_type(:)');
de_encoding_type = bi2de(encoding_type,'left-msb');
if de_encoding_type == 1 
    length_data_bin = data_corr_code(5:14);
    length_data = bi2de(length_data_bin,'left-msb');
    i_part = fix(length_data/3);
    r_part = mod(length_data,3);
    data_area_i = [15,14+i_part*10];  
     data_bin_i = data_corr_code(:,data_area_i(1):data_area_i(2));
    data_bin_mat_i = vec2mat(data_bin_i,10);
    data_de_i = bi2de(data_bin_mat_i,'left-msb');  
    data_de = data_de_i;
    times = 0:(i_part-1);
    ratio = 1000.^times;
    ratio = fliplr(ratio);
    data_str = ratio*data_de;
    
    if r_part == 1
        data_area_r = [data_area_i(2)+1,data_area_i(2)+4];
        data_bin_r = data_corr_code(:,data_area_r(1):data_area_r(2));
        data_de_r = bi2de(data_bin_r,'left-msb');
        data_str = data_str*10+data_de_r;
    elseif r_part == 2 
        data_area_r = [data_area_i(2)+1,data_area_i(2)+7];
        data_bin_r = data_corr_code(:,data_area_r(1):data_area_r(2));
        data_de_r = bi2de(data_bin_r,'left-msb');  
        data_str = data_str*100+data_de_r;
    end
    
elseif de_encoding_type == 2
    Alp_Table = {'0';'1';'2';'3';'4';'5';'6';'7';'8';'9';...
        'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z';...
        ' ';'$';'%';'*';'+';'_';'.';'/';':'};
    length_data_bin = data_corr_code(5:13);
    length_data = bi2de(length_data_bin,'left-msb');
    i_part = fix(length_data/2);
    r_part = mod(length_data,2);
    data_area_i = [14,13+i_part*11];   
    data_bin_i = data_corr_code(:,data_area_i(1):data_area_i(2));
    data_bin_mat_i = vec2mat(data_bin_i,11);
    data_de_i = bi2de(data_bin_mat_i,'left-msb');
    data_de_i_pair = [fix(data_de_i./45),mod(data_de_i,45)]+1;
    data_str_i = Alp_Table(data_de_i_pair)';
    data_str_i = data_str_i(:)';
    data_str = cell2mat(data_str_i);
    if r_part ~= 0
    data_area_r = [data_area_i(2)+1,data_area_i(2)+6];
    data_bin_r = data_corr_code(:,data_area_r(1):data_area_r(2));
    data_de_r = bi2de(data_bin_r,'left-msb')+1;
    data_str_r = Alp_Table(data_de_r);
    data_str = [data_str,cell2mat(data_str_r)];
    end
    
elseif de_encoding_type == 4
    length_data_bin = data_corr_code(5:12);
    length_data = bi2de(length_data_bin,'left-msb');
    data_area = [13,12+length_data*8];
    data_bin = data_corr_code(:,data_area(1):data_area(2));
    data_bin_mat = vec2mat(data_bin,8);
    data_de = bi2de(data_bin_mat,'left-msb');
    data_str = native2unicode(data_de)';
end

disp(data_str);
%% test area



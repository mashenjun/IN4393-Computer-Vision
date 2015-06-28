javaaddpath('.\zxing-core-1.6.jar');
javaaddpath('.\zxing-javase.jar');
import com.google.zxing.qrcode.*;
import com.google.zxing.client.j2se.BufferedImageLuminanceSource;
import com.google.zxing.BinaryBitmap;
import com.google.zxing.common.HybridBinarizer;
import com.google.zxing.Result;


QR = imread('input/input5.jpg');
level = 0.5;
QR_binary_invert = im2bw(255-QR,level);
QR_binary = im2bw(QR,level);
imshow(QR);

%% parameter setting here
Harris_corner_detect_num = 200;


%QR = imread('input.jpg');
%qr_reader = QRCodeReader();
%jimg = im2java2d(QR);
%source = BufferedImageLuminanceSource(jimg);
%bitmap = BinaryBitmap(HybridBinarizer(source));
%Text = bitmap.toString();
%result = qr_reader.decode(bitmap); 
%message = char(result.getText());

%% find corner points
QR_Corner = corner(QR_binary,'Harris',1000);
imshow(QR);
hold on
plot(QR_Corner(:,1), QR_Corner(:,2), 'r*');
[Temp,Corner_temp] = sort(QR_Corner(:,1));
Sorted_Corner_temp = QR_Corner(Corner_temp,:);

%% edge detection 
BW1 = edge(QR_binary_invert,'sobel');
BW2 = edge(QR_binary_invert,'canny'); % better
figure;
imshowpair(BW1,BW2,'montage')
title('Sobel Filter                                   Canny Filter');

%% contour detection
[Contour,h] = contour(QR);
clabel(Contour);

%% centroids /useless 
s = regionprops(QR,'centroid');
centroids = cat(1, s.Centroid);
imshow(QR)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')

%% after getting the edge 
[edge_point_Y,edge_point_X] = find(BW2);
edge_point = [edge_point_X,edge_point_Y];
[counts2_X,center_X] = hist( edge_point_X, unique(edge_point_X));
[counts2_Y,center_Y] = hist( edge_point_Y, unique(edge_point_Y));
count_X = [counts2_X;center_X'];
count_Y = [counts2_Y;center_Y';];
[Temp_X,index_X] = sort(count_X(1,:),'descend');
[Temp_Y,index_Y] = sort(count_Y(1,:),'descend');
Sorted_count_X = count_X(:,index_X);
Sorted_count_Y = count_Y(:,index_Y);
bound = Ma_findbound(QR_binary_invert);
for i = bound(:,1) : bound(:,1)
    line = QR_binary_invert(i,:);
    
end

%% useful here

QR = imread('input/input1.jpg');
level = 0.5;
QR_binary_invert = im2bw(255-QR,level);
QR_binary = im2bw(QR,level);
imshow(QR);

temp_QR_binary_invert = QR_binary_invert;
CC = bwconncomp(temp_QR_binary_invert);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
temp_QR_binary_invert(CC.PixelIdxList{idx}) = 0;

imshow(temp_QR_binary_invert);
QR_Corner = corner(temp_QR_binary_invert,'Harris',Harris_corner_detect_num);
hold on
plot(QR_Corner(:,1), QR_Corner(:,2), 'r*');

hold on;
L = bwlabel(temp_QR_binary_invert);
R = regionprops(L,'Area','BoundingBox','PixelList');
NR = numel(R);

maxArea = 0;
index = 1;
cell_R = struct2cell(R);
mat_R_Area = cell2mat(cell_R(1,:));
mat_R_Area = sort(mat_R_Area,'descend');
thread_Area = mat_R_Area(4);
kmax = zeros(1,3);
A = zeros(1,13);

for k = 1:NR
    A(k) = prod(R(k).BoundingBox(3:4));
    if R(k).Area > thread_Area
        maxArea = R(k).Area;
        kmax(index) = k;
        index = index+1;
    end
end

XYLIMS = zeros(3,4);
Point_outer = zeros(3,8);
Point_outer_format2 = zeros(4,6);
Point_outer_index = zeros(3,5);

for i=1:size(kmax,2)
BBmax = R(kmax(i)).BoundingBox;
DIAG1 = sum(R(kmax(i)).PixelList,2);
DIAG2 = diff(R(kmax(i)).PixelList,[],2);

[~,dUL] = min(DIAG1);    [~,dDR] = max(DIAG1);
[~,dDL] = min(DIAG2);    [m,dUR] = max(DIAG2);

pts = R(kmax(i)).PixelList([dUL dDL dDR dUR dUL],:);
Point_outer(i,:) = [pts(1,:),pts(2,:),pts(3,:),pts(4,:)];
Point_outer_format2(:,2*i-1:2*i) = pts(1:4,:);
Point_outer_index(i,:) = [kmax(i) dUL dDL dDR dUR];
h_pts = plot(pts(:,1),pts(:,2),'m','linewidth',3);

XYLIMS(i,:) = [BBmax(1) + [0 BBmax(3)] BBmax(2) + [0 BBmax(4)]];
end
%% try connected components 

imshow(BW2);

CC = bwconncomp(BW2);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
BW2(CC.PixelIdxList{idx}) = 0;
 
figure, imshow(BW2);

%% oeder 3 block
% some value needed 
total_count_QR_invert_b = sum(sum(QR_binary_invert));
P2X_QR_invert_b = sum(QR_binary_invert);
P2Y_QR_invert_b = sum(QR_binary_invert,2);
PrL_Y = find(P2Y_QR_invert_b,1,'last')- find(P2Y_QR_invert_b,1,'first');
PrL_X = find(P2X_QR_invert_b,1,'last')- find(P2X_QR_invert_b,1,'first');
Thread2X_find_boundary = sqrt(total_count_QR_invert_b/PrL_X);
Thread2Y_find_boundary = sqrt(total_count_QR_invert_b/PrL_Y);
F_bounadry_X = [find(P2X_QR_invert_b>Thread2X_find_boundary,1,'first'),find(P2X_QR_invert_b>Thread2X_find_boundary,1,'last')];
F_bounadry_Y = [find(P2Y_QR_invert_b>Thread2Y_find_boundary,1,'first'),find(P2Y_QR_invert_b>Thread2Y_find_boundary,1,'last')];
F_boundry_Center = [mean(F_bounadry_X),mean(F_bounadry_Y)];
Center_outer_block = sum(Point_outer_format2)/4;
% change row to colum
Center_outer_block = reshape(Center_outer_block,[2,3])';
dis = [pdist([mean(Center_outer_block(2:3,:));F_boundry_Center]);...
    pdist([mean([Center_outer_block(1,:);Center_outer_block(3,:)]);F_boundry_Center]);...
    pdist([mean(Center_outer_block(1:2,:));F_boundry_Center])];
%dis = [pdist(Center_outer_block(2:3,:));...
%    pdist([Center_outer_block(1,:);Center_outer_block(3,:)]);...
%    pdist(Center_outer_block(1:2,:))];
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

sum_QR_Corner = sum (QR_Corner,2);
new_QR_Corner = repmat(QR_Corner,1,4);
I_O_points_outer_block = zeros(4,6);
hold on ;
for i=1:3
temp_row = ordered_Point_outer_format2(:,2*i-1:2*i)';
temp_row = temp_row(:)';
new_4_point = repmat(temp_row,Harris_corner_detect_num,1);
result = abs(new_QR_Corner - new_4_point);
result=[sum(result(:,1:2),2),sum(result(:,3:4),2),sum(result(:,5:6),2),sum(result(:,7:8),2)];
[~,I1]=min(result(:,1));
[~,I2]=min(result(:,2));
[~,I3]=min(result(:,3));
[~,I4]=min(result(:,4));
new_points_box2 = [I1,I2,I3,I4,I1];
I_O_points_outer_block(:,2*i-1:2*i) = QR_Corner(new_points_box2(1:4),:);
temp_plot_data = QR_Corner(new_points_box2,:);
plot(temp_plot_data(:,1),temp_plot_data(:,2),'b','linewidth',1);
end




%% try to find the far point
%Hsum = sum(QR_Corner,2);
%Hdiff= diff(QR_Corner,[],2);
%[~,HUL] = min(Hsum);    [~,HDR] = max(Hsum);
%[~,HDL] = min(Hdiff);    [m,HUR] = max(Hdiff);
%Hts = QR_Corner([HUL HDL HDR HUR HUL],:);
%hold on;
%plot(Hts(:,1), Hts(:,2), 'b*');
I_O_center = sum(I_O_points_outer_block)/4;
I_O_center = sum([I_O_center(3:4);I_O_center(5:6)])/2;
x3 = [I_O_points_outer_block(1,1),I_O_points_outer_block(3,1)];
y3 = [I_O_points_outer_block(1,2),I_O_points_outer_block(3,2)];
P_farpoint = [F_bounadry_X(2),F_bounadry_Y(2)];
if QR_Orient ==1
    x1 = [I_O_points_outer_block(4,3),I_O_points_outer_block(3,3)];
    y1 = [I_O_points_outer_block(4,4),I_O_points_outer_block(3,4)];
    x2 = [I_O_points_outer_block(2,5),I_O_points_outer_block(3,5)];
    y2 = [I_O_points_outer_block(2,6),I_O_points_outer_block(3,6)];
    Lx = min(x3);
    Ly = min(y3);
    Bx = I_O_points_outer_block(4,3);
    By = I_O_points_outer_block(4,4);
    Rx = I_O_points_outer_block(2,5);
    Ry = I_O_points_outer_block(2,6);
elseif QR_Orient ==2
    x1 = [I_O_points_outer_block(1,3),I_O_points_outer_block(4,3)];
    y1 = [I_O_points_outer_block(1,4),I_O_points_outer_block(4,4)];
    x2 = [I_O_points_outer_block(3,5),I_O_points_outer_block(4,5)];
    y2 = [I_O_points_outer_block(3,6),I_O_points_outer_block(4,6)];
    Lx = max(x3);
    Ly = min(y3);
    Bx = I_O_points_outer_block(1,3);
    By = I_O_points_outer_block(1,4);
    Rx = I_O_points_outer_block(3,5);
    Ry = I_O_points_outer_block(3,6);
elseif QR_Orient ==3
    x1 = [I_O_points_outer_block(2,3),I_O_points_outer_block(1,3)];
    y1 = [I_O_points_outer_block(2,4),I_O_points_outer_block(1,4)];
    x2 = [I_O_points_outer_block(4,5),I_O_points_outer_block(1,5)];
    y2 = [I_O_points_outer_block(4,6),I_O_points_outer_block(1,6)];
    Lx = max(x3);
    Ly = max(y3);
    Bx = I_O_points_outer_block(2,3);
    By = I_O_points_outer_block(2,4);
    Rx = I_O_points_outer_block(4,5);
    Ry = I_O_points_outer_block(4,6);
elseif QR_Orient ==4
    x1 = [I_O_points_outer_block(3,3),I_O_points_outer_block(2,3)];
    y1 = [I_O_points_outer_block(3,4),I_O_points_outer_block(2,4)];
    x2 = [I_O_points_outer_block(1,5),I_O_points_outer_block(2,5)];
    y2 = [I_O_points_outer_block(1,6),I_O_points_outer_block(2,6)];
    Lx = min(x3);
    Ly = max(y3);
    Bx = I_O_points_outer_block(3,3);
    By = I_O_points_outer_block(3,4);
    Rx = I_O_points_outer_block(1,5);
    Ry = I_O_points_outer_block(1,6);
else 
    disp('direction wrong');
end
    

 %line1
%x1  = [Point_outer(3,3) Point_outer(3,5)];
%y1  = [Point_outer(3,4) Point_outer(3,6)];
%line2
%x2 = [Point_outer(2,5) Point_outer(2,7)];
%y2 = [Point_outer(2,6) Point_outer(2,8)];
%fit linear polynomial
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

x_intersect_5 = mean([x_intersect_1,x_intersect_2,x_intersect_3,x_intersect_4,P_farpoint(1)]);
y_intersect_5 = mean([y_intersect_1,y_intersect_2,y_intersect_3,y_intersect_4,P_farpoint(2)]);

x_intersect_list = [x_intersect_1;x_intersect_2;x_intersect_3;x_intersect_4;x_intersect_5];
y_intersect_list = [y_intersect_1;y_intersect_2;y_intersect_3;y_intersect_4;y_intersect_5];
intersect_list = [x_intersect_list,y_intersect_list];
dist_intersect_list_dist = [pdist([intersect_list(1,:);P_farpoint]);...
    pdist([intersect_list(2,:);P_farpoint]);...
    pdist([intersect_list(3,:);P_farpoint]);...
    pdist([intersect_list(4,:);P_farpoint]);...
    pdist([intersect_list(5,:);P_farpoint])];
[~,x_intersect_F_index] = min(dist_intersect_list_dist);
x_intersect_F = intersect_list(x_intersect_F_index,:);

line(x1,y1);
hold on;
line(x2,y2);
line(x3,y3)
hold on;
plot(x_intersect_F(1),x_intersect_F(2),'k*');

%% rotation correction 
old_homo_set = [[Lx;Ly;1], [Bx;By;1], [x_intersect_F(1);x_intersect_F(2);1],[Rx;Ry;1]];
new_homo_set = [[1;1;1],[1;200;1],[200;200;1],[200;1;1]];
H = homography2d(old_homo_set,new_homo_set);
old_set_T = old_homo_set(1:2,:)';
X_boundary = round([min(old_set_T(:,1))-1,max(old_set_T(:,1))+1]);
Y_boundary = round([min(old_set_T(:,2))-1,max(old_set_T(:,2))]+1);
MAX_boundary_dis = max(diff([X_boundary;Y_boundary],[],2));

T = projective2d(H');
%Correct_QR_b = imwarp(QR_binary_invert,T);
Correct_QR = imwarp(255-QR,T);

imshow(Correct_QR);
%% crop center image
sample_C_QR = Correct_QR(1:265,1:180);
total_Correct_QR = sum(sum(sum(Correct_QR/3)));
P2X_Correct_QR = sum(sum(Correct_QR,3))/3;
P2Y_Correct_QR = sum(sum(Correct_QR,3),2)/3;
P2X_Correct_QR_pos = [find(P2X_Correct_QR,1,'first'),find(P2X_Correct_QR,1,'last')];
P2Y_Correct_QR_pos = [find(P2Y_Correct_QR,1,'first'),find(P2Y_Correct_QR,1,'last')];
Correct_QR_area = diff(P2X_Correct_QR_pos,[],2)*diff(P2Y_Correct_QR_pos,[],2);
Thread_Correct_QR = total_Correct_QR/Correct_QR_area/255;
%% convert to bw again 
Correct_QR_b = im2bw(Correct_QR,Thread_Correct_QR);
temp_C_QR_binary_invert = Correct_QR_b;
CC = bwconncomp(temp_C_QR_binary_invert);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
temp_C_QR_binary_invert(CC.PixelIdxList{idx}) = 0;

imshow(temp_C_QR_binary_invert);

C_QR_Corner = corner(temp_C_QR_binary_invert,'Harris',Harris_corner_detect_num);
hold on
plot(C_QR_Corner(:,1), C_QR_Corner(:,2), 'r*');

hold on;
L = bwlabel(temp_C_QR_binary_invert);
R = regionprops(L,'Area','BoundingBox','PixelList');
NR = numel(R);

maxArea = 0;
index = 1;
cell_R = struct2cell(R);
mat_R_Area = cell2mat(cell_R(1,:));
mat_R_Area = sort(mat_R_Area,'descend');
thread_Area = mat_R_Area(4);
kmax = zeros(1,3);
A = zeros(1,13);

for k = 1:NR
    A(k) = prod(R(k).BoundingBox(3:4));
    if R(k).Area > thread_Area
        maxArea = R(k).Area;
        kmax(index) = k;
        index = index+1;
    end
end

C_Point_outer = zeros(3,8);
C_Point_outer_format2 = zeros(4,6);
C_Point_outer_index = zeros(3,5);

for i=1:size(kmax,2)
BBmax = R(kmax(i)).BoundingBox;
DIAG1 = sum(R(kmax(i)).PixelList,2);
DIAG2 = diff(R(kmax(i)).PixelList,[],2);

[~,dUL] = min(DIAG1);    [~,dDR] = max(DIAG1);
[~,dDL] = min(DIAG2);    [m,dUR] = max(DIAG2);

pts = R(kmax(i)).PixelList([dUL dDL dDR dUR dUL],:);
C_Point_outer(i,:) = [pts(1,:),pts(2,:),pts(3,:),pts(4,:)];
C_Point_outer_format2(:,2*i-1:2*i) = pts(1:4,:);
C_Point_outer_index(i,:) = [kmax(i) dUL dDL dDR dUR];
h_pts = plot(pts(:,1),pts(:,2),'m','linewidth',3);
end

%% detect the size of the gird
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
    C_QR_bound_box(3) = C_grid_num*C_grid_size;
end

%% generate grid
C_grid = zeros(C_grid_size);
C_grid(round(C_grid_size/2),round(C_grid_size/2)) = 1;
C_grid_sampler = repmat(C_grid,C_grid_num);

%total_Correct_QR_b = sum(sum(Correct_QR_b));
%P2X_Correct_QR_b = sum(Correct_QR_b);
%P2Y_Correct_QR_b = sum(Correct_QR_b,2);
%P2X_Correct_QR_b_pos = [find(P2X_Correct_QR_b,1,'first'),find(P2X_Correct_QR_b,1,'last')];
%P2Y_Correct_QR_b_pos = [find(P2Y_Correct_QR_b,1,'first'),find(P2Y_Correct_QR_b,1,'last')];
%Thread2X_Correct_QR_b = total_Correct_QR_b/(diff(P2X_Correct_QR_b_pos,[],2));
%Thread2Y_Correct_QR_b = total_Correct_QR_b/(diff(P2Y_Correct_QR_b_pos,[],2));
%Correct_QR_b_leangth = round(mean(diff(P2X_Correct_QR_b_pos,[],2),diff(P2Y_Correct_QR_b_pos,[],2)));
%F2X_Correct_QR_b_pos = [Center_Correct_QR_b_pos(1)-round(Correct_QR_b_leangth/2),P2X_Correct_QR_b_pos(1)+Correct_QR_b_leangth];
%F2Y_Correct_QR_b_pos = [Center_Correct_QR_b_pos(2)-round(Correct_QR_b_leangth/2),P2Y_Correct_QR_b_pos(1)+Correct_QR_b_leangth];
imshow(Correct_QR_b);
%% extract bit information 
offset = C_QR_bound_box(3)-1;
Croped_Correct_QR_b = Correct_QR_b(C_QR_bound_box(2):C_QR_bound_box(2)+offset,C_QR_bound_box(1):C_QR_bound_box(1)+offset);
imshow(Croped_Correct_QR_b);

S_mat = (Croped_Correct_QR_b+1).*C_grid_sampler;
new_S_mat = S_mat(any(S_mat,2),:);
final_S_mat = new_S_mat(:,any(S_mat))-1;
imshow(1-final_S_mat);

%==========================================================================
%==============================Introduction================================
%==========================================================================
% Subject:
%    Face Recognition using Orthogonal Difference Local Binary Pattern
% Author:
%    Pouya Hosseini <970131010>
%    hosseini.pouya7279@gmail.com
%    University of Yasuj
%    Summer 2022
%==========================================================================

% Close all figures, clear command window & workspace
close all;
clear;
clc;

%==========================================================================
%============================Database Images===============================
%==========================================================================
% Folder of image database
Folder = "Image Datasets\Yale\Database Images\";

% Get list of all jpg files in the folder
Filelist = dir(fullfile(Folder,'*.gif'));

% Get number of images in the folder
Num_of_files = numel(Filelist);

% Store images into a cell
Cell_folder = "Image Datasets\Yale\Database Images\";
Cell_filelist = dir(fullfile(Cell_folder,'*.gif'));
Cell_num_of_files = numel(Cell_filelist);

% Prelocate array & cell
Euclidean_distance = zeros(Num_of_files,1);
Cell_read = cell(Cell_num_of_files,1);
Cell_OD_LBP_1 = cell(Cell_num_of_files,1);
Cell_OD_LBP_2 = cell(Cell_num_of_files,1);
Cell_OD_LBP_3 = cell(Cell_num_of_files,1);
Cell_histogram = cell(Cell_num_of_files,1);
Cell_euclidean_distance = cell(Cell_num_of_files,1);
Cell_query = cell(Cell_num_of_files,1);
Cell_location = cell(Cell_num_of_files,1);

% Read images in  the cell
for k = 1 : Cell_num_of_files
    Cell_read{k} = imread(fullfile(Cell_folder,Cell_filelist(k).name));
end

% Read images in the database
for k = 1 : Num_of_files
	Base_file_name = Filelist(k).name;
	Full_file_name = fullfile(Filelist(k).folder,Base_file_name);
	fprintf('Processing image #%d of %d : "%s"\n', ...
             k, Num_of_files,Base_file_name);
	Base_image = imread(Full_file_name);

	% Convert image to double
    Image = im2double(Base_image);

    % Convert image to grayscale if RGB
    if size(Image,3) == 3
        Image = rgb2gray(Image);
    end

    % Get the dimensions of the image
    [Rows_1,Columns_1] = size(Image);

    % Preallocate array
    OD_LBP_1 = zeros(size(Image),'uint8');
    OD_LBP_2 = zeros(size(Image),'uint8');
    OD_LBP_3 = zeros(size(Image),'uint8');
    
    % Local binary patterning algorithm 
    % Size of window ====> 3*3
    for i_2 = 2 : Rows_1 - 1
        for j_2 = 2 : Columns_1 - 1

            % Find out center pixel
            T_c = Image(i_2,j_2);

            % Pixels around center
		    % Upper left
		    T_1 = Image(i_2-1,j_2-1);
            % Upper middle
		    T_2 = Image(i_2-1,j_2);
            % Upper right
		    T_3 = Image(i_2-1,j_2+1);
            % Middle right
		    T_4 = Image(i_2,j_2+1);
            % Lower Right
		    T_5 = Image(i_2+1,j_2+1);
            % Lower middle
		    T_6 = Image(i_2+1,j_2);
            % Lower left
		    T_7 = Image(i_2+1,j_2-1);
            % Middle left
		    T_8 = Image(i_2,j_2-1);

            % Gray level diffrence
            % First group
            D_1 = [(T_1 - T_2),(T_1 - T_c),(T_1 - T_8)];
            D_2 = [(T_7 - T_6),(T_7 - T_c),(T_7 - T_8)];
            D_3 = [(T_5 - T_4),(T_5 - T_c),(T_5 - T_6)];
            D_4 = [(T_3 - T_2),(T_3 - T_c),(T_3 - T_4)];
            % Second group
            D_5 = [(T_2 - T_3),(T_2 - T_c),(T_2 - T_1)];
            D_6 = [(T_8 - T_7),(T_8 - T_c),(T_8 - T_1)];
            D_7 = [(T_6 - T_5),(T_6 - T_c),(T_6 - T_7)];
            D_8 = [(T_4 - T_3),(T_4 - T_c),(T_4 - T_5)];
            
            % Thresholding
            % First group
            A_1 = (sum(D_1)) / (var(double(D_1)));
            A_2 = (sum(D_2)) / (var(double(D_2)));
            A_3 = (sum(D_3)) / (var(double(D_3)));
            A_4 = (sum(D_4)) / (var(double(D_4)));
            % Second group
            A_5 = (sum(D_5)) / (var(double(D_5)));
            A_6 = (sum(D_6)) / (var(double(D_6)));
            A_7 = (sum(D_7)) / (var(double(D_7)));
            A_8 = (sum(D_8)) / (var(double(D_8)));
            
            % Binary pattern for each orthogonal position
            % First group
            % Upper left
		    BP_1 = D_1 - A_1 >= 0;
            % Upper middle
		    BP_2 = D_2 - A_2 >= 0;
            % Upper right
		    BP_3 = D_3 - A_3 >= 0;
            % Middle right
		    BP_4 = D_4 - A_4 >= 0;
            % Second group
            % Lower Right
		    BP_5 = D_5 - A_5 >= 0;
            % Lower middle
		    BP_6 = D_6 - A_6 >= 0;
            % Lower left
		    BP_7 = D_7 - A_7 >= 0;
            % Middle left
		    BP_8 = D_8 - A_8 >= 0;
            
            % Concatenated results
            CC_1 = [BP_1,BP_2,BP_3,BP_4];
            CC_2 = [BP_5,BP_6,BP_7,BP_8];
            CC_3 = [CC_1,CC_2];
            
            % Split results
            BP_9 = (CC_3(1:8));
            BP_10 = (CC_3(9:16));
            BP_11 = (CC_3(17:24));

            % Decimal codes
            % Prelocate array
            LBP_1 = zeros(1);
            LBP_2 = zeros(1);
            LBP_3 = zeros(1);

            for n=1:8
                LBP_1(n) = (2^(n-1)) * (BP_9(1,n));
                LBP_2(n) = (2^(n-1)) * (BP_10(1,n));
                LBP_3(n) = (2^(n-1)) * (BP_11(1,n));
            end

            LBP_1(k) = sum(LBP_1);
            LBP_2(k) = sum(LBP_2);
            LBP_3(k) = sum(LBP_3);

            Database_OD_LBP = [LBP_1,LBP_2,LBP_3];

            % OD-LBP image
            OD_LBP_1(i_2,j_2) = LBP_1(k);
            OD_LBP_2(i_2,j_2) = LBP_2(k);
            OD_LBP_3(i_2,j_2) = LBP_3(k);

            Cell_OD_LBP_1{k} = OD_LBP_1;
            Cell_OD_LBP_2{k} = OD_LBP_2;
            Cell_OD_LBP_3{k} = OD_LBP_3;
        end
    end

%--------------------------------------------------------------------------
% Divide database LBP 1
% Figure out where to divide it
[Rows_2,Columns_2] = size(OD_LBP_1);
Rows_3 = round(Rows_2 / 3);
Column_3 = round(Columns_2 / 3);

% Divide image to 9 blocks (3*3)
% Upper left
Block_1_1 = OD_LBP_1(1:Rows_3,1:Column_3);

% Upper middle
Block_2_1 = OD_LBP_1(1:Rows_3,Column_3+1:2*Column_3);

% Upper right
Block_3_1 = OD_LBP_1(1:Rows_3,2*Column_3+1:end);

% Middle left
Block_4_1 = OD_LBP_1(Rows_3+1:2*Rows_3,1:Column_3);

% Center
Block_5_1 = OD_LBP_1(Rows_3+1:2*Rows_3,Column_3+1:2*Column_3);

% Middle right
Block_6_1 = OD_LBP_1(Rows_3+1:2*Rows_3,2*Column_3+1:end);

% Lower left
Block_7_1 = OD_LBP_1(2*Rows_3+1:end,1:Column_3);

% Lower middle
Block_8_1 = OD_LBP_1(2*Rows_3+1:end,Column_3+1:2*Column_3);

% Lower right
Block_9_1 = OD_LBP_1(2*Rows_3+1:end,2*Column_3+1:end);

% Histogram of each block
Histogram_1_1 = imhist(uint8(Block_1_1));
Histogram_2_1 = imhist(uint8(Block_2_1));
Histogram_3_1 = imhist(uint8(Block_3_1));
Histogram_4_1 = imhist(uint8(Block_4_1));
Histogram_5_1 = imhist(uint8(Block_5_1));
Histogram_6_1 = imhist(uint8(Block_6_1));
Histogram_7_1 = imhist(uint8(Block_7_1));
Histogram_8_1 = imhist(uint8(Block_8_1));
Histogram_9_1 = imhist(uint8(Block_9_1));

% Total histogram of LBP 1
Total_histogram_1 = [Histogram_1_1(:);Histogram_2_1(:);Histogram_3_1(:);
                     Histogram_4_1(:);Histogram_5_1(:);Histogram_6_1(:);
                     Histogram_7_1(:);Histogram_8_1(:);Histogram_9_1(:)];

%--------------------------------------------------------------------------
% Divide database LBP 2
% Figure out where to divide it
[Rows_3,Columns_3] = size(OD_LBP_2);
Rows_4 = round(Rows_3 / 3);
Column_4 = round(Columns_3 / 3);

% Divide image to 9 blocks (3*3)
% Upper left
Block_1_2 = OD_LBP_2(1:Rows_4,1:Column_4);

% Upper middle
Block_2_2 = OD_LBP_2(1:Rows_4,Column_4+1:2*Column_4);

% Upper right
Block_3_2 = OD_LBP_2(1:Rows_4,2*Column_4+1:end);

% Middle left
Block_4_2 = OD_LBP_2(Rows_4+1:2*Rows_4,1:Column_4);

% Center
Block_5_2 = OD_LBP_2(Rows_4+1:2*Rows_4,Column_4+1:2*Column_4);

% Middle right
Block_6_2 = OD_LBP_2(Rows_4+1:2*Rows_4,2*Column_4+1:end);

% Lower left
Block_7_2 = OD_LBP_2(2*Rows_4+1:end,1:Column_4);

% Lower middle
Block_8_2 = OD_LBP_2(2*Rows_4+1:end,Column_4+1:2*Column_4);

% Lower right
Block_9_2 = OD_LBP_2(2*Rows_4+1:end,2*Column_4+1:end);

% Histogram of each block
Histogram_1_2 = imhist(uint8(Block_1_2));
Histogram_2_2 = imhist(uint8(Block_2_2));
Histogram_3_2 = imhist(uint8(Block_3_2));
Histogram_4_2 = imhist(uint8(Block_4_2));
Histogram_5_2 = imhist(uint8(Block_5_2));
Histogram_6_2 = imhist(uint8(Block_6_2));
Histogram_7_2 = imhist(uint8(Block_7_2));
Histogram_8_2 = imhist(uint8(Block_8_2));
Histogram_9_2 = imhist(uint8(Block_9_2));

% Total histogram of LBP 2
Total_histogram_2 = [Histogram_1_2(:);Histogram_2_2(:);Histogram_3_2(:);
                     Histogram_4_2(:);Histogram_5_2(:);Histogram_6_2(:);
                     Histogram_7_2(:);Histogram_8_2(:);Histogram_9_2(:)];

%--------------------------------------------------------------------------
% Divide database LBP 3
% Figure out where to divide it
[Rows_5,Columns_5] = size(OD_LBP_3);
Rows_6 = round(Rows_5 / 3);
Column_6 = round(Columns_5 / 3);

% Divide image to 9 blocks (3*3)
% Upper left
Block_1_3 = OD_LBP_3(1:Rows_6,1:Column_6);

% Upper middle
Block_2_3 = OD_LBP_3(1:Rows_6,Column_6+1:2*Column_6);

% Upper right
Block_3_3 = OD_LBP_3(1:Rows_6,2*Column_6+1:end);

% Middle left
Block_4_3 = OD_LBP_3(Rows_6+1:2*Rows_6,1:Column_6);

% Center
Block_5_3 = OD_LBP_3(Rows_6+1:2*Rows_6,Column_6+1:2*Column_6);

% Middle right
Block_6_3 = OD_LBP_3(Rows_6+1:2*Rows_6,2*Column_6+1:end);

% Lower left
Block_7_3 = OD_LBP_3(2*Rows_6+1:end,1:Column_6);

% Lower middle
Block_8_3 = OD_LBP_3(2*Rows_6+1:end,Column_6+1:2*Column_6);

% Lower right
Block_9_3 = OD_LBP_3(2*Rows_6+1:end,2*Column_6+1:end);

% Histogram of each block
Histogram_1_3 = imhist(uint8(Block_1_3));
Histogram_2_3 = imhist(uint8(Block_2_3));
Histogram_3_3 = imhist(uint8(Block_3_3));
Histogram_4_3 = imhist(uint8(Block_4_3));
Histogram_5_3 = imhist(uint8(Block_5_3));
Histogram_6_3 = imhist(uint8(Block_6_3));
Histogram_7_3 = imhist(uint8(Block_7_3));
Histogram_8_3 = imhist(uint8(Block_8_3));
Histogram_9_3 = imhist(uint8(Block_9_3));

% Total histogram of LBP 3
Total_histogram_3 = [Histogram_1_3(:);Histogram_2_3(:);Histogram_3_3(:);
                     Histogram_4_3(:);Histogram_5_3(:);Histogram_6_3(:);
                     Histogram_7_3(:);Histogram_8_3(:);Histogram_9_3(:)];

% Total histogram of LBP
Total_histogram = [Total_histogram_1(:);Total_histogram_2(:);
                   Total_histogram_3(:)];

Cell_histogram{k} = Total_histogram';

% Concatenate arrays within a loop
All_total_histogram = vertcat(Cell_histogram{:});
save('OD_LBP_Features','All_total_histogram')

% Show images (if you want to show images uncomment it!)
% figure("Name",['Image Number ',num2str(k)]);
% set(gcf,'position',[365,110,800,600]);
% 
% subplot(4,3,(1:3));
% imshow(Base_image)
% title(['Image Number ',num2str(k)]);
% 
% subplot(4,3,4);
% imshow(OD_LBP_1)
% title('LBP 1 Image');
% 
% subplot(4,3,5);
% imshow(OD_LBP_2)
% title('LBP 2 Image');
% 
% subplot(4,3,6);
% imshow(OD_LBP_3)
% title('LBP 3 Image');
% 
% subplot(4,3,7);
% bar(Total_histogram_1,'BarWidth',1,'EdgeColor','none');
% grid on;
% title('LBP 1 Histogram');
% 
% subplot(4,3,8);
% bar(Total_histogram_2,'BarWidth',1,'EdgeColor','none');
% grid on;
% title('LBP 2 Histogram');
% 
% subplot(4,3,9);
% bar(Total_histogram_3,'BarWidth',1,'EdgeColor','none');
% grid on;
% title('LBP 3 Histogram');
% 
% subplot(4,3,(10:12));
% bar(Total_histogram,'BarWidth',1,'EdgeColor','none');
% grid on;
% title('Total LBP Histogram');

end

% Number of Top images
fprintf('\n')
Num_of_top = input('>>Please enter the number of top images: ');

% Calculate Precision
OD_LBP_Precision(Num_of_files,Num_of_top)
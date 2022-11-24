function OD_LBP_Precision(Num_of_files,Num_of_top)

% Prelocate array & cell
Euclidean_distance = zeros(Num_of_files,1);
Cell_euclidean_distance = cell(Num_of_files,1);
Cell_query = cell(Num_of_files,1);
Cell_location = cell(Num_of_files,1);
Pre = zeros(1,Num_of_top);

for Num_of_top = 1:Num_of_top

load('OD_LBP_Features.mat','All_total_histogram')

% Prelocate array & cell
Class = zeros(1,Num_of_files);
Cell_counter = zeros(1,Num_of_files);
Percentage_of_each_image = zeros(1,Num_of_files);
Location_of_top = cell(1,Num_of_files);

% Making a class for comparison
for c = 1:15
    Class(c) = c;
end

Class = repelem(Class,11);

% Find euclidean distance & accuracy
for m = 1 : Num_of_files
    Query = All_total_histogram(m,:);
    Cell_query{m,1} = Query;

    for n = 1 : Num_of_files
        Euclidean_distance(n,1) = sqrt(sum((Query-...
                                            All_total_histogram(n,:)).^2));
    end

    Cell_euclidean_distance{m,1} = Euclidean_distance;
    [~,Location] = sort(Euclidean_distance);
    Cell_location{m,1} = Location;

% Find & show matching images (if you want to show images uncomment it!)
%      for i = 1 : Num_of_top
%          figure("Name",['Matching Database Image Number ',num2str(i)]);
%          set(gcf,'position',[600,300,800,600]);
%          imshow(Cell{Location(i,1)})
%          title(['Matching Database Image Number ',num2str(i)]);
%      end

%==========================================================================
%========================Accuracy of Algorithm=============================
%==========================================================================
     % Put location of top images in a cell
     Location_of_top{m}=Cell_location{m,1}(1:Num_of_top)';
     
     % Set a counter for accuracy
     Counter = 0;

     for z = 1:Num_of_top

         if Class(Location_of_top{1,m}(z)) == Class(m)
             Counter = Counter + 1;
         end

         Cell_counter(1,m) = Counter;
     end

     Percentage_of_each_image(m) = (Counter/Num_of_top)*100;
end

Total_accuracy_percentage = mean(Percentage_of_each_image);

fprintf('Total Percentage of %d images is ===> %0.2f \n',Num_of_top, ...
         Total_accuracy_percentage);

Pre(Num_of_top) = Total_accuracy_percentage;
end
end
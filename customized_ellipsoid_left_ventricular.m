clear;
% Parametric equations for the ellipsoid
theta_ellipsoid = linspace(0, 2*pi, 100); % 经度
phi_ellipsoid = linspace(pi/2, pi, 100); % 纬度
latitude = -71:1:-2; % 纬度
latitude_num = length(latitude);
longitude_num = length(theta_ellipsoid);
fiber_turns = 1;

[theta_ellipsoid, phi_ellipsoid] = meshgrid(theta_ellipsoid, phi_ellipsoid);

a_base = 38; % Semi-minor axis length of the ellipsoid (shorter)
a_base_in =31; % Semi-minor axis of left ventricular inner surface
b_base = 38;  % Semi-minor axis length of the ellipsoid (shorter)
c_base = 76;  % Semi-major axis length of the ellipsoid
c_base_in = 69; % Semi-major axis of left ventricular inner surface
spiral_angle = pi/40:pi/40:2*pi; % this is to rotate spiral ralative to z axis
spiral_num = length(spiral_angle);
radius_ratio = 0.8250:0.0125:1; % the radius of inner layer is 31mm
layer_num = length(radius_ratio);
x_spiral_r = cell(layer_num,1);
y_spiral_r = cell(layer_num,1);
z_spiral_r = cell(layer_num,1);
x_spiral_key = cell(layer_num,1);
y_spiral_key = cell(layer_num,1);
z_spiral_key = cell(layer_num,1);
zero_angle_ind = zeros(layer_num, latitude_num);

figure;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');

hold on;

for layer_ind=1:layer_num
    ratio = radius_ratio(layer_ind);
    a = a_base*ratio;
    b = b_base*ratio;
    c = c_base + (ratio-radius_ratio(end))*(c_base_in-c_base)/(radius_ratio(1)-radius_ratio(end));
    % c = c_base*(ratio/(a_base_in/a_base)*(c_base_in/c_base));
    x_ellipsoid = a * cos(theta_ellipsoid) .* sin(phi_ellipsoid);
    y_ellipsoid = b * sin(theta_ellipsoid) .* sin(phi_ellipsoid);
    z_ellipsoid = c * cos(phi_ellipsoid);
    surf(x_ellipsoid, y_ellipsoid, z_ellipsoid,'FaceAlpha',0.1,'EdgeColor','none');

    % for i=1:5:100
    %     for j=1:5:100
    %         plot3(x_ellipsoid(i,j),y_ellipsoid(i,j),z_ellipsoid(i,j),'.','Color','g','MarkerSize',10);
    %     end
    % end
    % Generate the spiral points
    num_points = 1000;
    % if (ratio - 0.9375)<1e-5
    %     theta = linspace(0, fiber_turns*pi, num_points);
    % elseif (ratio - 0.95)<1e-5
    %     theta = linspace(0, fiber_turns*2*pi, num_points); % 经度，a larger up limit result in more turns
    % elseif (ratio - 0.9625)<1e-5
    %     theta = linspace(0, fiber_turns*3*pi, num_points);
    % elseif (ratio - 0.975)<1e-5
    %     theta = linspace(0, -fiber_turns*3*pi, num_points);
    % elseif (ratio - 0.9875)<1e-5
    %     theta = linspace(0, -fiber_turns*2*pi, num_points);
    % elseif (ratio - 1)<1e-5
    %     theta = linspace(0, -fiber_turns*pi, num_points);
    % end
    max_ang = -pi + (ratio-radius_ratio(end))*2*pi ...
        /(radius_ratio(1)-radius_ratio(end));
    % disp(max_ang/pi);
    theta = linspace(0, max_ang, num_points);

    radius = linspace(0, a, num_points); % Radius increases linearly from 0 to 10
    x_spiral = radius .* cos(theta);
    y_spiral = radius .* sin(theta);

    % Apply rotation to the spiral points
    x_spiral_key{layer_ind} = zeros(spiral_num,latitude_num); % x_spiral_key{j}的代表第j层的纤维，每一行代表一条纤维的x坐标,由latitude_num个点组成一根纤维
    y_spiral_key{layer_ind} = zeros(spiral_num,latitude_num);
    z_spiral_key{layer_ind} = zeros(spiral_num,latitude_num);
    x_spiral_r{layer_ind} = zeros(spiral_num,num_points); % x_spiral_r{j}的代表第j层的纤维，每一行代表一条纤维的x坐标
    y_spiral_r{layer_ind} = zeros(spiral_num,num_points); % y_spiral_r{j}每一行代表一条纤维的y坐标
    z_spiral_r{layer_ind} = zeros(spiral_num,num_points); % z_spiral_r{j}每一行代表一条纤维的z坐标
    for i=1:spiral_num
        x_spiral_r{layer_ind}(i,:) = x_spiral * cos(spiral_angle(i)) - y_spiral * sin(spiral_angle(i));
        y_spiral_r{layer_ind}(i,:) = x_spiral * sin(spiral_angle(i)) + y_spiral * cos(spiral_angle(i));
        z_spiral_r{layer_ind}(i,:) = -real(sqrt((1-x_spiral.^2/a^2-y_spiral.^2/b^2)*c^2));
        % Plot the spiral
        plot3(x_spiral_r{layer_ind}(i,:)', y_spiral_r{layer_ind}(i,:)', ...
            z_spiral_r{layer_ind}(i,:)', 'Color', [0.8, 0.2, 0.2, 0.3],'LineWidth', 1);
    end

    for i=1:latitude_num
        height = latitude(i);
        ind = find( abs(z_spiral_r{layer_ind}(1,:)-height) == ...
            min(abs(z_spiral_r{layer_ind}(1,:)-height)) );
        % 将不同层连接的点的索引记录下来
        x_spiral_key{layer_ind}(:,i) = x_spiral_r{layer_ind}(:,ind); 
        y_spiral_key{layer_ind}(:,i) = y_spiral_r{layer_ind}(:,ind); 
        z_spiral_key{layer_ind}(:,i) = z_spiral_r{layer_ind}(:,ind); 

        point_angle = atan2(y_spiral_key{layer_ind}(:,i),x_spiral_key{layer_ind}(:,i));
        [~, zero_angle_ind(layer_ind,i)] = min(abs(point_angle));
        plot3(x_spiral_key{layer_ind}(:,i), ...
            y_spiral_key{layer_ind}(:,i), ...
            z_spiral_key{layer_ind}(:,i),'.','MarkerSize',5);
        plot3(x_spiral_key{layer_ind}(zero_angle_ind(layer_ind,i),i), ...
            y_spiral_key{layer_ind}(zero_angle_ind(layer_ind,i),i), ...
            z_spiral_key{layer_ind}(zero_angle_ind(layer_ind,i),i),'.','MarkerSize',20,'Color','r');
    end
end
view([0,0]);
clear x_spiral_r  y_spiral_r  z_spiral_r;
% figure;
% plot(x_spiral, y_spiral, 'r', 'LineWidth', 2);
%%
fileID = fopen('lines.vtk', 'w');

% Write VTK header
fprintf(fileID, '# vtk DataFile Version 3.0\n');
fprintf(fileID, 'Lines example\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET POLYDATA\n');
numLines = layer_num*spiral_num;
numGaps = latitude_num*spiral_num;
numPointsPerLine = latitude_num;
lines = zeros(numPointsPerLine, 3, numLines);
gaps = zeros(layer_num, 3, numGaps);

figure;
gm = importGeometry("ventricular1.step");
pdegplot(gm,"FaceAlpha",0.3);
hold on;
i = 1;
for layer_ind=1:layer_num
    for spiral_i=1:length(spiral_angle)
        lines(:,1,i) = x_spiral_key{layer_ind}(spiral_i,:)';
        lines(:,2,i) = y_spiral_key{layer_ind}(spiral_i,:)';
        lines(:,3,i) = z_spiral_key{layer_ind}(spiral_i,:)';
        plot3(lines(:,1,i), lines(:,2,i), ...
            lines(:,3,i), 'Color', [0.8, 0.2, 0.2, 0.5],'LineWidth', 1);
        i = i+1;
    end
end
i = 1;
for spiral_ind = 1:spiral_num
    % longitude_interest = 2*pi/50;
    % longitude_ind = round(longitude_interest/theta_ellipsoid(1,2));
    for latitude_i=1:latitude_num
        for layer_ind=1:layer_num
            if zero_angle_ind(layer_ind,latitude_i)+spiral_ind <= spiral_num
                angle_ind = zero_angle_ind(layer_ind,latitude_i)+spiral_ind;
            else
                angle_ind = zero_angle_ind(layer_ind,latitude_i)+spiral_ind-spiral_num;
            end
            gaps(layer_ind,:,i) = [x_spiral_key{layer_ind}(angle_ind,latitude_i), ...
                y_spiral_key{layer_ind}(angle_ind,latitude_i) ...
                z_spiral_key{layer_ind}(angle_ind,latitude_i)];
            % plot3(x_spiral_key{layer_ind}(angle_ind,i), ...
            %     y_spiral_key{layer_ind}(angle_ind,i), ...
            %     z_spiral_key{layer_ind}(angle_ind,i),'.','MarkerSize',5,'Color','b');
        end
        plot3(gaps(:,1,i), gaps(:,2,i), gaps(:,3,i), 'Color', [0.2, 0.2, 0.8, 0.5],'LineWidth', 1);
        i = i+1;
    end
end
axis equal;
axis off;
view([20,18]);

%% Open a file for writing
fileID = fopen('lines.vtk', 'w');

% Write VTK header
fprintf(fileID, '# vtk DataFile Version 3.0\n');
fprintf(fileID, 'Lines example\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET POLYDATA\n');

% Write the points
fprintf(fileID, 'POINTS %d float\n', numLines * numPointsPerLine);
for i = 1:numLines
    for j = 1:numPointsPerLine
        fprintf(fileID, '%f %f %f\n', lines(j, :, i));
    end
end

% Write the fiber lines
fprintf(fileID, 'LINES %d %d\n', numLines, numLines * (numPointsPerLine + 1));
for i = 1:numLines
    fprintf(fileID, '%d ', numPointsPerLine);
    for j = 1:numPointsPerLine
        fprintf(fileID, '%d ', (i - 1) * numPointsPerLine + j - 1);
    end
    fprintf(fileID, '\n');
end

% Close the file
fclose(fileID);
%%
fileID = fopen('gaps.vtk', 'w');

% Write VTK header
fprintf(fileID, '# vtk DataFile Version 3.0\n');
fprintf(fileID, 'gaps example\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET POLYDATA\n');

% Write the points
fprintf(fileID, 'POINTS %d float\n', numGaps*layer_num);
for i = 1:numGaps
    for j = 1:layer_num
        fprintf(fileID, '%f %f %f\n', gaps(j, :, i));
    end
end

% Write the fiber lines
fprintf(fileID, 'LINES %d %d\n', numGaps, numGaps * (layer_num + 1));
for i = 1:numGaps
    fprintf(fileID, '%d ', layer_num);
    for j = 1:layer_num
        fprintf(fileID, '%d ', (i - 1) * layer_num + j - 1);
    end
    fprintf(fileID, '\n');
end

% Close the file
fclose(fileID);
%% interpolation
choice = 1;
if choice == 1
    addpath('.\func');
    x_spiral_interp = cell(layer_num,1);
    y_spiral_interp = cell(layer_num,1);
    z_spiral_interp = cell(layer_num,1);
    tic;
    parfor layer_ind=1:layer_num
        % x_spiral_interp{layer_ind} = zeros(spiral_num,latitude_num); % x_spiral_key{j}的代表第j层的纤维，每一行代表一条纤维的x坐标,由latitude_num个点组成一根纤维
        % y_spiral_interp{layer_ind} = zeros(spiral_num,latitude_num);
        % z_spiral_interp{layer_ind} = zeros(spiral_num,latitude_num);
        for spiral_ind = 1:spiral_num
            curve_data = [x_spiral_key{layer_ind}(spiral_ind,:)', ...
                y_spiral_key{layer_ind}(spiral_ind,:)', ...
                z_spiral_key{layer_ind}(spiral_ind,:)'];
            [z, index, ~] = unique(curve_data(:, 3), 'stable');  % z coordinate has unique value
            % 提取数据
            filteredData = curve_data(index, :);
            x = filteredData(:, 1);
            y = filteredData(:, 2);
            % 设定间隔长度
            interval = 0.1; % cell length=100um, gap length can be neglected
            % 这里生成的数据是先按原数据拟合样条曲线，然后等间隔取点
            totalLength = arclength(x, y, z, 's');
            n = round(totalLength / interval);
            tValues = linspace(0, 1, n);
            intep_points = interparc(tValues, x, y, z, 's');
            x_spiral_interp{layer_ind} = [x_spiral_interp{layer_ind} ;intep_points(:,1)'];
            y_spiral_interp{layer_ind} = [y_spiral_interp{layer_ind} ;intep_points(:,2)'];
            z_spiral_interp{layer_ind} = [z_spiral_interp{layer_ind} ;intep_points(:,3)'];
        end
    end
    toc;
    % 绘制拟合曲线和取点结果
    figure; axis equal; hold on;
    plot3(x_spiral_key{1}(1,:)', ...
                y_spiral_key{1}(1,:)', ...
                z_spiral_key{1}(1,:)', '-r', 'LineWidth', 1); % 原曲线
    scatter3(x_spiral_interp{1}(1,:)', ...
                y_spiral_interp{1}(1,:)', ...
                z_spiral_interp{1}(1,:)'); % 取点
    xlabel('X / mm'); ylabel('Y / mm'); zlabel('Z / mm');
    legend('原曲线','取点结果');
    
    save("customized_ventricular_fiber.mat","x_spiral_interp","y_spiral_interp", ...
        "z_spiral_interp","x_spiral_key","y_spiral_key","z_spiral_key","zero_angle_ind");
end

function maketable(pdfA,pdfV,pdfVabs, pdfAabs,folderout)
%Create some example variables and a table
Name = {'Mean', 'STD' ,'Skewness','Kurtosis'};

Vabs = [pdfVabs.mean, pdfVabs.std ,pdfVabs.skewness, pdfVabs.flatness];
Aabs = [pdfAabs.mean, pdfAabs.std, pdfAabs.skewness, pdfAabs.flatness];

vx = [pdfV(1).mean, pdfV(1).std ,pdfV(1).skewness, pdfV(1).flatness];
vy = [pdfV(2).mean, pdfV(2).std, pdfV(2).skewness, pdfV(2).flatness];
vz = [pdfV(3).mean, pdfV(3).std, pdfV(3).skewness, pdfV(3).flatness];

ax = [pdfA(1).mean, pdfA(1).std, pdfA(1).skewness, pdfA(1).flatness];
ay = [pdfA(2).mean, pdfA(2).std, pdfA(2).skewness, pdfA(2).flatness];
az = [pdfA(3).mean, pdfA(3).std, pdfA(3).skewness, pdfA(3).flatness];


% Create a new table with consistent data types and handling of missing values
T = table(Name', Vabs', Aabs' ,vx', vy', vz', ax', ay', az', 'VariableNames', {'Name', 'Vabs','Aabs','Vx','Vy','Vz','Ax','Ay','Az'});

% Create a figure and display only the table without any axes or labels
f = figure('Visible', 'on');
uitable(f, 'Data', table2cell(T), 'ColumnName', T.Properties.VariableNames, ...
    'Position', [000, 100, 1000, 100], 'FontSize', 12);

print(f,[folderout filesep 'stats.png'], '-dpng', '-r300'); % Adjust the DPI as needed

% Close the figure
%close(f);
end
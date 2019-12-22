% ZIPproject - create a zip file of the project

zipfilename = 'HiT_demo_20100820.zip';

%% specify files to put in ZIP
filelist = {};
count = 1;

% readme
filelist{count} = sprintf('README.txt');
disp([num2str(count) ': added ' filelist{count}]);
count = count + 1;

% .m files in current folder
names = dir('*.m');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% .c files in current folder
names = dir('*.c');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% .cpp files in current folder
names = dir('*.cpp');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% .h files in current folder
names = dir('*.h');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% images
paths = dir('./*');
for j = 1:length(paths)
	if paths(j).isdir && paths(j).name(1) ~= '.'
		folder = paths(j).name;
		names = dir([folder '/*']);
		for i = 1:length(names)
			if names(i).name(1) ~= '.'
				filelist{count} = sprintf('%s/%s',folder, names(i).name);
				disp([num2str(count) ': added ' filelist{count}]);
				count = count + 1;
			end
		end
	end
end


names = dir('*.mat');
for i = 1:length(names)
    if names(i).name(1) ~= '.'
        filelist{count} = sprintf('%s',names(i).name);
        disp([num2str(count) ': added ' filelist{count}]);
        count = count + 1;
    end
end

%% zip
zip(zipfilename,filelist,pwd);

